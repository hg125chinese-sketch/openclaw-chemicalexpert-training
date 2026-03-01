---
name: chem-scaffold-conditioned-gen
description: Generate molecules conditioned on target-relevant scaffolds and fragments, with built-in KL-collapse prevention. Covers scaffold-constrained VAE, fragment-conditioned generation, cyclical annealing, and hinge binder coverage KPIs. Use when unconditional generation produces candidates outside the target's chemical space.
homepage: https://arxiv.org/abs/1903.10145
metadata: { "openclaw": { "emoji": "🧬", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "numpy", "pandas", "selfies", "sklearn"] } } }
---

# Scaffold-Conditioned Generation with Anti-Collapse

Unconditional molecular generation (vanilla VAE on SELFIES/SMILES) explores broad chemical space but has no mechanism to stay near the target's pharmacophore. For kinase inhibitors, this means the generator happily produces molecules without hinge binder motifs — the one structural feature that makes a kinase inhibitor a kinase inhibitor. This skill fixes generation at the source instead of relying on post-hoc filtering.

## When to Use

- When hinge binder coverage ratio (skill 13) is CRITICAL (<0.5 of training set frequency)
- When unconditional VAE produces <5% kinase-like molecules in the candidate pool
- When KL/dim < 0.1 in a trained VAE (KL collapse confirmed)
- When starting Cycle 2+ of a DMTA loop where Cycle 1 showed generation quality problems
- When scaffold diversity in the generated pool is significantly lower than the training set

## Core Philosophy

1. **Fix generation, not filtering.** If 97% of your generated molecules lack kinase motifs, no amount of post-hoc filtering will save you. You need molecules that are born with the right chemistry.
2. **Condition on fragments, not whole molecules.** Conditioning on a complete scaffold is too restrictive (you just get the training set back). Conditioning on a hinge-binding fragment (e.g., aminopyrimidine) lets the generator explore diverse elaborations while anchoring to the pharmacophore.
3. **KL collapse is the default failure mode of molecular VAEs.** Assume it will happen and build prevention into every training run. Cyclical annealing is the minimum; free bits and word dropout are additional insurance.
4. **Coverage KPIs are hard gates, not suggestions.** If the generated pool fails the hinge binder coverage gate (skill 13), the model is broken. Don't proceed to docking.
5. **Every generation experiment is a controlled experiment.** Change one thing at a time, measure with the same KPIs, and commit results. No "I think this looks better."

## Phase 1: Diagnosing the Generation Problem

### 1.1 Confirm KL Collapse

Before conditioning, confirm whether the current VAE has collapsed:

```python
#!/opt/conda/envs/chem/bin/python
"""Diagnose VAE health: KL collapse, reconstruction quality, latent usage."""
import torch
import numpy as np


def diagnose_vae(model, dataloader, device="cuda"):
    """Run diagnostic metrics on a trained VAE.

    Returns:
        dict with kl_per_dim, active_units, au_ratio, recon_accuracy, diagnosis
    """
    model.eval()
    all_mu, all_logvar = [], []
    n_correct, n_total = 0, 0

    with torch.no_grad():
        for batch in dataloader:
            batch = batch.to(device)
            mu, logvar = model.encode(batch)
            recon = model.decode(model.reparameterize(mu, logvar))
            all_mu.append(mu.cpu())
            all_logvar.append(logvar.cpu())
            n_correct += (recon.argmax(dim=-1) == batch).all(dim=-1).sum().item()
            n_total += batch.size(0)

    mu_cat = torch.cat(all_mu, dim=0)
    logvar_cat = torch.cat(all_logvar, dim=0)

    # Per-dimension KL: 0.5 * (mu^2 + exp(logvar) - logvar - 1)
    kl_per_dim = 0.5 * (mu_cat.pow(2) + logvar_cat.exp() - logvar_cat - 1).mean(dim=0)
    kl_total = kl_per_dim.sum().item()
    active_units = (kl_per_dim > 0.01).sum().item()
    n_dims = kl_per_dim.shape[0]
    au_ratio = active_units / n_dims
    recon_acc = n_correct / n_total if n_total > 0 else 0

    if au_ratio > 0.5 and kl_total / n_dims > 0.1:
        diagnosis = "healthy"
    elif au_ratio > 0.1:
        diagnosis = "partial_collapse"
    else:
        diagnosis = "full_collapse"

    return {
        "kl_per_dim": kl_per_dim.numpy(),
        "kl_total": kl_total,
        "kl_per_dim_mean": kl_total / n_dims,
        "active_units": active_units,
        "total_dims": n_dims,
        "au_ratio": au_ratio,
        "recon_accuracy": recon_acc,
        "diagnosis": diagnosis,
    }
```

### 1.2 Decision Based on Diagnosis

```
diagnosis == "full_collapse" (AU ratio < 0.1)?
+-- YES -> Full retraining with cyclical annealing (Phase 2)
|         The latent space is useless; patching won't help.
|
+-- diagnosis == "partial_collapse" (AU 0.1-0.5)?
|   +-- Try fine-tuning with cyclical annealing first
|       If still failing after 2 attempts -> full retrain
|
+-- diagnosis == "healthy" but hinge coverage still low?
    +-- The VAE works but training data doesn't have enough
        kinase motifs. Fix the training set (Phase 3) and/or
        add scaffold conditioning (Phase 4).
```

## Phase 2: Anti-Collapse Training Strategies

### 2.1 Cyclical Annealing (Primary Defense)

This is the single most effective fix for KL collapse in molecular VAEs.

```python
class CyclicalAnnealer:
    """Cyclical KL annealing schedule.

    Reference: Fu et al., "Cyclical Annealing Schedule: A Simple Approach
    to Mitigating KL Vanishing" (NAACL 2019).

    The key insight: periodically resetting beta to 0 forces the model to
    re-learn the latent space from a better starting point each cycle.
    """

    def __init__(self, total_steps, n_cycles=4, ratio=0.5, max_beta=1.0):
        self.total_steps = total_steps
        self.n_cycles = n_cycles
        self.ratio = ratio
        self.max_beta = max_beta

    def __call__(self, step):
        cycle_length = self.total_steps / self.n_cycles
        tau = (step % cycle_length) / cycle_length
        if tau <= self.ratio:
            beta = self.max_beta * (tau / self.ratio)
        else:
            beta = self.max_beta
        return beta


def train_vae_with_cyclical_annealing(model, train_loader, val_loader,
                                       n_epochs=100, n_cycles=4,
                                       lr=1e-3, device="cuda",
                                       free_bits=0.0, word_dropout=0.0):
    """Train a VAE with cyclical annealing and optional free bits + word dropout.

    Returns dict with training history and final diagnostics.
    """
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    total_steps = n_epochs * len(train_loader)
    annealer = CyclicalAnnealer(total_steps, n_cycles=n_cycles)

    history = {"epoch": [], "train_loss": [], "recon_loss": [],
               "kl_loss": [], "beta": [], "val_loss": []}
    global_step = 0

    for epoch in range(n_epochs):
        model.train()
        epoch_recon, epoch_kl, n_batches = 0, 0, 0

        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()

            # Word dropout on decoder input
            if word_dropout > 0:
                decoder_input = apply_word_dropout(batch, word_dropout)
            else:
                decoder_input = batch

            recon, mu, logvar = model(batch)  # encode full, decode may use dropout
            recon_loss = model.reconstruction_loss(recon, batch)

            # KL with optional free bits
            kl_per_dim = 0.5 * (mu.pow(2) + logvar.exp() - logvar - 1)
            if free_bits > 0:
                kl_per_dim = torch.clamp(kl_per_dim, min=free_bits)
            kl_loss = kl_per_dim.sum(dim=-1).mean()

            beta = annealer(global_step)
            loss = recon_loss + beta * kl_loss
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_recon += recon_loss.item()
            epoch_kl += kl_loss.item()
            n_batches += 1
            global_step += 1

        # Log
        history["epoch"].append(epoch)
        history["recon_loss"].append(epoch_recon / n_batches)
        history["kl_loss"].append(epoch_kl / n_batches)
        history["beta"].append(beta)

        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{n_epochs} | "
                  f"Recon: {epoch_recon/n_batches:.4f} | "
                  f"KL: {epoch_kl/n_batches:.4f} | Beta: {beta:.4f}")

    final_diag = diagnose_vae(model, val_loader, device)
    history["final_diagnostics"] = final_diag
    return history


def apply_word_dropout(input_seq, dropout_rate=0.1, pad_idx=0):
    """Randomly replace input tokens with pad to force latent usage."""
    if dropout_rate == 0:
        return input_seq
    mask = torch.bernoulli(
        torch.full_like(input_seq, 1.0 - dropout_rate, dtype=torch.float)
    ).long()
    return input_seq * mask + pad_idx * (1 - mask)
```

### 2.2 Anti-Collapse Experiment Template

Every experiment follows this controlled template:

```python
EXPERIMENT_TEMPLATE = {
    "experiment_id": "ipf_cycle2_gen_exp01",
    "target": "ALK5/TGFBR1",
    "baseline": "Cycle 1 vanilla SELFIES VAE (KL/dim=0.031, AU=3/128)",
    "changes": {
        "annealing": "cyclical, 4 cycles",
        "free_bits": 0.25,
        "word_dropout": 0.1,
        # Change ONLY ONE per experiment for clean comparison
    },
    "metrics_to_report": [
        "kl_per_dim_mean",       # Target: > 0.1
        "active_units",          # Target: > 30% of dims
        "recon_accuracy",        # Target: > 80%
        "validity_pct",          # Target: > 90%
        "uniqueness_pct",        # Target: > 80%
        "hinge_binder_coverage", # Target: > 0.5 * training set (skill 13 KPI)
        "scaffold_diversity",    # Target: > 5 distinct scaffolds in top motifs
    ],
    "hard_gates": {
        "kl_per_dim_mean": (">=", 0.05),
        "hinge_binder_coverage_ratio": (">=", 0.3),
        "validity_pct": (">=", 85.0),
    },
}
```

## Phase 3: Training Set Enrichment

### 3.1 Curate a Kinase-Focused Training Set

If the training set itself lacks kinase motifs, no amount of annealing will help.

```python
from rdkit import Chem


def enrich_training_set(original_smiles, kinase_smiles,
                        target_kinase_fraction=0.3):
    """Enrich the training set with kinase-active molecules.

    Strategy:
    1. Start with the original ChEMBL dataset
    2. Add kinase-active molecules from broader kinase datasets
    3. Target a minimum kinase motif presence

    Args:
        original_smiles: list of SMILES from original training set
        kinase_smiles: list of SMILES from kinase-focused datasets
                       (e.g., ChEMBL kinase subset, KLIFS ligands)
        target_kinase_fraction: desired fraction of kinase-like molecules

    Returns:
        enriched training set + statistics
    """
    n_original = len(original_smiles)
    n_original_kinase = sum(1 for smi in original_smiles if _has_hinge_binder(smi))
    original_frac = n_original_kinase / n_original if n_original > 0 else 0

    if original_frac >= target_kinase_fraction:
        return {
            "enriched": original_smiles,
            "n_added": 0,
            "original_kinase_frac": original_frac,
            "final_kinase_frac": original_frac,
            "note": "Training set already meets kinase fraction target",
        }

    # Calculate how many kinase molecules to add
    n_needed = int(
        (target_kinase_fraction * n_original - n_original_kinase) /
        (1 - target_kinase_fraction)
    )
    n_needed = min(n_needed, len(kinase_smiles))

    # Deduplicate
    original_set = set(original_smiles)
    additions = [smi for smi in kinase_smiles if smi not in original_set][:n_needed]
    enriched = original_smiles + additions
    n_enriched_kinase = sum(1 for smi in enriched if _has_hinge_binder(smi))

    return {
        "enriched": enriched,
        "n_original": n_original,
        "n_added": len(additions),
        "n_total": len(enriched),
        "original_kinase_frac": original_frac,
        "final_kinase_frac": n_enriched_kinase / len(enriched),
    }


def _has_hinge_binder(smi):
    """Quick check for any hinge binder motif (from skill 13)."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    patterns = [
        Chem.MolFromSmarts("c1ncnc(N)c1"),    # 2-aminopyrimidine
        Chem.MolFromSmarts("[nH]1ccnc1"),       # pyrazole NH
        Chem.MolFromSmarts("c1ccnc(N)c1"),      # 2-aminopyridine
        Chem.MolFromSmarts("c1ncnc2[nH]cnc12"), # purine-like
    ]
    return any(mol.HasSubstructMatch(p) for p in patterns if p is not None)
```

### 3.2 Data Sources for Kinase Enrichment

```
Priority sources for kinase ligand SMILES:
1. KLIFS (https://klifs.net/) — kinase-ligand structures with standardized alignment
2. ChEMBL kinase subset — filter by target_class = "Kinase", pChEMBL >= 5.0
3. BindingDB kinase — Ki/IC50 data for kinases
4. PDB co-crystal ligands — from kinase structures in the PDB

For ALK5/TGFBR1 specifically:
- ChEMBL target CHEMBL4439 (TGFBR1)
- Known inhibitor scaffolds: galunisertib (LY2157299), SB-431542, SB-505124
- Extract hinge-binding fragments from these for conditioning
```

## Phase 4: Scaffold-Conditioned Generation

### 4.1 Fragment-Conditioned VAE Architecture

The key idea: encode a hinge-binding fragment as a condition vector, concatenate with the latent code, and decode the full molecule.

```python
import torch
import torch.nn as nn


class FragmentConditionedVAE(nn.Module):
    """VAE conditioned on a hinge-binding fragment.

    At generation time:
    1. Choose a fragment (e.g., aminopyrimidine)
    2. Encode it to get the fragment embedding
    3. Sample z ~ N(0, I)
    4. Decode concat(z, frag_embed) to get a full molecule
    """

    def __init__(self, vocab_size, embed_dim=128, latent_dim=128,
                 frag_embed_dim=64, hidden_dim=512, max_len=100):
        super().__init__()
        self.latent_dim = latent_dim
        self.frag_embed_dim = frag_embed_dim
        self.embedding = nn.Embedding(vocab_size, embed_dim)

        # Fragment encoder
        self.frag_gru = nn.GRU(embed_dim, frag_embed_dim, batch_first=True)

        # Molecule encoder
        self.mol_gru = nn.GRU(embed_dim, hidden_dim, batch_first=True)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

        # Decoder (conditioned on z + fragment)
        self.decoder_init = nn.Linear(latent_dim + frag_embed_dim, hidden_dim)
        self.decoder_gru = nn.GRU(embed_dim, hidden_dim, batch_first=True)
        self.decoder_out = nn.Linear(hidden_dim, vocab_size)
        self.max_len = max_len

    def encode_fragment(self, frag_seq):
        embedded = self.embedding(frag_seq)
        _, h = self.frag_gru(embedded)
        return h.squeeze(0)  # (batch, frag_embed_dim)

    def encode(self, mol_seq):
        embedded = self.embedding(mol_seq)
        _, h = self.mol_gru(embedded)
        h = h.squeeze(0)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    def decode(self, z, frag_embed, target_seq=None):
        combined = torch.cat([z, frag_embed], dim=-1)
        h0 = self.decoder_init(combined).unsqueeze(0)

        if target_seq is not None:
            # Teacher forcing
            embedded = self.embedding(target_seq)
            output, _ = self.decoder_gru(embedded, h0)
            return self.decoder_out(output)
        else:
            # Autoregressive
            batch_size = z.size(0)
            tok = torch.zeros(batch_size, 1, dtype=torch.long, device=z.device)
            outputs = []
            hidden = h0
            for _ in range(self.max_len):
                embedded = self.embedding(tok)
                output, hidden = self.decoder_gru(embedded, hidden)
                logits = self.decoder_out(output)
                tok = logits.argmax(dim=-1)
                outputs.append(tok)
            return torch.cat(outputs, dim=1)

    def forward(self, mol_seq, frag_seq):
        frag_embed = self.encode_fragment(frag_seq)
        mu, logvar = self.encode(mol_seq)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z, frag_embed, mol_seq)
        return recon, mu, logvar
```

### 4.2 Simpler Alternative: Rejection Sampling with Existing VAE

If rewriting the full VAE is too heavy, sample from a pre-trained unconditional VAE with immediate rejection filtering:

```python
def scaffold_constrained_sampling(model, detokenizer,
                                   required_smarts,
                                   n_target=500, max_attempts=50000,
                                   temperature=1.0, device="cuda"):
    """Generate molecules with rejection sampling for required fragments.

    Less efficient than true conditioning but works with any pre-trained VAE.

    Args:
        model: trained VAE with decode(z) method
        detokenizer: function to convert token tensor to SMILES
        required_smarts: list of SMARTS that must be present
        n_target: target number of accepted molecules
        max_attempts: cap on total generation attempts

    Returns:
        dict with generated molecules and efficiency stats
    """
    from rdkit import Chem

    compiled = [Chem.MolFromSmarts(s) for s in required_smarts]
    compiled = [f for f in compiled if f is not None]

    model.eval()
    accepted = []
    n_valid, n_attempted = 0, 0

    while len(accepted) < n_target and n_attempted < max_attempts:
        batch_size = min(1000, max_attempts - n_attempted)
        z = torch.randn(batch_size, model.latent_dim).to(device) * temperature

        with torch.no_grad():
            tokens = model.decode(z)

        smiles_list = [detokenizer(t) for t in tokens]
        n_attempted += batch_size

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            n_valid += 1
            if any(mol.HasSubstructMatch(f) for f in compiled):
                accepted.append(smi)
                if len(accepted) >= n_target:
                    break

    efficiency = len(accepted) / n_attempted if n_attempted > 0 else 0

    return {
        "generated": accepted,
        "n_accepted": len(accepted),
        "n_attempted": n_attempted,
        "n_valid": n_valid,
        "efficiency": efficiency,
        "efficiency_diagnosis": (
            "good (>5%)" if efficiency > 0.05 else
            "poor (0.5-5%)" if efficiency > 0.005 else
            "terrible (<0.5%) — retrain with conditioning"
        ),
    }
```

## Phase 5: Validation and KPIs

### 5.1 Post-Generation Quality Check

```python
def validate_generation(generated_smiles, training_smiles,
                         hinge_binder_smarts=None):
    """Comprehensive quality check for a generated pool.

    This is the hard gate. If it fails, do NOT proceed to docking.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.DataStructs import TanimotoSimilarity
    import numpy as np

    # Basic quality
    n_total = len(generated_smiles)
    valid_mols = [Chem.MolFromSmiles(s) for s in generated_smiles]
    valid_mols = [m for m in valid_mols if m is not None]
    n_valid = len(valid_mols)
    valid_smiles = [Chem.MolToSmiles(m) for m in valid_mols]
    unique_smiles = set(valid_smiles)
    n_unique = len(unique_smiles)

    validity = 100.0 * n_valid / n_total if n_total > 0 else 0
    uniqueness = 100.0 * n_unique / n_valid if n_valid > 0 else 0

    # Novelty
    train_set = set(training_smiles)
    novel = [s for s in unique_smiles if s not in train_set]
    novelty = 100.0 * len(novel) / n_unique if n_unique > 0 else 0

    # Hinge binder coverage
    n_with_hinge = 0
    if hinge_binder_smarts:
        compiled = [Chem.MolFromSmarts(s) for s in hinge_binder_smarts]
        compiled = [p for p in compiled if p is not None]
        for mol in valid_mols:
            if any(mol.HasSubstructMatch(p) for p in compiled):
                n_with_hinge += 1
    hinge_coverage = 100.0 * n_with_hinge / n_valid if n_valid > 0 else 0

    # Training set hinge coverage for comparison
    train_mols = [Chem.MolFromSmiles(s) for s in training_smiles[:5000]]
    train_mols = [m for m in train_mols if m is not None]
    n_train_hinge = 0
    if hinge_binder_smarts:
        for mol in train_mols:
            if any(mol.HasSubstructMatch(p) for p in compiled):
                n_train_hinge += 1
    train_hinge_pct = 100.0 * n_train_hinge / len(train_mols) if train_mols else 0
    coverage_ratio = hinge_coverage / train_hinge_pct if train_hinge_pct > 0 else 0

    # Internal diversity (pairwise Tanimoto distance on a sample)
    sample = valid_mols[:500]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in sample]
    dists = []
    for i in range(min(100, len(fps))):
        for j in range(i + 1, min(100, len(fps))):
            dists.append(1.0 - TanimotoSimilarity(fps[i], fps[j]))
    int_diversity = float(np.mean(dists)) if dists else 0

    # Hard gate evaluation
    gates_passed = True
    gate_results = {}

    gate_results["validity"] = {"value": validity, "threshold": 85.0,
                                 "pass": validity >= 85.0}
    gate_results["hinge_coverage_ratio"] = {"value": coverage_ratio,
                                             "threshold": 0.3,
                                             "pass": coverage_ratio >= 0.3}
    if not all(g["pass"] for g in gate_results.values()):
        gates_passed = False

    return {
        "n_generated": n_total,
        "n_valid": n_valid,
        "n_unique": n_unique,
        "validity_pct": round(validity, 1),
        "uniqueness_pct": round(uniqueness, 1),
        "novelty_pct": round(novelty, 1),
        "hinge_coverage_pct": round(hinge_coverage, 1),
        "train_hinge_pct": round(train_hinge_pct, 1),
        "coverage_ratio": round(coverage_ratio, 3),
        "internal_diversity": round(int_diversity, 3),
        "gate_results": gate_results,
        "all_gates_passed": gates_passed,
    }
```

### 5.2 Cycle 2 Readiness Checklist

Before proceeding from generation to docking in Cycle 2:

```
Hard gates (ALL must pass):
[x] validity_pct >= 85%
[x] hinge_binder_coverage_ratio >= 0.3 (vs training set)
[x] KL diagnosis != "full_collapse"

Soft targets (track for improvement):
[ ] uniqueness_pct >= 80%
[ ] novelty_pct >= 50%
[ ] internal_diversity >= 0.6
[ ] coverage_ratio >= 0.5 (ideally matching training set)

If ANY hard gate fails:
-> Do NOT proceed to docking
-> Diagnose: is it collapse (Phase 1-2) or data (Phase 3) or conditioning (Phase 4)?
-> Fix and regenerate
-> Commit the failed experiment results for the audit trail
```

## Phase 6: Integration with DMTA Pipeline

### 6.1 Decision Tree for Cycle 2+

```
Starting Cycle 2 (after Cycle 1 showed generation problems):
|
+-- Step 0: Diagnose current VAE -> diagnose_vae()
|   +-- full_collapse -> Phase 2 (retrain with cyclical annealing)
|   +-- partial_collapse -> Phase 2 (fine-tune with cyclical annealing)
|   +-- healthy -> Phase 3 (check training data) or Phase 4 (add conditioning)
|
+-- Step 1: Choose strategy
|   +-- Strategy A: Retrain unconditional VAE with anti-collapse
|   |   (fastest, try first)
|   |
|   +-- Strategy B: Enrich training set + retrain
|   |   (if Strategy A passes collapse gate but fails hinge coverage)
|   |
|   +-- Strategy C: Fragment-conditioned VAE
|   |   (most powerful, requires more engineering)
|   |
|   +-- Strategy D: Rejection sampling from existing VAE
|       (quick fix, use if efficiency > 5%)
|
+-- Step 2: Generate candidate pool (n=1000-5000)
|
+-- Step 3: validate_generation() -> hard gate check
|   +-- FAIL -> go back to Step 1, try next strategy
|   +-- PASS -> proceed
|
+-- Step 4: Run chem-kinase-sar (skill 13) for detailed motif analysis
+-- Step 5: Run chem-reactivity-safety (skill 14)
+-- Step 6: Run chem-protonation-tautomer (skill 15) for docking-ready inputs
+-- Step 7: Dock (skill 9)
+-- Step 8: Run chem-docking-interactions (skill 16) for pose validation
+-- Step 9: Final ranking with composite score
```

### 6.2 Strategy Selection Guide

```
Q: Is the VAE latent space collapsed?
|
+-- YES (AU < 10%) -> Strategy A first
|   |
|   +-- After retrain, hinge coverage still < 30%?
|       +-- YES -> Strategy B (enrich data)
|       +-- NO -> proceed
|
+-- NO (AU > 10%)
    |
    +-- Hinge coverage < 30%?
    |   +-- YES, rejection sampling efficiency > 5%? -> Strategy D
    |   +-- YES, rejection sampling efficiency < 5%? -> Strategy B or C
    |   +-- NO -> proceed (generation is fine)
    |
    +-- Coverage OK but diversity low?
        +-- Increase temperature or latent sampling radius
```

## Failure Modes

1. **Fixing collapse but not fixing the training set.** A perfectly healthy VAE trained on random drug-like molecules will generate random drug-like molecules — not kinase inhibitors. Anti-collapse is necessary but not sufficient; you also need kinase-relevant training data.

2. **Over-conditioning.** If you condition on a very specific scaffold (e.g., the entire galunisertib core), you'll get back variations of galunisertib. Condition on the MINIMAL pharmacophoric fragment (e.g., just the aminopyridine that binds the hinge).

3. **Changing multiple variables at once.** Each experiment must change exactly ONE thing: annealing schedule, or free_bits value, or training data composition, or conditioning fragment. Otherwise you can't attribute improvement.

4. **Not committing failed experiments.** A failed experiment (hinge coverage still 2.5% after cyclical annealing) is valuable data. Commit it with clear metrics. The audit trail is the evidence.

5. **Using rejection sampling as a permanent solution.** If efficiency < 1%, you're wasting 99% of compute. Rejection sampling is a temporary patch, not architecture.

6. **Ignoring reconstruction accuracy.** Anti-collapse strategies can hurt reconstruction. If recon_accuracy drops below 60%, the model isn't learning molecular structure properly. Balance KL health with reconstruction quality.

## Relationship to Other Skills

| Skill | Relationship |
|-------|-------------|
| chem-kinase-sar (13) | Provides the hinge binder SMARTS and coverage KPI that this skill targets. The coverage_ratio is the primary success metric. |
| chem-reactivity-safety (14) | Run AFTER generation to filter unsafe molecules. Generation should produce chemically sensible molecules; safety screening catches the rest. |
| chem-protonation-tautomer (15) | Run AFTER generation, BEFORE docking. Generation operates on canonical SMILES; protomer enumeration is a downstream step. |
| chem-docking-interactions (16) | Run AFTER docking to validate that generated molecules actually make the expected interactions. Closes the feedback loop. |
| chem-experiment (12) | Generation experiments follow the same controlled experiment framework. Commit all results (success and failure). |

## One-Sentence Rule

**If the generator can't produce molecules with the right pharmacophore, everything downstream is wasted — fix generation first, fix generation with evidence, and never proceed to docking with a broken generator.**
