---
name: chem-molgen
description: Design, train, evaluate, and debug molecular generative models (VAE, Transformer, Diffusion) for de novo molecule design. Covers SMILES/SELFIES representations, training stability, evaluation metrics, and the common failure modes that make most published generative chemistry results misleading.
homepage: https://github.com/aspuru-guzik-group/selfies
metadata: { "openclaw": { "emoji": "🧪", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "selfies", "numpy", "pandas"] } } }
---

# Molecular Generative Models

Design, train, and rigorously evaluate models that generate novel molecules. This skill covers string-based (VAE, Transformer, Diffusion on SMILES/SELFIES) and graph-based approaches, with emphasis on what separates toy demos from useful generative chemistry.

## When to Use

- User asks to generate novel molecules with desired properties
- User asks to train a molecular VAE, Transformer, or Diffusion model
- User asks to evaluate a generative model (validity, novelty, diversity, property targeting)
- User asks about SMILES vs SELFIES, latent space optimization, or KL collapse
- User wants to reproduce or critique a generative chemistry paper

## Core Philosophy

1. **Validity is table stakes, not a metric.** Any model using SELFIES gets 100% validity for free. The real questions are: are the molecules novel, diverse, synthesizable, and do they hit the target property?
2. **KL collapse is the default.** Without explicit countermeasures, a VAE's decoder will ignore the latent space. Expect it, detect it, fix it — in that order.
3. **Evaluation requires a reference.** "We generated 10,000 molecules" means nothing without comparing to training set overlap, property distributions, and structural diversity against a baseline.
4. **Representation shapes the landscape.** SMILES is fragile (one wrong character → invalid molecule). SELFIES guarantees validity but has its own quirks (many-to-one mapping, dead regions). The choice is not neutral.
5. **Latent space quality > sample quality.** A smooth, interpolatable latent space that responds to property optimization is more valuable than cherry-picked pretty molecules.

## Phase 1: Representation & Data

### 1.1 SMILES vs SELFIES Decision

```
Choose SELFIES when:
├── Training a VAE (latent perturbation must produce valid molecules)
├── Doing latent space optimization (every decode must succeed)
└── Benchmarking (eliminates validity as confound)

Choose SMILES when:
├── Using a pretrained Transformer (most are SMILES-trained)
├── Doing autoregressive generation with validity filtering
└── Need precise control over stereochemistry (SELFIES support varies)
```

### 1.2 SMILES Tokenization

```python
import re

SMILES_REGEX = r"(\%\([0-9]{3}\)|\[[^\]]+\]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"

def tokenize_smiles(smi):
    return re.findall(SMILES_REGEX, smi)

def build_vocab(smiles_list):
    tokens = set()
    for smi in smiles_list:
        tokens.update(tokenize_smiles(smi))
    vocab = {"<pad>": 0, "<bos>": 1, "<eos>": 2, "<unk>": 3}
    for i, t in enumerate(sorted(tokens)):
        vocab[t] = i + 4
    return vocab

def encode_smiles(smi, vocab, max_len=120):
    tokens = tokenize_smiles(smi)
    ids = [vocab.get("<bos>")] + [vocab.get(t, vocab["<unk>"]) for t in tokens] + [vocab.get("<eos>")]
    ids = ids[:max_len]
    ids += [vocab["<pad>"]] * (max_len - len(ids))
    return ids
```

### 1.3 SELFIES Encoding

```python
import selfies as sf

def smiles_to_selfies(smi):
    """Convert SMILES to SELFIES. Returns None on failure."""
    try:
        return sf.encoder(smi)
    except sf.EncoderError:
        return None

def selfies_to_smiles(sel):
    """Convert SELFIES back to SMILES. Always produces a valid SMILES."""
    return sf.decoder(sel)

def tokenize_selfies(sel):
    """Split SELFIES into individual tokens like [C], [=O], [Branch1]."""
    return list(sf.split_selfies(sel))

def build_selfies_vocab(selfies_list):
    alphabet = sf.get_semantic_robust_alphabet()
    vocab = {"<pad>": 0, "<bos>": 1, "<eos>": 2}
    for i, sym in enumerate(sorted(alphabet)):
        vocab[sym] = i + 3
    return vocab
```

### 1.4 Dataset Preparation

```python
#!/opt/conda/envs/chem/bin/python
"""Prepare ZINC-250K subset for generative modeling."""
import pandas as pd
from rdkit import Chem

# Load (ZINC-250K or your dataset)
df = pd.read_csv("zinc250k.csv")  # column: smiles

# Canonicalize + filter
def process(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    canon = Chem.MolToSmiles(mol)
    if len(canon) > 120:  # skip very long SMILES
        return None
    return canon

df["canon"] = df["smiles"].apply(process)
df = df.dropna(subset=["canon"]).drop_duplicates(subset=["canon"])
print(f"Clean dataset: {len(df)} molecules")

# Train/val/test split (random for generative models is OK — we're not predicting properties)
from sklearn.model_selection import train_test_split
train, temp = train_test_split(df, test_size=0.1, random_state=42)
val, test = train_test_split(temp, test_size=0.5, random_state=42)
print(f"Train: {len(train)}, Val: {len(val)}, Test: {len(test)}")
```

> **Note**: For generative models, random split is acceptable because we're modeling the distribution P(molecule), not predicting a label. Scaffold split is for QSAR (Skill 2).

## Phase 2: Model Architectures

### 2.1 String VAE (SMILES/SELFIES)

```python
import torch
import torch.nn as nn

class MolVAE(nn.Module):
    """Minimal SMILES/SELFIES VAE with GRU encoder/decoder."""

    def __init__(self, vocab_size, emb_dim=128, hidden=256, latent=64, max_len=120):
        super().__init__()
        self.max_len = max_len
        self.latent = latent

        # Encoder
        self.emb = nn.Embedding(vocab_size, emb_dim, padding_idx=0)
        self.enc_gru = nn.GRU(emb_dim, hidden, batch_first=True, bidirectional=True)
        self.fc_mu = nn.Linear(hidden * 2, latent)
        self.fc_logvar = nn.Linear(hidden * 2, latent)

        # Decoder
        self.dec_input = nn.Linear(latent, hidden)
        self.dec_gru = nn.GRU(emb_dim + latent, hidden, batch_first=True)
        self.dec_out = nn.Linear(hidden, vocab_size)

    def encode(self, x):
        emb = self.emb(x)
        _, h = self.enc_gru(emb)
        h = torch.cat([h[0], h[1]], dim=-1)  # bidirectional
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, x):
        h = self.dec_input(z).unsqueeze(0)  # initial hidden
        emb = self.emb(x)
        z_expand = z.unsqueeze(1).expand(-1, emb.size(1), -1)
        dec_in = torch.cat([emb, z_expand], dim=-1)
        out, _ = self.dec_gru(dec_in, h)
        return self.dec_out(out)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        logits = self.decode(z, x)
        return logits, mu, logvar
```

### 2.2 VAE Loss with KL Annealing

```python
def vae_loss(logits, targets, mu, logvar, beta=1.0):
    """VAE ELBO loss with tunable beta for KL weight."""
    # Reconstruction: cross-entropy per token
    recon = nn.functional.cross_entropy(
        logits.view(-1, logits.size(-1)),
        targets.view(-1),
        ignore_index=0,  # pad
        reduction="mean"
    )
    # KL divergence
    kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return recon + beta * kl, recon, kl
```

### 2.3 KL Annealing Schedules

This is the **most critical hyperparameter** for molecular VAEs.

```python
def kl_anneal_cyclical(step, total_steps, n_cycles=4, ratio=0.5):
    """Cyclical annealing (Fu et al. 2019). Recommended default."""
    cycle_len = total_steps // n_cycles
    pos = step % cycle_len
    ramp = min(1.0, pos / (cycle_len * ratio))
    return ramp

def kl_anneal_monotonic(step, warmup_steps=5000):
    """Linear warmup annealing (Bowman et al. 2015)."""
    return min(1.0, step / warmup_steps)

def kl_anneal_freebits(mu, logvar, free_bits=2.0):
    """Free bits (Kingma et al. 2016). Per-dimension KL floor."""
    kl_per_dim = -0.5 * (1 + logvar - mu.pow(2) - logvar.exp())
    kl_clamped = torch.clamp(kl_per_dim, min=free_bits)
    return kl_clamped.mean()
```

**Decision tree for KL annealing:**

```
Is KL collapsing? (KL → 0, reconstruction looks good but latent is useless)
├── YES, during early training
│   └── Use cyclical annealing (n_cycles=4, ratio=0.5)
├── YES, persists after annealing
│   └── Add free bits (lambda=2.0) OR add BoW auxiliary loss
├── NO, but reconstruction is poor
│   └── Increase beta > 1.0 cautiously (beta-VAE)
└── NO, both look fine
    └── Don't touch it
```

### 2.4 KL Collapse Detection

```python
def diagnose_latent_space(model, dataloader, device):
    """Check if latent space is being used. Run after each epoch."""
    model.eval()
    all_mu, all_logvar = [], []
    with torch.no_grad():
        for batch in dataloader:
            x = batch.to(device)
            mu, logvar = model.encode(x)
            all_mu.append(mu.cpu())
            all_logvar.append(logvar.cpu())

    mu = torch.cat(all_mu)
    logvar = torch.cat(all_logvar)
    std = torch.exp(0.5 * logvar)

    # Active units: dimensions where variance of mu > 0.01
    au = (mu.var(dim=0) > 0.01).sum().item()
    avg_std = std.mean().item()
    avg_kl = (-0.5 * (1 + logvar - mu.pow(2) - logvar.exp())).mean().item()

    print(f"Active units: {au}/{mu.size(1)}")
    print(f"Avg posterior std: {avg_std:.4f}")
    print(f"Avg KL/dim: {avg_kl:.4f}")

    # Red flags:
    # - Active units < 5 → severe collapse
    # - Avg std ≈ 1.0 → posterior = prior, decoder ignoring z
    # - Avg KL/dim < 0.1 → latent not encoding information
    return {"active_units": au, "avg_std": avg_std, "avg_kl": avg_kl}
```

## Phase 3: Training Loop

```python
def train_vae(model, train_loader, val_loader, device, epochs=50,
              lr=3e-4, anneal="cyclical"):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5)

    total_steps = epochs * len(train_loader)
    step = 0
    best_val_loss = float("inf")

    for epoch in range(epochs):
        model.train()
        epoch_recon, epoch_kl = 0, 0

        for batch in train_loader:
            x = batch.to(device)
            logits, mu, logvar = model(x)

            # Teacher forcing: predict next token from current
            target = x[:, 1:]  # shift by 1
            logits = logits[:, :-1]

            # KL weight
            if anneal == "cyclical":
                beta = kl_anneal_cyclical(step, total_steps)
            elif anneal == "monotonic":
                beta = kl_anneal_monotonic(step)
            else:
                beta = 1.0

            loss, recon, kl = vae_loss(logits, target, mu, logvar, beta=beta)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_recon += recon.item()
            epoch_kl += kl.item()
            step += 1

        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch in val_loader:
                x = batch.to(device)
                logits, mu, logvar = model(x)
                target = x[:, 1:]
                logits = logits[:, :-1]
                loss, _, _ = vae_loss(logits, target, mu, logvar, beta=1.0)
                val_loss += loss.item()

        val_loss /= len(val_loader)
        scheduler.step(val_loss)

        # Diagnostics
        diag = diagnose_latent_space(model, val_loader, device)

        avg_recon = epoch_recon / len(train_loader)
        avg_kl = epoch_kl / len(train_loader)
        print(f"Epoch {epoch+1}: recon={avg_recon:.3f} kl={avg_kl:.3f} "
              f"beta={beta:.3f} val={val_loss:.3f} AU={diag['active_units']}")

        # Save best
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), "best_vae.pt")

    return model
```

## Phase 4: Generation & Decoding

### 4.1 Sample from Prior

```python
def sample_from_prior(model, vocab, n=100, device="cuda", max_len=120, temperature=1.0):
    """Sample molecules by decoding random z ~ N(0,I)."""
    model.eval()
    idx_to_token = {v: k for k, v in vocab.items()}
    molecules = []

    with torch.no_grad():
        z = torch.randn(n, model.latent).to(device) * temperature

        # Autoregressive decode
        input_ids = torch.full((n, 1), vocab["<bos>"], dtype=torch.long).to(device)

        for t in range(max_len - 1):
            logits = model.decode(z, input_ids)
            next_logits = logits[:, -1, :] / temperature
            next_token = torch.multinomial(torch.softmax(next_logits, -1), 1)
            input_ids = torch.cat([input_ids, next_token], dim=1)

    # Decode to strings
    for seq in input_ids.cpu().numpy():
        tokens = []
        for idx in seq[1:]:  # skip <bos>
            tok = idx_to_token.get(idx, "")
            if tok in ("<eos>", "<pad>"):
                break
            tokens.append(tok)
        molecules.append("".join(tokens))

    return molecules
```

### 4.2 Latent Space Interpolation

```python
def interpolate(model, smi_a, smi_b, vocab, device, n_points=10):
    """Linear interpolation between two molecules in latent space."""
    model.eval()
    with torch.no_grad():
        x_a = torch.tensor([encode_smiles(smi_a, vocab)]).to(device)
        x_b = torch.tensor([encode_smiles(smi_b, vocab)]).to(device)
        mu_a, _ = model.encode(x_a)
        mu_b, _ = model.encode(x_b)

        results = []
        for alpha in torch.linspace(0, 1, n_points):
            z = (1 - alpha) * mu_a + alpha * mu_b
            # decode z → molecule (greedy or sample)
            # ... (use sample_from_prior logic with fixed z)
            results.append(z)

    return results
```

## Phase 5: Evaluation Metrics

### 5.1 Core Metrics Suite

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np

def evaluate_generated(generated_smiles, training_smiles):
    """Comprehensive evaluation of generated molecules."""

    # 1. Validity
    valid_mols = []
    for smi in generated_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid_mols.append(Chem.MolToSmiles(mol))
    validity = len(valid_mols) / len(generated_smiles)

    # 2. Uniqueness (among valid)
    unique = set(valid_mols)
    uniqueness = len(unique) / len(valid_mols) if valid_mols else 0

    # 3. Novelty (not in training set)
    train_set = set(training_smiles)
    novel = [s for s in unique if s not in train_set]
    novelty = len(novel) / len(unique) if unique else 0

    # 4. Internal diversity (Tanimoto distance)
    if len(valid_mols) > 1:
        fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 2048)
               for s in list(unique)[:1000]]  # cap at 1000 for speed
        from rdkit import DataStructs
        sims = []
        for i in range(len(fps)):
            for j in range(i+1, min(i+50, len(fps))):  # sample pairs
                sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
        diversity = 1.0 - np.mean(sims) if sims else 0
    else:
        diversity = 0

    # 5. Scaffold diversity
    scaffolds = set()
    for smi in list(unique)[:5000]:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            scaf = MurckoScaffold.MakeScaffoldGeneric(
                MurckoScaffold.GetScaffoldForMol(mol))
            scaffolds.add(Chem.MolToSmiles(scaf))
    scaffold_diversity = len(scaffolds) / min(len(unique), 5000) if unique else 0

    results = {
        "n_generated": len(generated_smiles),
        "validity": validity,
        "uniqueness": uniqueness,
        "novelty": novelty,
        "diversity": diversity,
        "scaffold_diversity": scaffold_diversity,
    }

    print(f"\n=== Generation Evaluation ===")
    for k, v in results.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

    return results
```

### 5.2 Property Distribution Comparison

```python
def compare_property_distributions(gen_smiles, train_smiles, props=None):
    """Compare property distributions between generated and training molecules."""
    if props is None:
        props = {
            "MolWt": Descriptors.MolWt,
            "LogP": Descriptors.MolLogP,
            "TPSA": Descriptors.TPSA,
            "QED": lambda m: __import__("rdkit.Chem.QED", fromlist=["qed"]).qed(m),
            "NumRings": Descriptors.RingCount,
        }

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(props), figsize=(4 * len(props), 3))
    if len(props) == 1:
        axes = [axes]

    for ax, (name, fn) in zip(axes, props.items()):
        gen_vals = [fn(Chem.MolFromSmiles(s)) for s in gen_smiles
                    if Chem.MolFromSmiles(s) is not None]
        train_vals = [fn(Chem.MolFromSmiles(s)) for s in train_smiles[:5000]
                      if Chem.MolFromSmiles(s) is not None]

        ax.hist(train_vals, bins=50, alpha=0.5, density=True, label="Train")
        ax.hist(gen_vals, bins=50, alpha=0.5, density=True, label="Gen")
        ax.set_title(name)
        ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig("property_distributions.png", dpi=150)
    print("Saved: property_distributions.png")
    plt.close(fig)
```

### 5.3 Evaluation Reporting Table

Always produce this summary:

```markdown
| Metric | Value | Baseline (random SMILES) | Notes |
|--------|-------|-------------------------|-------|
| Validity | X.XX | 0.00 (SMILES) / 1.00 (SELFIES) | |
| Uniqueness | X.XX | ~1.00 | < 0.9 → mode collapse |
| Novelty | X.XX | ~1.00 | < 0.5 → memorizing training set |
| Diversity (1-Tanimoto) | X.XX | ~0.85 | < 0.7 → low diversity |
| Scaffold diversity | X.XX | ~0.50 | |
| Active units | X/Y | — | < 5 → KL collapse |
| Property KL divergence | X.XX | — | vs training distribution |
```

## Phase 6: Common Failure Modes & Fixes

| Symptom | Diagnosis | Fix |
|---------|-----------|-----|
| KL → 0, reconstruction perfect | KL collapse: decoder ignoring z | Cyclical annealing + free bits |
| Active units < 5 | Severe collapse | Increase free_bits to 4.0; try BoW loss |
| Validity < 50% (SMILES) | Bad tokenization or too short training | Check tokenizer; increase epochs; switch to SELFIES |
| Uniqueness < 80% | Mode collapse | Lower temperature; add diversity penalty; check data |
| Novelty < 30% | Memorization | Increase latent noise; reduce model capacity |
| Generated molecules unrealistic | Bad prior assumptions | Compare property distributions; add property predictor constraint |
| Loss NaN/explodes | Gradient issues | Clip gradients to 1.0; lower learning rate; check for empty sequences |

## Checklist Before Reporting

- [ ] **KL collapse checked**: Active units > 10? Avg KL/dim > 0.1?
- [ ] **Full metrics reported**: Validity, uniqueness, novelty, diversity, scaffold diversity
- [ ] **Property distributions compared**: Generated vs training for MolWt, LogP, TPSA, QED
- [ ] **Baseline included**: What does random sampling from training set give?
- [ ] **Representation stated**: SMILES or SELFIES? Tokenization method?
- [ ] **KL annealing described**: Which schedule? Parameters?
- [ ] **Sample molecules shown**: 10-20 random samples (not cherry-picked)
- [ ] **Latent space tested**: Interpolation between known molecules smooth?

## Integration with Knowledge Base

- **Save models** to `research/ai4chem/models/<model-name>/`
- **Save experiments** to `research/ai4chem/experiments/molgen/<date>-<model>.md`
- **Cross-reference** with paper notes in `research/ai4chem/papers/molecule-generation/`
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "molgen: <description>"`
