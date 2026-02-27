---
name: chem-admet
description: Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties of drug candidates. Covers multi-endpoint prediction, Lipinski/Veber rules, ADMET-specific datasets, and integration with generative models as a go/no-go filter.
homepage: https://tdcommons.ai
metadata: { "openclaw": { "emoji": "💊", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "numpy", "pandas", "sklearn", "deepchem"] } } }
---

# ADMET Prediction

Predict drug-likeness and pharmacokinetic properties. ADMET is the bridge between molecular design (skills 2-5) and real-world drug development. A beautiful molecule is worthless if it can't be absorbed, distributed, metabolized, excreted safely, and is non-toxic.

## When to Use

- User generates molecules (skill 3/4) and needs to filter by drug-likeness
- User asks about drug-likeness, bioavailability, toxicity, or pharmacokinetics
- User wants to predict specific ADMET endpoints (Caco-2, CYP450, hERG, etc.)
- User asks to apply Lipinski/Veber/lead-likeness rules
- User needs a multi-endpoint ADMET profile for a candidate set

## Core Philosophy

1. **ADMET is a filter, not an optimizer.** The goal is to eliminate bad candidates early, not to maximize ADMET scores. A molecule that passes all ADMET filters but has no activity is still useless.
2. **Rules first, ML second.** Lipinski/Veber rules are fast, interpretable, and backed by decades of data. Use them as pre-filter before running ML models.
3. **Multi-endpoint matters.** A molecule that passes hERG but fails CYP3A4 is still a problem. Always predict multiple endpoints together and report the full profile.
4. **Confidence is as important as prediction.** ADMET models are trained on specific chemical spaces. Flag when a query molecule is far from the training distribution (applicability domain).
5. **Toxicity is not binary.** "Toxic" depends on dose, route, target organ, and exposure duration. Report predictions with context, not just yes/no.

## Phase 1: Rule-Based Filters (Always Run First)

### 1.1 Lipinski's Rule of Five

```python
#!/opt/conda/envs/chem/bin/python
"""Drug-likeness rule-based filters."""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def lipinski_ro5(smi):
    """Lipinski's Rule of Five — oral bioavailability filter.
    Violation of >= 2 rules suggests poor oral absorption."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    props = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
    }

    violations = []
    if props["MW"] > 500: violations.append("MW > 500")
    if props["LogP"] > 5: violations.append("LogP > 5")
    if props["HBD"] > 5: violations.append("HBD > 5")
    if props["HBA"] > 10: violations.append("HBA > 10")

    props["violations"] = violations
    props["n_violations"] = len(violations)
    props["pass"] = len(violations) < 2
    return props
```

### 1.2 Veber Rules (Oral Bioavailability)

```python
def veber_rules(smi):
    """Veber rules — rotatable bonds and polar surface area."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    props = {
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
    }

    violations = []
    if props["RotBonds"] > 10: violations.append("RotBonds > 10")
    if props["TPSA"] > 140: violations.append("TPSA > 140")

    props["violations"] = violations
    props["pass"] = len(violations) == 0
    return props
```

### 1.3 Extended Drug-Likeness Profile

```python
def drug_likeness_profile(smi):
    """Comprehensive drug-likeness assessment."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    from rdkit.Chem import QED

    profile = {
        "smiles": smi,
        # Basic properties
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "Rings": Descriptors.RingCount(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        # Composite scores
        "QED": QED.qed(mol),  # 0-1, higher = more drug-like
        # Rule checks
        "Lipinski_pass": lipinski_ro5(smi)["pass"],
        "Veber_pass": veber_rules(smi)["pass"],
    }

    # Lead-likeness (Teague & Davis)
    profile["LeadLike"] = (
        250 <= profile["MW"] <= 350 and
        profile["LogP"] <= 3.5 and
        profile["RotBonds"] <= 7
    )

    return profile
```

### 1.4 Batch Profiling

```python
import pandas as pd

def profile_molecules(smiles_list):
    """Profile a list of molecules with all rule-based filters."""
    results = []
    for smi in smiles_list:
        p = drug_likeness_profile(smi)
        if p is not None:
            results.append(p)
    df = pd.DataFrame(results)

    print(f"\n=== Drug-Likeness Summary (N={len(df)}) ===")
    print(f"  Lipinski pass: {df['Lipinski_pass'].sum()} / {len(df)} ({df['Lipinski_pass'].mean():.1%})")
    print(f"  Veber pass:    {df['Veber_pass'].sum()} / {len(df)} ({df['Veber_pass'].mean():.1%})")
    print(f"  Lead-like:     {df['LeadLike'].sum()} / {len(df)} ({df['LeadLike'].mean():.1%})")
    print(f"  QED (mean):    {df['QED'].mean():.3f}")
    print(f"  MW range:      [{df['MW'].min():.0f}, {df['MW'].max():.0f}]")
    print(f"  LogP range:    [{df['LogP'].min():.1f}, {df['LogP'].max():.1f}]")

    return df
```

## Phase 2: ADMET Endpoint Prediction (ML)

### 2.1 Key ADMET Endpoints

| Category | Endpoint | Task Type | Dataset | Relevance |
|----------|----------|-----------|---------|-----------|
| **A**bsorption | Caco-2 permeability | regression | TDC | Intestinal absorption proxy |
| | HIA (human intestinal absorption) | classification | TDC | Oral bioavailability |
| | Pgp inhibitor | classification | TDC | Drug efflux |
| **D**istribution | BBB penetration | classification | TDC/MoleculeNet | CNS drug design |
| | PPB (plasma protein binding) | regression | TDC | Free fraction |
| | VDss | regression | TDC | Volume of distribution |
| **M**etabolism | CYP2D6 inhibitor | classification | TDC | Drug-drug interaction |
| | CYP3A4 inhibitor | classification | TDC | Most common CYP |
| | CYP2C9 inhibitor | classification | TDC | Warfarin metabolism |
| **E**xcretion | Half-life | regression | TDC | Dosing frequency |
| | Clearance | regression | TDC | Elimination rate |
| **T**oxicity | hERG inhibitor | classification | TDC | Cardiac toxicity risk |
| | AMES mutagenicity | classification | TDC/MoleculeNet | Genotoxicity |
| | LD50 | regression | TDC | Acute toxicity |
| | DILI (drug-induced liver injury) | classification | TDC | Hepatotoxicity |

### 2.2 Loading ADMET Datasets

```python
# Option A: via DeepChem / MoleculeNet
import deepchem as dc

# BBBP (blood-brain barrier)
tasks, datasets, transformers = dc.molnet.load_bbbp(featurizer="ECFP")
train, valid, test = datasets

# Tox21 (multi-label toxicity)
tasks, datasets, transformers = dc.molnet.load_tox21(featurizer="ECFP")
```

```python
# Option B: via TDC (Therapeutics Data Commons) — more endpoints
# pip install PyTDC (if not installed)
try:
    from tdc.single_pred import ADME, Tox
    # Caco-2
    data = ADME(name="Caco2_Wang")
    split = data.get_split()  # train/valid/test DataFrames

    # hERG
    data = Tox(name="hERG")
    split = data.get_split()
except ImportError:
    print("PyTDC not installed. Install with: /opt/conda/envs/chem/bin/pip install PyTDC")
```

### 2.3 Multi-Endpoint ADMET Model

```python
import torch
import torch.nn as nn
from rdkit.Chem import AllChem
import numpy as np

class ADMETPredictor(nn.Module):
    """Multi-task model for predicting multiple ADMET endpoints."""

    def __init__(self, input_dim=2048, hidden=512, endpoints=None):
        super().__init__()
        if endpoints is None:
            endpoints = {
                "Lipophilicity": "regression",
                "BBBP": "classification",
                "hERG": "classification",
                "AMES": "classification",
            }
        self.endpoints = endpoints

        # Shared backbone
        self.backbone = nn.Sequential(
            nn.Linear(input_dim, hidden),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden, hidden // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
        )

        # Per-endpoint heads
        self.heads = nn.ModuleDict()
        for name, task_type in endpoints.items():
            self.heads[name] = nn.Linear(hidden // 2, 1)

    def forward(self, x):
        shared = self.backbone(x)
        outputs = {}
        for name in self.endpoints:
            outputs[name] = self.heads[name](shared).squeeze(-1)
        return outputs

    def predict(self, smiles_list, device="cpu"):
        """Predict ADMET profile for a list of SMILES."""
        self.eval()
        fps = []
        valid_idx = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(np.array(fp))
                valid_idx.append(i)

        if not fps:
            return []

        X = torch.tensor(np.array(fps), dtype=torch.float32).to(device)
        with torch.no_grad():
            outputs = self(X)

        results = []
        for j, idx in enumerate(valid_idx):
            profile = {"smiles": smiles_list[idx]}
            for name, task_type in self.endpoints.items():
                val = outputs[name][j].item()
                if task_type == "classification":
                    profile[name] = {"prob": torch.sigmoid(torch.tensor(val)).item(),
                                     "label": int(val > 0)}
                else:
                    profile[name] = val
            results.append(profile)
        return results
```

## Phase 3: ADMET Reporting

### 3.1 ADMET Profile Table

Always produce this for candidate molecules:

```markdown
| SMILES | MW | LogP | QED | HBD | HBA | TPSA | Ro5 | Veber | BBB | hERG_risk | AMES | SA |
|--------|-----|------|-----|-----|-----|------|-----|-------|-----|-----------|------|----|
| ... | ... | ... | ... | ... | ... | ... | ✓/✗ | ✓/✗ | +/- | low/high | +/- | ... |
```

### 3.2 Traffic Light Summary

```python
def admet_traffic_light(profile):
    """Summarize ADMET profile as traffic light (green/yellow/red)."""
    signals = {}

    # Lipinski
    signals["Oral absorption"] = "🟢" if profile.get("Lipinski_pass") else "🔴"

    # QED
    qed = profile.get("QED", 0)
    if qed > 0.6: signals["Drug-likeness"] = "🟢"
    elif qed > 0.3: signals["Drug-likeness"] = "🟡"
    else: signals["Drug-likeness"] = "🔴"

    # TPSA (BBB penetration proxy)
    tpsa = profile.get("TPSA", 999)
    if tpsa < 90: signals["BBB penetration"] = "🟢"
    elif tpsa < 140: signals["BBB penetration"] = "🟡"
    else: signals["BBB penetration"] = "🔴"

    # LogP (solubility/absorption tradeoff)
    logp = profile.get("LogP", 999)
    if -0.4 <= logp <= 5: signals["Lipophilicity"] = "🟢"
    elif logp <= 6: signals["Lipophilicity"] = "🟡"
    else: signals["Lipophilicity"] = "🔴"

    return signals
```

## Phase 4: Integration with Other Skills

### 4.1 Filter Generated Molecules (connects to skill 3: chem-molgen)

```python
def admet_filter_generated(generated_smiles, strict=True):
    """Apply ADMET filters to molecules from generative models.

    Pipeline: generated → valid → drug-like → ADMET-pass
    """
    from rdkit import Chem

    # Step 1: Validity (should already be done in skill 3)
    valid = [s for s in generated_smiles if Chem.MolFromSmiles(s) is not None]

    # Step 2: Rule-based drug-likeness
    profiles = [drug_likeness_profile(s) for s in valid]
    profiles = [p for p in profiles if p is not None]

    if strict:
        # Both Lipinski + Veber
        passed = [p for p in profiles if p["Lipinski_pass"] and p["Veber_pass"]]
    else:
        # Lipinski only
        passed = [p for p in profiles if p["Lipinski_pass"]]

    # Step 3: QED filter
    final = [p for p in passed if p["QED"] > 0.3]

    print(f"\nADMET Filter Pipeline:")
    print(f"  Generated:    {len(generated_smiles)}")
    print(f"  Valid:         {len(valid)}")
    print(f"  Drug-like:    {len(passed)}")
    print(f"  QED > 0.3:    {len(final)}")
    print(f"  Pass rate:    {len(final) / max(len(generated_smiles), 1):.1%}")

    return final
```

### 4.2 Synthesizability + ADMET Combined Score (connects to skill 4)

```python
def combined_drugability_score(smi):
    """Combined score: drug-likeness + synthesizability.
    Use after generating candidates (skill 3) before retrosynthesis (skill 4).
    """
    profile = drug_likeness_profile(smi)
    if profile is None:
        return None

    from rdkit.Chem import RDConfig
    import os, sys
    sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_path not in sys.path:
        sys.path.insert(0, sa_path)
    try:
        import sascorer
        sa = sascorer.calculateScore(Chem.MolFromSmiles(smi))
    except ImportError:
        sa = 5.0  # default mid-range

    # Normalize SA to 0-1 (lower SA = easier = higher score)
    sa_norm = max(0, 1 - (sa - 1) / 9)

    profile["SA_score"] = sa
    profile["SA_normalized"] = sa_norm
    # Combined: 60% QED + 40% SA (tunable)
    profile["drugability"] = 0.6 * profile["QED"] + 0.4 * sa_norm

    return profile
```

## Checklist Before Reporting

- [ ] **Rule-based filters applied first**: Lipinski + Veber before any ML
- [ ] **QED reported**: Quantitative drug-likeness score
- [ ] **Multiple endpoints**: Not just one ADMET prediction
- [ ] **Context for toxicity**: Dose-dependent? Target organ?
- [ ] **Applicability domain checked**: Query molecule in model's training space?
- [ ] **Integration with pipeline**: Connected to generation (skill 3) and retrosynthesis (skill 4)?
- [ ] **SA Score included**: From skill 4 for synthesizability context

## Integration with Knowledge Base

- **Save profiles** to `research/ai4chem/admet/<date>-<candidate-set>.md`
- **Cross-reference** with generated molecules from chem-molgen
- **Cross-reference** with retrosynthesis assessments from chem-retrosynthesis
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "admet: <candidate-set> profiling"`
