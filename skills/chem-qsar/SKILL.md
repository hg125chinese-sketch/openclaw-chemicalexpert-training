---
name: chem-qsar
description: Build and evaluate QSAR / molecular property prediction models using RDKit descriptors/fingerprints + sklearn/PyTorch. Covers featurization, splitting, training, evaluation, and the pitfalls that make most published QSAR results unreproducible.
homepage: https://www.rdkit.org
metadata: { "openclaw": { "emoji": "🧬", "requires": { "bins": ["python3"], "python": ["rdkit", "numpy", "pandas", "sklearn", "torch"] } } }
---

# QSAR / Molecular Property Prediction

Build rigorous, reproducible molecular property prediction models. This skill emphasizes **what most tutorials get wrong**: splitting, leakage, evaluation, and honest reporting.

## When to Use

- User asks to predict a molecular property (solubility, logP, toxicity, HOMO/LUMO, etc.)
- User provides a dataset of SMILES + labels
- User asks to compare fingerprints, descriptors, or model architectures
- User asks to evaluate or critique an existing QSAR model

## Core Philosophy

1. **The split is the model.** A QSAR result is only as trustworthy as its train/test split. Scaffold split by default; random split only as a sanity check.
2. **Leakage kills silently.** Feature scaling, augmentation, and hyperparameter tuning must happen AFTER splitting. Never touch test data during training.
3. **Baselines before deep learning.** Always run Random Forest on Morgan fingerprints first. If a GNN can't beat RF+Morgan by > 5% on scaffold split, it's not worth the complexity.
4. **Report what matters.** RMSE alone is meaningless without the data range. Always report RMSE, MAE, R², and the ratio RMSE/std(y). For classification: AUROC, AUPRC, and confusion matrix.
5. **Molecules are not i.i.d.** Chemical series create hidden correlations. Scaffold split approximates real-world deployment; random split approximates memorization.

## Phase 1: Data Preparation

### 1.1 Load and Validate SMILES

```python
#!/opt/conda/envs/chem/bin/python
"""Load a SMILES dataset and validate molecules."""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

df = pd.read_csv("dataset.csv")  # expects columns: smiles, y

# Validate and canonicalize
def validate_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)  # canonical SMILES

df["canon_smi"] = df["smiles"].apply(validate_smiles)
n_invalid = df["canon_smi"].isna().sum()
print(f"Total: {len(df)}, Valid: {len(df) - n_invalid}, Invalid: {n_invalid}")
df = df.dropna(subset=["canon_smi"])

# Deduplicate by canonical SMILES (keep first)
n_before = len(df)
df = df.drop_duplicates(subset=["canon_smi"], keep="first")
print(f"Deduplicated: {n_before} → {len(df)}")
```

### 1.2 Exploratory Data Analysis

Before modeling, always check:

```python
import numpy as np

y = df["y"].values
print(f"N = {len(y)}")
print(f"mean = {y.mean():.3f}, std = {y.std():.3f}")
print(f"min = {y.min():.3f}, max = {y.max():.3f}")
print(f"range = {y.max() - y.min():.3f}")

# Check for class imbalance (classification tasks)
if df["y"].nunique() <= 10:
    print(df["y"].value_counts())
```

**Red flags to watch for:**
- std(y) very small → model might look good but just predicts the mean
- Heavy class imbalance → accuracy is misleading, use AUROC/AUPRC
- Duplicates with conflicting labels → data quality issue

## Phase 2: Featurization

### 2.1 Morgan Fingerprints (default choice)

```python
from rdkit.Chem import AllChem
import numpy as np

def smiles_to_morgan(smi, radius=2, nbits=2048):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    return np.array(fp)

X = np.array([smiles_to_morgan(s) for s in df["canon_smi"]])
y = df["y"].values
print(f"X shape: {X.shape}, y shape: {y.shape}")
```

### 2.2 RDKit 2D Descriptors

```python
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

descriptor_names = [name for name, _ in Descriptors.descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

def smiles_to_descriptors(smi):
    mol = Chem.MolFromSmiles(smi)
    return np.array(calc.CalcDescriptors(mol))

X_desc = np.array([smiles_to_descriptors(s) for s in df["canon_smi"]])

# Handle NaN/inf (some descriptors fail on certain molecules)
X_desc = np.nan_to_num(X_desc, nan=0.0, posinf=0.0, neginf=0.0)
print(f"Descriptors shape: {X_desc.shape} ({len(descriptor_names)} descriptors)")
```

### 2.3 Feature Selection Guide

| Feature Type | Dim | Best For | Compute Cost |
|-------------|-----|----------|-------------|
| Morgan FP (r=2, 2048) | 2048 | General QSAR, activity cliffs | Low |
| Morgan FP (r=3, 4096) | 4096 | SAR with larger substructures | Low |
| RDKit 2D descriptors | ~210 | Physical property prediction | Low |
| Morgan + descriptors | ~2260 | Ensembles, stacking | Low |
| MACCS keys | 167 | Quick substructure screen | Very low |

**Decision rule:** Start with Morgan r=2 2048-bit. If R² < 0.6 on scaffold split, try adding 2D descriptors. If still poor, the problem may need graph-level features (GNN territory).

## Phase 3: Splitting

### 3.1 Scaffold Split (DEFAULT — always use this)

```python
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict

def scaffold_split(smiles_list, y, test_frac=0.2, seed=42):
    """Butina-style scaffold split. Molecules with same Murcko scaffold
    go to the same fold. Smaller scaffolds go to test first."""
    scaffolds = defaultdict(list)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        scaffold = MurckoScaffold.MakeScaffoldGeneric(
            MurckoScaffold.GetScaffoldForMol(mol))
        scaffolds[Chem.MolToSmiles(scaffold)].append(i)

    # Sort scaffold groups by size (smallest first → test)
    scaffold_groups = sorted(scaffolds.values(), key=len)

    test_idx, train_idx = [], []
    test_cutoff = int(len(smiles_list) * test_frac)

    for group in scaffold_groups:
        if len(test_idx) + len(group) <= test_cutoff:
            test_idx.extend(group)
        else:
            train_idx.extend(group)

    return np.array(train_idx), np.array(test_idx)

train_idx, test_idx = scaffold_split(df["canon_smi"].tolist(), y)
print(f"Train: {len(train_idx)}, Test: {len(test_idx)}")
print(f"Train scaffolds ∩ Test scaffolds should be ∅")

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]
```

### 3.2 Random Split (sanity check only)

```python
from sklearn.model_selection import train_test_split

X_train_r, X_test_r, y_train_r, y_test_r = train_test_split(
    X, y, test_size=0.2, random_state=42)
```

### 3.3 Split Quality Check

```python
# Always report distribution shift between train and test
print(f"Train y: mean={y_train.mean():.3f} std={y_train.std():.3f}")
print(f"Test  y: mean={y_test.mean():.3f} std={y_test.std():.3f}")

# Large shift is EXPECTED for scaffold split — that's the point.
# If train and test have identical distributions, your split is too easy.
```

**The split decision tree:**

```
Is this for a paper / production model?
├── YES → scaffold split (report both scaffold + random for comparison)
└── NO (quick experiment / tutorial)
    └── random split is OK, but label it clearly
```

## Phase 4: Modeling

### 4.1 Baseline: Random Forest (ALWAYS run this first)

```python
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=42)
rf.fit(X_train, y_train)

y_pred_rf = rf.predict(X_test)
rmse = mean_squared_error(y_test, y_pred_rf, squared=False)
mae = mean_absolute_error(y_test, y_pred_rf)
r2 = r2_score(y_test, y_pred_rf)

print(f"RF Baseline — RMSE: {rmse:.3f}, MAE: {mae:.3f}, R²: {r2:.3f}")
print(f"RMSE/std(y_test): {rmse / y_test.std():.3f}")  # < 0.7 is decent
```

### 4.2 Neural Network (PyTorch, GPU)

```python
import torch
import torch.nn as nn

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

class MolMLP(nn.Module):
    def __init__(self, input_dim, hidden=[512, 128], dropout=0.2):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)

# Scale features (fit on train only!)
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)  # fit on TRAIN
X_test_s = scaler.transform(X_test)        # transform TEST (no fit!)

# To tensors
X_tr = torch.tensor(X_train_s, dtype=torch.float32).to(device)
y_tr = torch.tensor(y_train, dtype=torch.float32).to(device)
X_te = torch.tensor(X_test_s, dtype=torch.float32).to(device)
y_te = torch.tensor(y_test, dtype=torch.float32).to(device)

model = MolMLP(X_train.shape[1]).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
criterion = nn.MSELoss()

# Training loop
for epoch in range(200):
    model.train()
    optimizer.zero_grad()
    loss = criterion(model(X_tr), y_tr)
    loss.backward()
    optimizer.step()

    if (epoch + 1) % 50 == 0:
        model.eval()
        with torch.no_grad():
            pred = model(X_te).cpu().numpy()
        rmse = mean_squared_error(y_test, pred, squared=False)
        print(f"Epoch {epoch+1}: train_loss={loss.item():.4f}, test_RMSE={rmse:.3f}")

# Final evaluation
model.eval()
with torch.no_grad():
    y_pred_nn = model(X_te).cpu().numpy()

rmse_nn = mean_squared_error(y_test, y_pred_nn, squared=False)
mae_nn = mean_absolute_error(y_test, y_pred_nn)
r2_nn = r2_score(y_test, y_pred_nn)
print(f"MLP — RMSE: {rmse_nn:.3f}, MAE: {mae_nn:.3f}, R²: {r2_nn:.3f}")
```

### 4.3 Model Comparison Table Template

Always produce this table at the end:

```markdown
| Model | Split | RMSE | MAE | R² | RMSE/std(y) | Notes |
|-------|-------|------|-----|----|-------------|-------|
| RF-500 | scaffold | X.XX | X.XX | X.XX | X.XX | baseline |
| RF-500 | random | X.XX | X.XX | X.XX | X.XX | sanity check |
| MLP-512-128 | scaffold | X.XX | X.XX | X.XX | X.XX | |
| MLP-512-128 | random | X.XX | X.XX | X.XX | X.XX | |
```

**Interpretation guide:**
- RMSE/std(y) < 0.7 → useful model
- RMSE/std(y) 0.7–1.0 → marginal, better than mean but barely
- RMSE/std(y) > 1.0 → worse than predicting the mean
- Scaffold R² much lower than random R² → model memorizes series, not generalizable

## Phase 5: Evaluation & Reporting

### 5.1 Regression Metrics

```python
def full_regression_report(y_true, y_pred, label="Model"):
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    ratio = rmse / y_true.std()
    print(f"\n=== {label} ===")
    print(f"N = {len(y_true)}")
    print(f"RMSE  = {rmse:.4f}")
    print(f"MAE   = {mae:.4f}")
    print(f"R²    = {r2:.4f}")
    print(f"RMSE/std(y) = {ratio:.4f}")
    print(f"y range: [{y_true.min():.2f}, {y_true.max():.2f}]")
    return {"rmse": rmse, "mae": mae, "r2": r2, "ratio": ratio}
```

### 5.2 Classification Metrics

```python
from sklearn.metrics import roc_auc_score, average_precision_score, confusion_matrix

def full_classification_report(y_true, y_score, label="Model", threshold=0.5):
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)
    y_pred = (y_score >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    print(f"\n=== {label} ===")
    print(f"AUROC = {auroc:.4f}")
    print(f"AUPRC = {auprc:.4f}")
    print(f"Confusion matrix:\n{cm}")
    return {"auroc": auroc, "auprc": auprc}
```

### 5.3 Scatter Plot (pred vs true)

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_pred_vs_true(y_true, y_pred, title="Predicted vs True", save_path="pred_vs_true.png"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_true, y_pred, alpha=0.5, s=20)
    lo = min(y_true.min(), y_pred.min())
    hi = max(y_true.max(), y_pred.max())
    ax.plot([lo, hi], [lo, hi], "r--", lw=1)
    ax.set_xlabel("True")
    ax.set_ylabel("Predicted")
    ax.set_title(title)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"Saved: {save_path}")
    plt.close(fig)
```

## Phase 6: Common Pitfalls Checklist

Before reporting any QSAR result, verify:

- [ ] **No leakage**: Scaler fit on train only? Augmentation on train only? No test molecules seen during tuning?
- [ ] **Scaffold split used**: If not, explain why and flag as "optimistic estimate"
- [ ] **Duplicates removed**: Same canonical SMILES with different labels?
- [ ] **Baseline reported**: RF result included for comparison?
- [ ] **Metrics in context**: RMSE/std(y) reported? Data range shown?
- [ ] **Scaffold vs random gap**: If scaffold R² << random R², discuss why
- [ ] **Sample size adequate**: N < 500 → results are high-variance, report with caution
- [ ] **Stereochemistry handled**: If relevant (e.g., drug potency), was chirality preserved in SMILES?

## Quick Reference: Common Datasets

| Dataset | Task | Size | Source |
|---------|------|------|--------|
| ESOL | aqueous solubility (logS) | ~1128 | DeepChem / MoleculeNet |
| FreeSolv | hydration free energy | ~643 | DeepChem / MoleculeNet |
| Lipophilicity | logD | ~4200 | DeepChem / MoleculeNet |
| BBBP | blood-brain barrier (binary) | ~2039 | MoleculeNet |
| Tox21 | toxicity (multi-label) | ~7831 | MoleculeNet |
| QM9 | quantum properties | ~134k | Various |

Load via DeepChem (already installed):

```python
import deepchem as dc
tasks, datasets, transformers = dc.molnet.load_esol(featurizer="ECFP")
train, valid, test = datasets
```

> **Warning**: DeepChem's built-in splits may differ from your scaffold split. Always verify or redo the split yourself for reproducibility.

## Integration with Knowledge Base

- **Save results** to `research/ai4chem/experiments/<dataset>/<date>-<model>.md`
- **Template** for experiment notes:

```markdown
# Experiment: <dataset> — <model>

| Field | Value |
|-------|-------|
| Date | YYYY-MM-DD |
| Dataset | name (N=...) |
| Features | Morgan r=2 2048 |
| Split | scaffold 80/20 |
| Model | RF-500 / MLP-512-128 |

## Results

| Model | Split | RMSE | MAE | R² | RMSE/std |
|-------|-------|------|-----|----|----------|
| ... | ... | ... | ... | ... | ... |

## Observations
- ...

## Next Steps
- ...
```

- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "exp: <dataset> <model> scaffold-split"`
