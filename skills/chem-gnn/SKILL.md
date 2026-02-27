---
name: chem-gnn
description: Build and evaluate Graph Neural Network models for molecular property prediction. Covers molecular graph construction from SMILES, GCN/GAT/MPNN architectures, scaffold splitting for graphs, and when GNNs actually beat fingerprint baselines (and when they don't).
homepage: https://pytorch-geometric.readthedocs.io
metadata: { "openclaw": { "emoji": "🕸️", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "torch_geometric", "numpy", "pandas", "sklearn"] } } }
---

# GNN Molecular Property Prediction

Use Graph Neural Networks to learn directly from molecular graphs instead of precomputed fingerprints. This skill covers when and how to use GNNs, with emphasis on honest comparison against fingerprint baselines.

## When to Use

- User wants to predict molecular properties using graph-level learning
- Fingerprint-based models (skill 2: chem-qsar) plateau at R² < 0.6 on scaffold split
- User wants to capture 3D/topological features that fingerprints miss
- User asks about GCN, GAT, MPNN, AttentiveFP, SchNet, or similar architectures
- User wants to benchmark GNN vs fingerprint approaches

## Core Philosophy

1. **Graphs are not magic.** GNNs learn from topology + atom/bond features. If the property doesn't depend on topology (e.g., molecular weight), a GNN won't help.
2. **Beat the baseline or go home.** RF + Morgan FP is the baseline (skill 2). If GNN doesn't improve by > 5% RMSE on scaffold split, the complexity isn't justified.
3. **Message passing = learned fingerprints.** A GNN with K layers aggregates K-hop neighborhood information — conceptually similar to circular fingerprints of radius K, but learned rather than hashed.
4. **Overfitting is the default on small datasets.** Most molecular datasets are < 10K molecules. GNNs have far more parameters than RF. Regularize aggressively.
5. **Graph construction is a modeling choice.** What you encode as node features, edge features, and which bonds you include all affect results. Document every choice.

## Phase 1: Molecular Graph Construction

### 1.1 SMILES → PyG Data Object

```python
#!/opt/conda/envs/chem/bin/python
"""Convert SMILES to PyTorch Geometric Data objects."""
import torch
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from torch_geometric.data import Data

# Atom features
ATOM_FEATURES = {
    "atomic_num": list(range(1, 119)),
    "degree": [0, 1, 2, 3, 4, 5],
    "formal_charge": [-2, -1, 0, 1, 2],
    "hybridization": [
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2,
    ],
    "is_aromatic": [False, True],
}

def one_hot(value, options):
    vec = [0] * len(options)
    if value in options:
        vec[options.index(value)] = 1
    return vec

def atom_features(atom):
    """Encode atom as feature vector."""
    return (
        one_hot(atom.GetAtomicNum(), ATOM_FEATURES["atomic_num"]) +
        one_hot(atom.GetDegree(), ATOM_FEATURES["degree"]) +
        one_hot(atom.GetFormalCharge(), ATOM_FEATURES["formal_charge"]) +
        one_hot(atom.GetHybridization(), ATOM_FEATURES["hybridization"]) +
        one_hot(atom.GetIsAromatic(), ATOM_FEATURES["is_aromatic"]) +
        [atom.GetMass() / 100.0]  # normalized mass
    )

ATOM_FEAT_DIM = 118 + 6 + 5 + 5 + 2 + 1  # = 137

# Bond features
def bond_features(bond):
    """Encode bond as feature vector."""
    bt = bond.GetBondType()
    return [
        bt == Chem.rdchem.BondType.SINGLE,
        bt == Chem.rdchem.BondType.DOUBLE,
        bt == Chem.rdchem.BondType.TRIPLE,
        bt == Chem.rdchem.BondType.AROMATIC,
        bond.GetIsConjugated(),
        bond.IsInRing(),
    ]

BOND_FEAT_DIM = 6

def smiles_to_graph(smiles, y=None):
    """Convert a SMILES string to a PyG Data object.

    Returns None for invalid molecules.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Node features
    x = []
    for atom in mol.GetAtoms():
        x.append(atom_features(atom))
    x = torch.tensor(x, dtype=torch.float)

    # Edge index + edge features (undirected: add both directions)
    edge_index = []
    edge_attr = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bf = bond_features(bond)
        edge_index.extend([[i, j], [j, i]])
        edge_attr.extend([bf, bf])

    if len(edge_index) == 0:
        # Single atom molecule
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr = torch.zeros((0, BOND_FEAT_DIM), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

    if y is not None:
        data.y = torch.tensor([y], dtype=torch.float)

    data.smiles = smiles
    return data
```

### 1.2 Build Dataset

```python
from torch_geometric.data import InMemoryDataset
import pandas as pd

def build_graph_dataset(df, smiles_col="smiles", y_col="y"):
    """Convert a DataFrame to a list of PyG Data objects."""
    data_list = []
    skipped = 0
    for _, row in df.iterrows():
        data = smiles_to_graph(row[smiles_col], y=row[y_col])
        if data is not None:
            data_list.append(data)
        else:
            skipped += 1
    print(f"Built {len(data_list)} graphs, skipped {skipped} invalid SMILES")
    return data_list
```

### 1.3 Scaffold Split for Graphs

```python
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict

def scaffold_split_graphs(data_list, test_frac=0.2):
    """Scaffold split that works with PyG Data objects."""
    scaffolds = defaultdict(list)
    for i, data in enumerate(data_list):
        mol = Chem.MolFromSmiles(data.smiles)
        if mol is None:
            continue
        scaffold = MurckoScaffold.MakeScaffoldGeneric(
            MurckoScaffold.GetScaffoldForMol(mol))
        scaffolds[Chem.MolToSmiles(scaffold)].append(i)

    scaffold_groups = sorted(scaffolds.values(), key=len)
    test_idx, train_idx = [], []
    test_cutoff = int(len(data_list) * test_frac)

    for group in scaffold_groups:
        if len(test_idx) + len(group) <= test_cutoff:
            test_idx.extend(group)
        else:
            train_idx.extend(group)

    train_data = [data_list[i] for i in train_idx]
    test_data = [data_list[i] for i in test_idx]
    print(f"Scaffold split: train={len(train_data)}, test={len(test_data)}")
    return train_data, test_data
```

## Phase 2: GNN Architectures

### 2.1 GCN (Graph Convolutional Network) — Simplest Baseline

```python
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool

class MolGCN(nn.Module):
    """Simple GCN for graph-level regression/classification."""

    def __init__(self, in_dim, hidden=128, n_layers=3, dropout=0.2):
        super().__init__()
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(in_dim, hidden))
        for _ in range(n_layers - 1):
            self.convs.append(GCNConv(hidden, hidden))
        self.dropout = dropout
        self.head = nn.Sequential(
            nn.Linear(hidden, hidden // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden // 2, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv in self.convs:
            x = conv(x, edge_index)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        x = global_mean_pool(x, batch)  # graph-level readout
        return self.head(x).squeeze(-1)
```

### 2.2 GAT (Graph Attention Network)

```python
from torch_geometric.nn import GATConv

class MolGAT(nn.Module):
    """GAT with multi-head attention for molecular graphs."""

    def __init__(self, in_dim, hidden=128, heads=4, n_layers=3, dropout=0.2):
        super().__init__()
        self.convs = nn.ModuleList()
        self.convs.append(GATConv(in_dim, hidden // heads, heads=heads, dropout=dropout))
        for _ in range(n_layers - 1):
            self.convs.append(GATConv(hidden, hidden // heads, heads=heads, dropout=dropout))
        self.dropout = dropout
        self.head = nn.Sequential(
            nn.Linear(hidden, hidden // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden // 2, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv in self.convs:
            x = conv(x, edge_index)
            x = F.elu(x)
        x = global_mean_pool(x, batch)
        return self.head(x).squeeze(-1)
```

### 2.3 MPNN (Message Passing Neural Network)

```python
from torch_geometric.nn import NNConv, Set2Set

class MolMPNN(nn.Module):
    """MPNN with edge features and Set2Set readout."""

    def __init__(self, node_dim, edge_dim, hidden=128, n_layers=3):
        super().__init__()
        self.project = nn.Linear(node_dim, hidden)
        self.convs = nn.ModuleList()
        for _ in range(n_layers):
            edge_nn = nn.Sequential(
                nn.Linear(edge_dim, hidden * hidden),
            )
            self.convs.append(NNConv(hidden, hidden, edge_nn, aggr="mean"))
        self.set2set = Set2Set(hidden, processing_steps=3)
        self.head = nn.Sequential(
            nn.Linear(hidden * 2, hidden),
            nn.ReLU(),
            nn.Linear(hidden, 1)
        )

    def forward(self, data):
        x = self.project(data.x)
        for conv in self.convs:
            x = F.relu(conv(x, data.edge_index, data.edge_attr))
        x = self.set2set(x, data.batch)
        return self.head(x).squeeze(-1)
```

### 2.4 Architecture Decision Tree

```
What's your dataset size?
├── N < 1000
│   └── Use GCN (fewest params) + aggressive dropout (0.3-0.5)
│       WARNING: GNN may not beat RF here. Always compare.
├── 1000 < N < 10000
│   ├── No edge features needed → GCN or GAT
│   └── Bond type matters → MPNN (uses edge_attr)
└── N > 10000
    └── GAT or MPNN with deeper layers (4-6)
        Consider AttentiveFP or pretrained models

Always start with GCN. Only go to GAT/MPNN if GCN plateaus.
```

## Phase 3: Training

```python
from torch_geometric.loader import DataLoader
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

def train_gnn(model, train_data, test_data, device, epochs=100, lr=1e-3, batch_size=64):
    """Full training loop with validation monitoring."""
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=batch_size)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=10, factor=0.5, min_lr=1e-6)
    criterion = nn.MSELoss()

    best_test_rmse = float("inf")
    patience_counter = 0

    for epoch in range(epochs):
        # Train
        model.train()
        total_loss = 0
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            pred = model(batch)
            loss = criterion(pred, batch.y.squeeze())
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            total_loss += loss.item() * batch.num_graphs

        avg_loss = total_loss / len(train_data)

        # Evaluate
        model.eval()
        y_true, y_pred = [], []
        with torch.no_grad():
            for batch in test_loader:
                batch = batch.to(device)
                pred = model(batch)
                y_true.extend(batch.y.squeeze().cpu().tolist())
                y_pred.extend(pred.cpu().tolist())

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)
        rmse = mean_squared_error(y_true, y_pred, squared=False)
        scheduler.step(rmse)

        if rmse < best_test_rmse:
            best_test_rmse = rmse
            torch.save(model.state_dict(), "best_gnn.pt")
            patience_counter = 0
        else:
            patience_counter += 1

        if (epoch + 1) % 20 == 0:
            r2 = r2_score(y_true, y_pred)
            print(f"Epoch {epoch+1}: loss={avg_loss:.4f} RMSE={rmse:.3f} R²={r2:.3f}")

        if patience_counter >= 30:
            print(f"Early stopping at epoch {epoch+1}")
            break

    # Load best model
    model.load_state_dict(torch.load("best_gnn.pt"))
    return model
```

## Phase 4: Evaluation & Comparison

### 4.1 GNN vs Fingerprint Comparison Table

This is the **most important output**. Always produce it.

```markdown
| Model | Features | Split | RMSE | MAE | R² | RMSE/std(y) | Params |
|-------|----------|-------|------|-----|----|-------------|--------|
| RF-500 | Morgan 2048 | scaffold | X.XX | X.XX | X.XX | X.XX | — |
| RF-500 | Morgan 2048 | random | X.XX | X.XX | X.XX | X.XX | — |
| GCN-3L-128 | atom feats | scaffold | X.XX | X.XX | X.XX | X.XX | XXk |
| GCN-3L-128 | atom feats | random | X.XX | X.XX | X.XX | X.XX | XXk |
| GAT-3L-128 | atom feats | scaffold | X.XX | X.XX | X.XX | X.XX | XXk |
| MPNN-3L-128 | atom+bond | scaffold | X.XX | X.XX | X.XX | X.XX | XXk |
```

### 4.2 When GNN Wins vs Loses

```
GNN likely BEATS fingerprints when:
├── Dataset > 5K molecules
├── Property depends on 3D topology / ring systems / stereo
├── Scaffold split R² of RF < 0.6 (fingerprints aren't capturing enough)
└── Multiple related tasks (multi-task GNN can share representations)

GNN likely LOSES to fingerprints when:
├── Dataset < 1K molecules (overfitting dominates)
├── Property is mostly additive (e.g., logP ≈ sum of fragment contributions)
├── Simple linear QSAR suffices (R² > 0.8 with RF)
└── Training budget is limited (GNN needs GPU + tuning time)
```

### 4.3 Parameter Count

```python
def count_params(model):
    total = sum(p.numel() for p in model.parameters())
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Total params: {total:,} | Trainable: {trainable:,}")
    return trainable
```

## Phase 5: Advanced Techniques

### 5.1 Readout Variants

```python
from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool

# Mean pool: average over all atoms (default, good for most tasks)
# Max pool: take max per feature (captures "any atom has this feature")
# Add pool: sum over atoms (sensitive to molecule size)
# Set2Set: learnable attention-based readout (best but slowest)

def combined_readout(x, batch):
    """Concatenate mean + max pooling."""
    return torch.cat([global_mean_pool(x, batch), global_max_pool(x, batch)], dim=-1)
# Note: doubles the input dim to the prediction head
```

### 5.2 Virtual Node (improves long-range message passing)

```python
from torch_geometric.nn import global_mean_pool

class VirtualNode(nn.Module):
    """Virtual node connected to all atoms for global information flow."""

    def __init__(self, hidden):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden)
        )

    def forward(self, x, batch):
        # Aggregate all node features to virtual node
        vn = global_mean_pool(x, batch)  # [num_graphs, hidden]
        # Broadcast back to all nodes
        vn_expanded = vn[batch]  # [num_nodes, hidden]
        x = x + self.mlp(vn_expanded)
        return x
```

### 5.3 Multi-Task Learning

```python
class MultiTaskGNN(nn.Module):
    """GNN with shared backbone and separate heads per task."""

    def __init__(self, in_dim, hidden=128, n_tasks=3, n_layers=3):
        super().__init__()
        self.convs = nn.ModuleList()
        self.convs.append(GCNConv(in_dim, hidden))
        for _ in range(n_layers - 1):
            self.convs.append(GCNConv(hidden, hidden))
        # Separate head per task
        self.heads = nn.ModuleList([
            nn.Sequential(nn.Linear(hidden, hidden // 2), nn.ReLU(), nn.Linear(hidden // 2, 1))
            for _ in range(n_tasks)
        ])

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv in self.convs:
            x = F.relu(conv(x, edge_index))
        x = global_mean_pool(x, batch)
        return torch.cat([head(x) for head in self.heads], dim=-1)
```

## Checklist Before Reporting

- [ ] **RF baseline included**: Always compare GNN against RF + Morgan (skill 2)
- [ ] **Scaffold split used**: Not just random split
- [ ] **Both splits reported**: Scaffold AND random for gap analysis
- [ ] **Parameter count shown**: GNN has how many params vs RF's implicit complexity?
- [ ] **RMSE/std(y) reported**: In context of data range
- [ ] **Graph construction documented**: What atom/bond features? Directed or undirected?
- [ ] **Overfitting checked**: Train loss vs test loss gap reasonable?
- [ ] **Architecture justified**: Why GCN/GAT/MPNN for this task?

## Integration with Knowledge Base

- **Save experiments** to `research/ai4chem/experiments/gnn/<date>-<dataset>-<model>.md`
- **Compare with skill 2 results**: Cross-reference QSAR baselines
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "gnn: <dataset> <model> scaffold-split"`
