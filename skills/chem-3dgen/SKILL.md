---
name: chem-3dgen
description: Generate 3D molecular structures using E(3)-equivariant diffusion models. Complements chem-molgen (skill 3, SMILES/SELFIES-based) by producing molecules with explicit 3D coordinates, enabling direct docking (skill 9) and conformational analysis without separate 3D embedding.
homepage: https://github.com/ehoogeboom/e3_diffusion_for_molecules
metadata: { "openclaw": { "emoji": "🧊", "requires": { "bins": ["python3"], "python": ["torch", "torch_geometric", "e3nn", "rdkit", "numpy"] } } }
---

# 3D Molecular Generation

Generate molecules directly in 3D space using E(3)-equivariant diffusion models. Unlike SMILES/SELFIES generators (skill 3), these produce atom positions and types simultaneously — the output is a point cloud that respects rotational/translational symmetry.

## When to Use

- User needs 3D molecular structures (not just SMILES graphs)
- User wants to generate molecules conditioned on a protein pocket (structure-based drug design)
- User's pipeline goes directly from generation → docking (skill 9) without 3D embedding step
- User wants to explore conformational diversity of generated molecules
- Skill 3 (SELFIES VAE) produces valid SMILES but 3D embedding fails for complex structures

## Core Philosophy

1. **3D is not an afterthought.** Generating in 3D avoids the lossy SMILES→3D conversion. Binding depends on shape, not topology alone.
2. **Equivariance is non-negotiable.** The model must produce the same distribution regardless of how you rotate/translate the input. E(3) equivariance guarantees this by construction.
3. **Diffusion = iterative refinement.** Start from noise, denoise step by step. Each step makes the molecule slightly more realistic. This is fundamentally different from one-shot VAE generation.
4. **Validity requires post-processing.** Raw point clouds need bond assignment, sanitization, and filtering. A "generated molecule" isn't real until RDKit can parse it.
5. **Compare to skill 3, don't replace it.** SELFIES VAE is faster and easier to train. 3D diffusion is better when shape matters (pocket-conditioned design, scaffold hopping in 3D).

## Phase 1: Equivariant Building Blocks

### 1.1 E(3)-Equivariant Layers (e3nn)

```python
#!/opt/conda/envs/chem/bin/python
"""Core equivariant building blocks using e3nn."""
import torch
import torch.nn as nn
from e3nn import o3
from e3nn.nn import FullyConnectedNet, Gate
from torch_geometric.nn import radius_graph

class EquivariantBlock(nn.Module):
    """Single equivariant message-passing block.

    Uses e3nn tensor products for SO(3)-equivariant updates.
    Scalar features (l=0) carry chemical information.
    Vector features (l=1) carry directional information.
    """
    def __init__(self, irreps_in, irreps_out, irreps_edge="1x0e + 1x1o",
                 num_radial=16, cutoff=5.0):
        super().__init__()
        self.cutoff = cutoff

        # Tensor product: node features ⊗ edge features → output
        self.tp = o3.FullyConnectedTensorProduct(
            irreps_in, irreps_edge, irreps_out, shared_weights=False
        )
        # Radial network predicts tensor product weights from distance
        self.radial_net = FullyConnectedNet(
            [num_radial, 64, self.tp.weight_numel], act=torch.nn.SiLU()
        )
        # Radial basis functions
        self.num_radial = num_radial

    def radial_basis(self, dist):
        """Bessel radial basis functions."""
        n = torch.arange(1, self.num_radial + 1, device=dist.device, dtype=dist.dtype)
        return torch.sin(n * torch.pi * dist.unsqueeze(-1) / self.cutoff) / dist.unsqueeze(-1)

    def forward(self, node_features, pos, edge_index):
        src, dst = edge_index
        rel_pos = pos[dst] - pos[src]
        dist = rel_pos.norm(dim=-1, keepdim=False)

        # Edge features: spherical harmonics of direction + radial basis
        edge_sh = o3.spherical_harmonics("1x0e + 1x1o", rel_pos, normalize=True)
        radial_features = self.radial_basis(dist)
        weights = self.radial_net(radial_features)

        # Message: tensor product with learned weights
        messages = self.tp(node_features[src], edge_sh, weights)

        # Aggregate
        out = torch.zeros_like(node_features)
        out.index_add_(0, dst, messages)
        return out
```

### 1.2 Atom Type Encoding

```python
# Standard atom types for drug-like molecules
ATOM_TYPES = ['H', 'C', 'N', 'O', 'F', 'S', 'Cl', 'Br', 'P']
NUM_ATOM_TYPES = len(ATOM_TYPES)

def one_hot_atoms(atom_indices, num_types=NUM_ATOM_TYPES):
    """One-hot encode atom types."""
    return torch.nn.functional.one_hot(atom_indices, num_types).float()

def atoms_from_one_hot(one_hot):
    """Decode one-hot to atom type indices (argmax)."""
    return one_hot.argmax(dim=-1)
```

## Phase 2: Diffusion Process

### 2.1 Equivariant Diffusion Model (EDM) Core

```python
class EquivariantDiffusion(nn.Module):
    """E(3)-equivariant diffusion model for 3D molecules.

    Jointly diffuses:
    - Continuous: atom positions (x ∈ R^{N×3})
    - Categorical: atom types (h ∈ R^{N×K}, one-hot)

    The model predicts noise ε for positions and types at each timestep.
    """
    def __init__(self, num_atom_types=NUM_ATOM_TYPES, hidden_dim=128,
                 n_layers=6, num_timesteps=1000, cutoff=5.0):
        super().__init__()
        self.num_timesteps = num_timesteps
        self.num_atom_types = num_atom_types

        # Noise schedule (polynomial)
        self.register_buffer('betas', self._polynomial_schedule(num_timesteps))
        alphas = 1.0 - self.betas
        self.register_buffer('alphas_cumprod', torch.cumprod(alphas, dim=0))

        # Input embedding: atom types → hidden
        self.atom_embed = nn.Linear(num_atom_types, hidden_dim)
        self.time_embed = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # Equivariant layers
        irreps_hidden = f"{hidden_dim}x0e"
        self.layers = nn.ModuleList([
            EquivariantBlock(irreps_hidden, irreps_hidden, cutoff=cutoff)
            for _ in range(n_layers)
        ])

        # Output heads
        self.pos_head = nn.Linear(hidden_dim, 3)  # predict position noise
        self.type_head = nn.Linear(hidden_dim, num_atom_types)  # predict type noise

    def _polynomial_schedule(self, T, s=1e-5, power=2):
        """Polynomial noise schedule (EDM-style)."""
        steps = torch.linspace(0, 1, T + 1)
        alphas = (1 - steps**power) * (1 - s) + s
        betas = 1 - alphas[1:] / alphas[:-1]
        return betas.clamp(0, 0.999)

    def forward_diffusion(self, x, h, t):
        """Add noise at timestep t.

        q(x_t | x_0) = N(sqrt(alpha_bar_t) * x_0, (1-alpha_bar_t) * I)
        """
        alpha_bar = self.alphas_cumprod[t].view(-1, 1, 1)
        noise_x = torch.randn_like(x)
        noise_h = torch.randn_like(h)

        x_t = torch.sqrt(alpha_bar) * x + torch.sqrt(1 - alpha_bar) * noise_x
        h_t = torch.sqrt(alpha_bar) * h + torch.sqrt(1 - alpha_bar) * noise_h

        # Center positions (translation invariance)
        x_t = x_t - x_t.mean(dim=-2, keepdim=True)

        return x_t, h_t, noise_x, noise_h

    def predict_noise(self, x_t, h_t, t, batch=None):
        """Predict noise from noisy input."""
        # Build graph from positions
        if batch is None:
            batch = torch.zeros(x_t.size(0), dtype=torch.long, device=x_t.device)
        edge_index = radius_graph(x_t, r=5.0, batch=batch)

        # Embed
        node_feat = self.atom_embed(h_t)
        t_emb = self.time_embed(t.float().view(-1, 1) / self.num_timesteps)
        # Broadcast time to all atoms in batch
        node_feat = node_feat + t_emb[batch]

        # Message passing
        for layer in self.layers:
            node_feat = node_feat + layer(node_feat, x_t, edge_index)

        # Predict noise
        eps_x = self.pos_head(node_feat)
        eps_h = self.type_head(node_feat)

        # Center position noise (maintain zero CoM)
        eps_x = eps_x - eps_x.mean(dim=0, keepdim=True)

        return eps_x, eps_h
```

### 2.2 Training Loop

```python
def train_edm(model, dataloader, epochs=100, lr=1e-4, device="cuda"):
    """Train equivariant diffusion model."""
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-12)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    model.to(device)
    model.train()

    for epoch in range(epochs):
        total_loss = 0
        for batch in dataloader:
            x = batch.pos.to(device)       # [N_total, 3]
            h = batch.x.to(device)          # [N_total, K]
            b = batch.batch.to(device)      # [N_total]

            # Sample random timesteps per molecule
            n_mols = b.max().item() + 1
            t = torch.randint(0, model.num_timesteps, (n_mols,), device=device)
            t_per_atom = t[b]

            # Forward diffusion
            x_t, h_t, noise_x, noise_h = model.forward_diffusion(x, h, t_per_atom)

            # Predict noise
            eps_x_pred, eps_h_pred = model.predict_noise(x_t, h_t, t_per_atom, b)

            # Loss: MSE on noise prediction
            loss_x = ((eps_x_pred - noise_x) ** 2).mean()
            loss_h = ((eps_h_pred - noise_h) ** 2).mean()
            loss = loss_x + 0.4 * loss_h  # position loss dominates

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            total_loss += loss.item()

        scheduler.step()
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{epochs}, Loss: {total_loss/len(dataloader):.4f}")

    return model
```

## Phase 3: Sampling (Reverse Diffusion)

### 3.1 Generate 3D Molecules

```python
@torch.no_grad()
def sample_molecules(model, n_mols, n_atoms_per_mol=20, device="cuda"):
    """Sample molecules via reverse diffusion.

    Args:
        model: trained EquivariantDiffusion
        n_mols: number of molecules to generate
        n_atoms_per_mol: atoms per molecule (can vary)
    Returns:
        list of (positions, atom_types) tuples
    """
    model.eval()
    N = n_mols * n_atoms_per_mol
    batch = torch.arange(n_mols, device=device).repeat_interleave(n_atoms_per_mol)

    # Start from noise
    x = torch.randn(N, 3, device=device)
    x = x - x.mean(dim=0, keepdim=True)  # center
    h = torch.randn(N, model.num_atom_types, device=device)

    # Reverse diffusion
    for t in reversed(range(model.num_timesteps)):
        t_tensor = torch.full((N,), t, device=device, dtype=torch.long)
        eps_x, eps_h = model.predict_noise(x, h, t_tensor, batch)

        alpha = 1 - model.betas[t]
        alpha_bar = model.alphas_cumprod[t]

        # DDPM update
        x = (1 / torch.sqrt(alpha)) * (x - (1 - alpha) / torch.sqrt(1 - alpha_bar) * eps_x)
        h = (1 / torch.sqrt(alpha)) * (h - (1 - alpha) / torch.sqrt(1 - alpha_bar) * eps_h)

        # Add noise (except at t=0)
        if t > 0:
            noise_x = torch.randn_like(x) * torch.sqrt(model.betas[t])
            noise_h = torch.randn_like(h) * torch.sqrt(model.betas[t])
            x = x + noise_x
            h = h + noise_h

        # Recenter
        x = x - x.mean(dim=0, keepdim=True)

    # Decode atom types
    atom_types = atoms_from_one_hot(h)

    # Split into individual molecules
    molecules = []
    for i in range(n_mols):
        mask = batch == i
        mol_pos = x[mask].cpu().numpy()
        mol_types = atom_types[mask].cpu().numpy()
        molecules.append((mol_pos, mol_types))

    return molecules
```

## Phase 4: Point Cloud → RDKit Molecule

### 4.1 Bond Assignment from 3D Coordinates

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial.distance import cdist

# Covalent radii (Angstrom) for bond detection
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
    'F': 0.57, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'P': 1.07,
}

def point_cloud_to_mol(positions, atom_type_indices, bond_margin=0.3):
    """Convert 3D point cloud to RDKit molecule with bonds.

    Bond assignment: distance < (r_i + r_j) * (1 + margin)
    """
    atom_symbols = [ATOM_TYPES[i] for i in atom_type_indices]
    n_atoms = len(positions)

    # Distance matrix
    dists = cdist(positions, positions)

    # Build editable molecule
    mol = Chem.RWMol()
    for symbol in atom_symbols:
        mol.AddAtom(Chem.Atom(symbol))

    # Assign bonds
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            r_i = COVALENT_RADII.get(atom_symbols[i], 0.77)
            r_j = COVALENT_RADII.get(atom_symbols[j], 0.77)
            threshold = (r_i + r_j) * (1 + bond_margin)
            if dists[i, j] < threshold:
                mol.AddBond(i, j, Chem.BondType.SINGLE)

    # Try to sanitize and assign bond orders
    try:
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except:
        return None

def validate_generated_3d(molecules):
    """Validate a batch of generated 3D molecules."""
    results = {"total": len(molecules), "valid": 0, "unique_smiles": set()}

    valid_mols = []
    for pos, types in molecules:
        mol = point_cloud_to_mol(pos, types)
        if mol is not None:
            smi = Chem.MolToSmiles(mol)
            results["valid"] += 1
            results["unique_smiles"].add(smi)
            valid_mols.append({"mol": mol, "smiles": smi, "pos": pos})

    results["validity"] = results["valid"] / results["total"]
    results["uniqueness"] = len(results["unique_smiles"]) / max(results["valid"], 1)
    results["unique_smiles"] = list(results["unique_smiles"])

    print(f"3D Generation Metrics:")
    print(f"  Valid: {results['valid']}/{results['total']} ({results['validity']:.1%})")
    print(f"  Unique: {len(results['unique_smiles'])} ({results['uniqueness']:.1%})")

    return results, valid_mols
```

## Phase 5: Integration with Pipeline

### 5.1 3D Gen → Docking (skip 3D embedding)

```python
def gen3d_to_docking(model, receptor_pdbqt, box, n_generate=100, device="cuda"):
    """Generate 3D molecules and dock directly (no SMILES→3D conversion).

    This is the key advantage over skill 3: molecules are born in 3D.
    """
    # Generate
    molecules = sample_molecules(model, n_generate, device=device)
    results, valid_mols = validate_generated_3d(molecules)

    # For each valid molecule, write PDB and convert to PDBQT for docking
    docking_ready = []
    for i, mol_info in enumerate(valid_mols):
        mol = mol_info["mol"]
        try:
            # Add hydrogens and set 3D coords from generated positions
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())  # re-embed for clean geometry
            AllChem.MMFFOptimizeMolecule(mol)

            docking_ready.append({
                "smiles": mol_info["smiles"],
                "mol": mol,
                "source": "3d_diffusion",
            })
        except:
            continue

    print(f"  Docking-ready: {len(docking_ready)}/{len(valid_mols)}")
    return docking_ready
```

### 5.2 When to Use 3D Gen vs SELFIES VAE (Skill 3)

```
Decision tree: Which generator?
├── Need molecules for a specific protein pocket?
│   └── Use 3D diffusion (pocket-conditioned if available) or 3D gen + dock
├── Need diverse scaffold exploration (no pocket)?
│   └── Use SELFIES VAE (skill 3) — faster, easier
├── RDKit 3D embedding fails for > 10% of skill 3 output?
│   └── Use 3D diffusion — generates valid 3D by construction
├── Training data < 1000 molecules?
│   └── Use SELFIES VAE — diffusion needs more data
└── Default?
    └── Use SELFIES VAE for speed, then 3D diffusion for top candidates
```

## Checklist Before Reporting

- [ ] **Equivariance verified**: Model output invariant to rotation/translation?
- [ ] **Bond assignment validated**: Point cloud → mol → RDKit sanitize?
- [ ] **Metrics reported**: Validity / uniqueness / novelty / atom stability?
- [ ] **Compared to skill 3**: Is 3D gen actually better for this task?
- [ ] **CoM centered**: Positions centered at zero throughout diffusion?
- [ ] **Training logged**: Loss curves, n_epochs, noise schedule documented?
