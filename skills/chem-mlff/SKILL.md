---
name: chem-mlff
description: Use machine-learned force fields (MACE-OFF, MACE-MP-0) for fast and accurate molecular geometry optimization, conformer generation, and energy evaluation. Replaces or supplements classical MMFF/UFF with near-DFT accuracy at classical speed. Integrates with ASE for molecular dynamics.
homepage: https://github.com/ACEsuit/mace
metadata: { "openclaw": { "emoji": "⚛️", "requires": { "bins": ["python3"], "python": ["mace-torch", "ase", "torch", "rdkit", "numpy"] } } }
---

# Machine-Learned Force Fields

Use pre-trained MACE models for molecular energy evaluation, geometry optimization, and conformer ranking at near-DFT accuracy. This replaces MMFF/UFF (RDKit defaults) when accuracy matters — binding pose refinement, strain energy calculation, or conformational analysis.

## When to Use

- User needs accurate geometry optimization (beyond MMFF)
- User wants to rank conformers by energy
- User asks about strain energy or conformational preference
- User needs to refine docking poses (skill 9 output)
- User wants molecular dynamics at near-DFT accuracy
- MMFF fails or gives suspicious geometries for complex molecules

## Core Philosophy

1. **MACE-OFF for organics, MACE-MP-0 for everything else.** MACE-OFF is trained specifically on organic molecules (CCSD(T) level). MACE-MP-0 covers 89 elements but is trained on DFT.
2. **Foundation models first, fine-tune if needed.** Start with pre-trained MACE. Only fine-tune if your chemistry is far outside training distribution.
3. **ML force fields are interpolators.** They're accurate within the training domain. Far outside it (exotic elements, extreme pressures, highly charged species), they may silently fail. Always sanity-check.
4. **ASE is the universal interface.** Atomic Simulation Environment connects MACE to optimizers, MD engines, and analysis tools. Learn ASE, not MACE internals.
5. **Compare to classical baseline.** Always report MMFF energy alongside MACE energy. If they disagree by > 50 kcal/mol, investigate.

## Phase 1: Setup & Model Loading

### 1.1 Load Pre-trained MACE Models

```python
#!/opt/conda/envs/chem/bin/python
"""Load and use MACE force fields via ASE."""
import torch
import numpy as np
from ase import Atoms
from ase.optimize import BFGS, LBFGS
from ase.io import read, write

def load_mace_calculator(model_type="mace-off", device="cpu"):
    """Load a pre-trained MACE calculator.

    Args:
        model_type: "mace-off" (organic molecules) or "mace-mp-0" (universal)
        device: "cpu" or "cuda"
    Returns:
        ASE Calculator
    """
    from mace.calculators import mace_off, mace_mp

    if model_type == "mace-off":
        # MACE-OFF: trained on organic molecules at CCSD(T) level
        # Best for drug-like molecules (H, C, N, O, F, S, Cl, Br)
        calc = mace_off(model="medium", device=device, default_dtype="float64")
        print("Loaded MACE-OFF (medium) — organic molecules, CCSD(T) accuracy")
    elif model_type == "mace-mp-0":
        # MACE-MP-0: universal, 89 elements, DFT accuracy
        calc = mace_mp(model="medium", device=device, default_dtype="float64")
        print("Loaded MACE-MP-0 (medium) — universal, DFT accuracy")
    else:
        raise ValueError(f"Unknown model: {model_type}. Use 'mace-off' or 'mace-mp-0'.")

    return calc
```

### 1.2 RDKit Mol → ASE Atoms

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def mol_to_ase(mol, conf_id=0):
    """Convert RDKit molecule (with 3D coords) to ASE Atoms."""
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no 3D conformer. Run EmbedMolecule first.")

    conf = mol.GetConformer(conf_id)
    positions = conf.GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms

def ase_to_mol(atoms, template_mol=None):
    """Convert ASE Atoms back to RDKit mol (update positions from ASE).

    Requires a template mol to preserve bond information.
    """
    if template_mol is None:
        raise ValueError("Need template mol to preserve bonds")

    mol = Chem.RWMol(template_mol)
    conf = mol.GetConformer()
    for i, pos in enumerate(atoms.positions):
        conf.SetAtomPosition(i, pos.tolist())

    return mol.GetMol()

def smiles_to_ase(smiles, optimize_mmff=True):
    """SMILES → 3D → ASE Atoms (convenience function)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        return None, None

    if optimize_mmff:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    atoms = mol_to_ase(mol)
    return atoms, mol
```

## Phase 2: Geometry Optimization

### 2.1 Optimize with MACE

```python
def optimize_geometry(smiles_or_atoms, model_type="mace-off", fmax=0.01,
                      max_steps=500, device="cpu"):
    """Optimize molecular geometry using MACE.

    Args:
        smiles_or_atoms: SMILES string or ASE Atoms
        model_type: "mace-off" or "mace-mp-0"
        fmax: force convergence criterion (eV/Å)
        max_steps: maximum optimization steps
    Returns:
        dict with optimized atoms, energy, forces, n_steps
    """
    # Load calculator
    calc = load_mace_calculator(model_type, device)

    # Prepare atoms
    if isinstance(smiles_or_atoms, str):
        atoms, mol = smiles_to_ase(smiles_or_atoms)
        if atoms is None:
            return {"error": f"Failed to embed: {smiles_or_atoms}"}
    else:
        atoms = smiles_or_atoms
        mol = None

    atoms.calc = calc

    # Get initial energy
    e_initial = atoms.get_potential_energy()

    # Optimize
    opt = LBFGS(atoms, logfile=None)
    converged = opt.run(fmax=fmax, steps=max_steps)

    e_final = atoms.get_potential_energy()
    forces = atoms.get_forces()
    max_force = np.max(np.linalg.norm(forces, axis=1))

    result = {
        "energy_initial_eV": round(float(e_initial), 4),
        "energy_final_eV": round(float(e_final), 4),
        "energy_diff_eV": round(float(e_final - e_initial), 4),
        "energy_diff_kcal": round(float((e_final - e_initial) * 23.06), 2),
        "max_force_eV_A": round(float(max_force), 6),
        "converged": bool(converged),
        "n_steps": opt.nsteps,
        "atoms": atoms,
        "mol": mol,
    }

    print(f"Optimization: {result['n_steps']} steps, "
          f"ΔE = {result['energy_diff_kcal']:.2f} kcal/mol, "
          f"max_force = {result['max_force_eV_A']:.4f} eV/Å, "
          f"converged = {result['converged']}")

    return result
```

### 2.2 Compare MACE vs MMFF

```python
def compare_mmff_vs_mace(smiles, device="cpu"):
    """Compare MMFF and MACE geometry optimization.

    Reports energy difference and RMSD between optimized geometries.
    """
    from rdkit.Chem import rdMolAlign

    # MMFF optimization
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mmff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    mmff_mol = Chem.Mol(mol)

    # MACE optimization (starting from MMFF geometry)
    atoms = mol_to_ase(mol)
    mace_result = optimize_geometry(atoms, device=device)
    mace_atoms = mace_result["atoms"]

    # Update mol with MACE positions
    mace_mol = ase_to_mol(mace_atoms, template_mol=mol)

    # RMSD between MMFF and MACE geometries
    try:
        rmsd = rdMolAlign.AlignMol(mace_mol, mmff_mol)
    except:
        rmsd = None

    comparison = {
        "smiles": smiles,
        "mmff_converged": mmff_result == 0,
        "mace_converged": mace_result["converged"],
        "mace_energy_kcal": mace_result["energy_diff_kcal"],
        "mace_n_steps": mace_result["n_steps"],
        "rmsd_mmff_vs_mace": round(rmsd, 3) if rmsd else None,
    }

    print(f"\nMMFF vs MACE comparison for: {smiles}")
    print(f"  MMFF converged: {comparison['mmff_converged']}")
    print(f"  MACE converged: {comparison['mace_converged']}")
    print(f"  RMSD (MMFF→MACE): {comparison['rmsd_mmff_vs_mace']} Å")
    print(f"  MACE ΔE: {comparison['mace_energy_kcal']} kcal/mol")

    return comparison
```

## Phase 3: Conformer Analysis

### 3.1 Generate & Rank Conformers

```python
def generate_and_rank_conformers(smiles, n_confs=50, model_type="mace-off",
                                  device="cpu"):
    """Generate conformers with RDKit, rank by MACE energy.

    Workflow: RDKit ETKDG → MMFF pre-filter → MACE single-point → rank
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate conformers
    params = AllChem.ETKDGv3()
    params.numThreads = 0
    params.pruneRmsThresh = 0.5  # Remove near-duplicates
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    print(f"Generated {len(cids)} conformers for {smiles}")

    if not cids:
        return None

    # MMFF pre-optimize all conformers
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)

    # Score with MACE
    calc = load_mace_calculator(model_type, device)
    energies = []

    for cid in cids:
        atoms = mol_to_ase(mol, conf_id=cid)
        atoms.calc = calc
        try:
            e = atoms.get_potential_energy()
            energies.append({"conf_id": cid, "energy_eV": float(e),
                           "energy_kcal": float(e * 23.06)})
        except:
            energies.append({"conf_id": cid, "energy_eV": None, "error": True})

    # Sort by energy
    valid = [e for e in energies if e.get("energy_eV") is not None]
    valid.sort(key=lambda x: x["energy_eV"])

    if valid:
        e_min = valid[0]["energy_kcal"]
        for v in valid:
            v["rel_energy_kcal"] = round(v["energy_kcal"] - e_min, 3)

    print(f"Ranked {len(valid)} conformers by MACE energy")
    if valid:
        print(f"  Lowest: {valid[0]['energy_kcal']:.2f} kcal/mol (conf {valid[0]['conf_id']})")
        print(f"  Spread: {valid[-1]['rel_energy_kcal']:.2f} kcal/mol")
        n_low = sum(1 for v in valid if v["rel_energy_kcal"] < 3.0)
        print(f"  Within 3 kcal/mol of min: {n_low}")

    return {"smiles": smiles, "n_total": len(cids), "n_ranked": len(valid),
            "conformers": valid, "mol": mol}
```

### 3.2 Strain Energy Calculation

```python
def strain_energy(smiles_or_mol, bound_conf_positions=None, model_type="mace-off",
                  device="cpu"):
    """Calculate strain energy of a molecule in a given conformation.

    Strain = E(given_conformation) - E(global_minimum_conformation)
    Useful for evaluating docking poses — high strain = unlikely binding mode.
    """
    # Get global minimum
    if isinstance(smiles_or_mol, str):
        conf_result = generate_and_rank_conformers(smiles_or_mol, n_confs=30,
                                                     model_type=model_type,
                                                     device=device)
        if conf_result is None or not conf_result["conformers"]:
            return None
        e_min = conf_result["conformers"][0]["energy_eV"]
    else:
        # Assume pre-computed minimum
        e_min = None

    if bound_conf_positions is not None:
        # Score the bound conformation
        calc = load_mace_calculator(model_type, device)
        mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
        mol = Chem.AddHs(mol)

        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        bound_atoms = Atoms(symbols=symbols, positions=bound_conf_positions)
        bound_atoms.calc = calc
        e_bound = bound_atoms.get_potential_energy()

        if e_min is not None:
            strain = (e_bound - e_min) * 23.06  # eV → kcal/mol
            print(f"Strain energy: {strain:.2f} kcal/mol")
            return {"strain_kcal": round(strain, 2), "e_bound_eV": e_bound, "e_min_eV": e_min}

    return None
```

## Phase 4: Molecular Dynamics (Short Simulations)

### 4.1 Quick MD for Conformational Sampling

```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units

def run_md(smiles, temperature_K=300, n_steps=1000, timestep_fs=1.0,
           model_type="mace-off", device="cpu"):
    """Run short MD simulation for conformational sampling.

    Returns trajectory of positions and energies.
    """
    atoms, mol = smiles_to_ase(smiles, optimize_mmff=True)
    if atoms is None:
        return None

    calc = load_mace_calculator(model_type, device)
    atoms.calc = calc

    # Set initial velocities
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K)

    # Langevin thermostat
    dyn = Langevin(atoms, timestep_fs * units.fs,
                   temperature_K=temperature_K,
                   friction=0.01 / units.fs)

    trajectory = []
    energies = []

    def record():
        trajectory.append(atoms.positions.copy())
        energies.append(float(atoms.get_potential_energy()))

    dyn.attach(record, interval=10)

    print(f"Running MD: {n_steps} steps, {temperature_K}K, {timestep_fs}fs timestep")
    dyn.run(n_steps)

    print(f"  Frames: {len(trajectory)}")
    print(f"  Energy range: [{min(energies):.2f}, {max(energies):.2f}] eV")

    return {"trajectory": trajectory, "energies": energies,
            "temperature_K": temperature_K, "n_steps": n_steps}
```

## Phase 5: Integration

### 5.1 Refine Docking Poses (skill 9 → skill 11)

```python
def refine_docking_pose(smiles, docked_positions, model_type="mace-off", device="cpu"):
    """Refine a docking pose using MACE optimization.

    Workflow: take docked pose → MACE optimize (ligand only) → report strain
    """
    # Calculate strain of docked pose
    strain_result = strain_energy(smiles, docked_positions, model_type, device)

    # Optimize from docked pose
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Set positions to docked pose
    conf = mol.GetConformer()
    for i, pos in enumerate(docked_positions[:mol.GetNumAtoms()]):
        conf.SetAtomPosition(i, pos.tolist())

    atoms = mol_to_ase(mol)
    opt_result = optimize_geometry(atoms, model_type, device=device)

    return {
        "strain_before_opt": strain_result,
        "optimization": opt_result,
    }
```

### 5.2 Decision: MMFF vs MACE

```
When to use MACE over MMFF:
├── Conformer ranking where < 1 kcal/mol matters → MACE
├── Strain energy of docking poses → MACE
├── Molecules with unusual bonding (hypervalent, strained rings) → MACE
├── Need energies comparable to DFT → MACE
├── Quick 3D embedding for 1000+ molecules → MMFF (faster)
└── Simple property calculation (no energy needed) → MMFF
```

## Checklist Before Reporting

- [ ] **Model stated**: MACE-OFF or MACE-MP-0? Why?
- [ ] **Convergence confirmed**: fmax < threshold? n_steps reasonable?
- [ ] **MMFF comparison**: Reported alongside for sanity check?
- [ ] **Strain interpreted**: < 3 kcal/mol = acceptable, > 10 = suspicious pose
- [ ] **Units correct**: eV vs kcal/mol clearly labeled (1 eV = 23.06 kcal/mol)?
- [ ] **Device stated**: CPU or GPU? Affects reproducibility of float precision.
