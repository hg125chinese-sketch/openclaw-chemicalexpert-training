---
name: chem-docking
description: Perform molecular docking and virtual screening using AutoDock Vina. Covers receptor/ligand preparation, docking box setup, scoring interpretation, virtual screening pipelines, and integration with generative models as a binding affinity filter.
homepage: https://github.com/ccsb-scripps/AutoDock-Vina
metadata: { "openclaw": { "emoji": "🔬", "requires": { "bins": ["python3", "vina"], "python": ["rdkit", "numpy", "pandas"] } } }
---

# Molecular Docking & Virtual Screening

Predict how small molecules bind to protein targets using AutoDock Vina. This skill connects molecular design (skills 2-6) to biological activity — a molecule's predicted binding pose and affinity can validate or reject candidates before synthesis.

## When to Use

- User has a protein target (PDB ID or structure) and wants to dock candidate molecules
- User wants to virtual-screen a library of generated molecules against a target
- User asks about binding affinity, binding pose, or protein-ligand interactions
- User wants to filter generated molecules (skill 3) by predicted binding
- User needs to set up a docking workflow from scratch

## Core Philosophy

1. **Docking scores are ranks, not energies.** Vina scores correlate roughly with binding affinity but are NOT accurate free energy predictions. Use them to rank compounds, not to predict Kd.
2. **Garbage in, garbage out.** Receptor preparation (protonation, missing residues, water molecules) determines 80% of docking quality. Spend time here.
3. **The binding site must be defined.** Vina needs a search box. Too small → misses poses. Too large → noise. Use known ligand or literature to guide box placement.
4. **Poses matter more than scores.** A compound with score -8.0 in a wrong pose is worse than -7.0 in the right pose. Always visually inspect top hits.
5. **Virtual screening is a funnel, not a filter.** Dock many → rank by score → inspect top 50-100 poses → apply ADMET (skill 6) → select 5-10 for synthesis.

## Phase 1: Receptor Preparation

### 1.1 Fetch PDB Structure

```python
#!/opt/conda/envs/chem/bin/python
"""Fetch and prepare a protein structure for docking."""
import os
import urllib.request

def fetch_pdb(pdb_id, output_dir="."):
    """Download a PDB file from RCSB."""
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filepath = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.exists(filepath):
        print(f"Already exists: {filepath}")
        return filepath

    print(f"Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, filepath)
    print(f"Saved: {filepath}")
    return filepath
```

### 1.2 Clean Receptor (Remove Water, Ligands, Select Chain)

```python
def clean_receptor(pdb_path, output_path=None, chain="A", keep_het=False):
    """Clean a PDB file for docking.

    - Keep only specified chain
    - Remove water (HOH)
    - Optionally remove all HETATM (ligands, ions)
    """
    if output_path is None:
        base = os.path.splitext(pdb_path)[0]
        output_path = f"{base}_clean.pdb"

    kept_lines = []
    removed_ligands = []

    with open(pdb_path) as f:
        for line in f:
            record = line[:6].strip()
            if record == "HETATM" and line[17:20].strip() == "HOH":
                continue
            if record in ("ATOM", "HETATM"):
                if line[21] != chain:
                    continue
            if record == "HETATM" and not keep_het:
                lig_name = line[17:20].strip()
                if lig_name not in removed_ligands:
                    removed_ligands.append(lig_name)
                continue
            if record in ("ATOM", "TER", "END"):
                kept_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(kept_lines)

    print(f"Cleaned receptor: {output_path}")
    print(f"  Chain: {chain}, Removed ligands: {removed_ligands}")
    return output_path
```

### 1.3 Extract Co-crystallized Ligand (for box definition)

```python
def extract_ligand(pdb_path, ligand_name, chain="A", output_path=None):
    """Extract a co-crystallized ligand from PDB for box definition."""
    if output_path is None:
        output_path = f"{ligand_name}_from_pdb.pdb"

    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() == "HETATM":
                if line[17:20].strip() == ligand_name and line[21] == chain:
                    lines.append(line)

    if not lines:
        print(f"Ligand {ligand_name} not found in chain {chain}")
        return None

    with open(output_path, "w") as f:
        f.writelines(lines)
        f.write("END\n")

    print(f"Extracted ligand {ligand_name}: {output_path} ({len(lines)} atoms)")
    return output_path
```

### 1.4 Convert Receptor to PDBQT

```python
def receptor_to_pdbqt(pdb_path, output_path=None):
    """Convert cleaned receptor PDB to PDBQT format for Vina."""
    if output_path is None:
        output_path = pdb_path.replace(".pdb", ".pdbqt")

    import subprocess
    mk_prep = os.path.expanduser("~/.local/bin/mk_prepare_receptor.py")

    if os.path.exists(mk_prep):
        cmd = ["/opt/conda/envs/chem/bin/python", mk_prep,
               "-i", pdb_path, "-o", output_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Receptor PDBQT: {output_path}")
            return output_path
        else:
            print(f"mk_prepare_receptor failed: {result.stderr[:200]}")

    # Fallback: minimal PDBQT conversion
    print("Falling back to minimal PDBQT conversion...")
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line[:4] == "ATOM":
                atom_name = line[12:16].strip()
                element = atom_name[0] if atom_name[0].isalpha() else atom_name[1]
                pdbqt_line = line[:54] + "  0.00  0.00" + f"    {element:>2s}" + "\n"
                lines.append(pdbqt_line)

    with open(output_path, "w") as f:
        f.writelines(lines)

    print(f"Receptor PDBQT (minimal): {output_path}")
    return output_path
```

## Phase 2: Ligand Preparation

### 2.1 SMILES → 3D → PDBQT (via Meeko)

```python
from rdkit import Chem
from rdkit.Chem import AllChem
# Optional: pip install meeko. Fallback: use obabel CLI for PDBQT conversion
import meeko

def prepare_ligand(smiles, output_path="ligand.pdbqt"):
    """Convert SMILES to docking-ready PDBQT via RDKit 3D + Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return None

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)
        if result != 0:
            print(f"3D embedding failed for: {smiles}")
            return None

    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    preparator = meeko.MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        print(f"Meeko preparation failed for: {smiles}")
        return None

    pdbqt_string = meeko.PDBQTWriterLegacy.write_string(mol_setups[0])
    with open(output_path, "w") as f:
        f.write(pdbqt_string[0])

    print(f"Ligand PDBQT: {output_path}")
    return output_path
```

### 2.2 Batch Ligand Preparation

```python
def prepare_ligands_batch(smiles_list, output_dir="ligands"):
    """Prepare multiple ligands for virtual screening."""
    os.makedirs(output_dir, exist_ok=True)
    prepared, failed = [], []

    for i, smi in enumerate(smiles_list):
        out_path = os.path.join(output_dir, f"lig_{i:04d}.pdbqt")
        result = prepare_ligand(smi, out_path)
        if result:
            prepared.append({"idx": i, "smiles": smi, "pdbqt": out_path})
        else:
            failed.append({"idx": i, "smiles": smi, "error": "preparation_failed"})

    print(f"\nLigand preparation: {len(prepared)} success, {len(failed)} failed")
    return prepared, failed
```

## Phase 3: Docking Box Setup

### 3.1 Box from Co-crystallized Ligand

```python
import numpy as np

def box_from_ligand(ligand_pdb_path, padding=8.0):
    """Define docking box centered on a co-crystallized ligand."""
    coords = []
    with open(ligand_pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    if not coords:
        return None

    coords = np.array(coords)
    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)
    size = extent + 2 * padding

    box = {
        "center_x": round(float(center[0]), 2),
        "center_y": round(float(center[1]), 2),
        "center_z": round(float(center[2]), 2),
        "size_x": round(float(size[0]), 1),
        "size_y": round(float(size[1]), 1),
        "size_z": round(float(size[2]), 1),
    }
    print(f"Docking box (padding={padding}Å): center=({box['center_x']}, {box['center_y']}, {box['center_z']}), size=({box['size_x']}, {box['size_y']}, {box['size_z']})")
    return box
```

### 3.2 Box from Residue List

```python
def box_from_residues(pdb_path, residue_ids, chain="A", padding=8.0):
    """Define docking box centered on specific residues."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:4] == "ATOM" and line[21] == chain:
                resid = int(line[22:26])
                if resid in residue_ids:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])

    if not coords:
        return None

    coords = np.array(coords)
    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)
    size = extent + 2 * padding

    # Ensure minimum box size
    box = {
        "center_x": round(float(center[0]), 2),
        "center_y": round(float(center[1]), 2),
        "center_z": round(float(center[2]), 2),
        "size_x": round(max(float(size[0]), 25.0), 1),
        "size_y": round(max(float(size[1]), 25.0), 1),
        "size_z": round(max(float(size[2]), 25.0), 1),
    }
    print(f"Docking box (from {len(residue_ids)} residues): center=({box['center_x']}, {box['center_y']}, {box['center_z']}), size=({box['size_x']}, {box['size_y']}, {box['size_z']})")
    return box
```

### 3.3 Box Guidelines

```
Box sizing rules:
├── Minimum: 15×15×15 Å (very small binding pocket)
├── Typical: 22×22×22 Å (standard drug-like pocket)
├── Maximum recommended: 30×30×30 Å (large / flexible site)
└── > 40 Å on any axis → too large, docking noise increases dramatically

Padding from reference ligand:
├── 6 Å: tight — only catches similar binding modes
├── 8 Å: standard (DEFAULT) — allows some flexibility
└── 10-12 Å: loose — for exploring alternative poses
```

## Phase 4: Running Docking

**Default workflow:** use **vina CLI** (robust, easy to time out, stable in batch mode).

**Optional / advanced:** Python bindings (`vina` package). Use only if you need in-process control and you have the dependency installed.

### 4.1 Single Docking (Python API)

```python
# Optional: pip install vina. Default workflow uses vina CLI instead
from vina import Vina

def dock_single(receptor_pdbqt, ligand_pdbqt, box, exhaustiveness=32, n_poses=9):
    """Dock a single ligand to a receptor."""
    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(
        center=[box["center_x"], box["center_y"], box["center_z"]],
        box_size=[box["size_x"], box["size_y"], box["size_z"]]
    )

    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = v.energies()

    results = []
    for i, energy_row in enumerate(energies):
        results.append({
            "pose": i + 1,
            "score_kcal_mol": round(energy_row[0], 2),
            "inter": round(energy_row[1], 2) if len(energy_row) > 1 else None,
            "intra": round(energy_row[2], 2) if len(energy_row) > 2 else None,
        })

    output_path = ligand_pdbqt.replace(".pdbqt", "_docked.pdbqt")
    v.write_poses(output_path, n_poses=n_poses, overwrite=True)

    print(f"Docking: {len(results)} poses, best={results[0]['score_kcal_mol']} kcal/mol")
    return results, output_path
```

### 4.2 Virtual Screening (Batch Docking)

```python
import time

def virtual_screen(receptor_pdbqt, ligand_list, box, exhaustiveness=16, n_poses=3):
    """Screen a library of ligands against a receptor.

    Args:
        receptor_pdbqt: path to receptor PDBQT
        ligand_list: list of dicts with 'smiles' and 'pdbqt' keys
        box: docking box dict
        exhaustiveness: 16 for VS speed, 32 for accuracy
    Returns:
        sorted list of results (best score first)
    """
    results = []
    t0 = time.time()

    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(
        center=[box["center_x"], box["center_y"], box["center_z"]],
        box_size=[box["size_x"], box["size_y"], box["size_z"]]
    )

    for i, lig in enumerate(ligand_list):
        try:
            v.set_ligand_from_file(lig["pdbqt"])
            v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
            energies = v.energies()
            best_score = energies[0][0]

            out_path = lig["pdbqt"].replace(".pdbqt", "_docked.pdbqt")
            v.write_poses(out_path, n_poses=n_poses, overwrite=True)

            results.append({
                "idx": lig.get("idx", i),
                "smiles": lig["smiles"],
                "best_score": round(best_score, 2),
                "n_poses": len(energies),
                "docked_pdbqt": out_path,
            })
        except Exception as e:
            results.append({
                "idx": lig.get("idx", i),
                "smiles": lig["smiles"],
                "best_score": None,
                "error": str(e)[:100],
            })

        if (i + 1) % 50 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed * 3600
            print(f"  [{i+1}/{len(ligand_list)}] {elapsed:.0f}s elapsed, ~{rate:.0f} lig/h")

    elapsed = time.time() - t0
    scored = [r for r in results if r["best_score"] is not None]
    scored.sort(key=lambda x: x["best_score"])

    print(f"\nVirtual screening complete:")
    print(f"  Total: {len(ligand_list)}, Scored: {len(scored)}, Failed: {len(results)-len(scored)}")
    print(f"  Time: {elapsed:.1f}s ({len(scored)/elapsed*3600:.0f} lig/h)")
    if scored:
        print(f"  Best: {scored[0]['best_score']} kcal/mol ({scored[0]['smiles'][:50]})")
        print(f"  Worst: {scored[-1]['best_score']} kcal/mol")

    return scored
```

### 4.3 Batch docking (robust): checkpoint/resume + error taxonomy + env self-check

Batch docking fails in the real world: toolchains break, long runs get interrupted, and some ligands are simply not dockable.
This section provides a **production-grade** batch docking pattern with:

- **Checkpoint/resume**: restart without losing progress
- **Partial results**: append/update a CSV after each ligand
- **Standard error classification**: consistent failure accounting
- **Environment self-check**: fail fast if required binaries are missing

#### Error taxonomy (standardized)

| code | meaning |
|---|---|
| `EMBED_FAIL` | RDKit 3D embedding failed |
| `PREP_FAIL` | ligand preparation to PDBQT failed (meeko/obabel) |
| `VINA_FAIL` | Vina internal error (e.g., `tree.h`), or docking failed |
| `TIMEOUT` | per-ligand docking exceeded time limit |

All failures should be written to `errors.json` with the code + reason.

#### Robust batch docking template (copy/paste)

```python
#!/opt/conda/envs/chem/bin/python
import json
import os
import shutil
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def env_check(required_bins=("obabel", "vina")):
    missing = [b for b in required_bins if shutil.which(b) is None]
    if missing:
        raise RuntimeError(
            "Missing required binaries: " + ", ".join(missing)
            + ". Install them (or ensure they are on PATH) before batch docking."
        )


def classify_exception(msg: str) -> str:
    m = (msg or "").lower()
    if "embed" in m or "etkdg" in m:
        return "EMBED_FAIL"
    if "meeko" in m or "obabel" in m or "pdbqt" in m or "prepare" in m:
        return "PREP_FAIL"
    if "timeout" in m:
        return "TIMEOUT"
    return "VINA_FAIL"


def rdkit_embed_mmff(smiles: str, seed: int = 0xC0FFEE, max_iters: int = 500) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid_smiles")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = int(seed)
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        raise RuntimeError("embed_failed")
    res = AllChem.MMFFOptimizeMolecule(mol, maxIters=int(max_iters))
    if res != 0:
        raise RuntimeError(f"mmff_not_converged(res={res})")
    return mol


@dataclass
class DockConfig:
    receptor_pdbqt: str
    box_center: tuple[float, float, float]
    box_size: tuple[float, float, float]
    exhaustiveness: int = 16
    n_poses: int = 3
    timeout_s: int = 300


def dock_one_with_vina_cli(ligand_pdbqt: Path, out_pdbqt: Path, cfg: DockConfig) -> float:
    # Uses vina CLI for robust timeout control.
    # You can also use Python API, but CLI makes TIMEOUT enforcement simpler.
    cmd = [
        "vina",
        "--receptor", cfg.receptor_pdbqt,
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cfg.box_center[0]),
        "--center_y", str(cfg.box_center[1]),
        "--center_z", str(cfg.box_center[2]),
        "--size_x", str(cfg.box_size[0]),
        "--size_y", str(cfg.box_size[1]),
        "--size_z", str(cfg.box_size[2]),
        "--exhaustiveness", str(cfg.exhaustiveness),
        "--num_modes", str(cfg.n_poses),
        "--out", str(out_pdbqt),
    ]

    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=int(cfg.timeout_s))
    except subprocess.TimeoutExpired:
        raise TimeoutError("vina_timeout")

    if p.returncode != 0:
        raise RuntimeError((p.stderr or p.stdout or "vina_failed")[:400])

    # Parse best score from stdout (vina prints a table).
    # As a fallback, return NaN and keep the docked pose file.
    best = None
    for line in (p.stdout or "").splitlines():
        if line.strip().startswith("1 "):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    best = float(parts[1])
                except Exception:
                    best = None
            break
    if best is None:
        raise RuntimeError("vina_score_parse_failed")
    return float(best)


def load_done_set(partial_csv: Path) -> set[str]:
    if not partial_csv.exists():
        return set()
    df = pd.read_csv(partial_csv)
    if "smiles" not in df.columns:
        return set()
    # done = rows with a numeric vina_score
    done = set(df[df["vina_score"].notna()]["smiles"].astype(str).tolist())
    return done


def append_row_atomic(partial_csv: Path, row: dict):
    partial_csv.parent.mkdir(parents=True, exist_ok=True)
    df_new = pd.DataFrame([row])
    if partial_csv.exists():
        df_old = pd.read_csv(partial_csv)
        df = pd.concat([df_old, df_new], ignore_index=True)
    else:
        df = df_new
    # Best-effort de-duplication: keep the first successful score per SMILES.
    if "vina_score" in df.columns:
        df = df.sort_values(by=["smiles", "vina_score"], ascending=[True, True])
        df = df.drop_duplicates(subset=["smiles"], keep="first")
    tmp = partial_csv.with_suffix(partial_csv.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, partial_csv)


def save_errors(errors_json: Path, err_obj: dict):
    errors_json.parent.mkdir(parents=True, exist_ok=True)
    data = []
    if errors_json.exists():
        try:
            data = json.loads(errors_json.read_text(encoding="utf-8"))
        except Exception:
            data = []
    data.append(err_obj)
    tmp = errors_json.with_suffix(errors_json.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, errors_json)


def robust_virtual_screen(
    smiles_list: list[str],
    out_csv: Path,
    errors_json: Path,
    cfg: DockConfig,
    workdir: Path,
):
    env_check(required_bins=("obabel", "vina"))

    done = load_done_set(out_csv)

    for i, smi in enumerate(smiles_list, 1):
        if smi in done:
            continue

        t0 = time.time()
        row = {
            "idx": i,
            "smiles": smi,
            "vina_score": None,
            "error_code": "",
            "error": "",
            "walltime_s": None,
        }

        try:
            # (1) RDKit 3D embed (pre-flight; useful for catching geometry issues early)
            _ = rdkit_embed_mmff(smi)

            # (2) Prepare ligand pdbqt — project-specific; plug in your own function.
            # Here we assume ligand_pdbqt already exists or is created externally.
            ligand_pdbqt = workdir / f"lig_{i:04d}.pdbqt"
            if not ligand_pdbqt.exists():
                raise RuntimeError("ligand_pdbqt_missing (implement preparation step)")

            # (3) Dock with Vina CLI (timeout enforced)
            out_pdbqt = workdir / f"lig_{i:04d}_dock.pdbqt"
            score = dock_one_with_vina_cli(ligand_pdbqt, out_pdbqt, cfg)
            row["vina_score"] = float(score)

        except Exception as e:
            msg = str(e)
            row["error_code"] = classify_exception(msg)
            row["error"] = msg[:400]
            save_errors(errors_json, {
                "smiles": smi,
                "code": row["error_code"],
                "error": row["error"],
            })

        row["walltime_s"] = round(time.time() - t0, 3)
        append_row_atomic(out_csv, row)

    return pd.read_csv(out_csv)
```

**How to use it:**
- Provide `smiles_list` and a ligand preparation step that writes a `lig_XXXX.pdbqt` per molecule.
- Run once; if interrupted, rerun and it will skip completed SMILES.

---

### 4.4 Multi-seed docking robustness (pose stability + hinge robustness)

Single-run docking can be a **silent failure**: different random seeds can produce different top poses, and hinge H-bonds may only appear in a subset of runs.

Use this on a **shortlist** (typically Top50–Top100 after gates like safety / basic plausibility) to quantify robustness.

#### Protocol (3–5 seeds per molecule)

1) For each molecule, run docking **K times** with different seeds (K=3–5).
   - If using vina CLI, pass `--seed <int>`.
   - Keep all other settings identical.
2) Collect the **top pose** from each run.
3) RMSD-cluster these top poses (within-molecule) to measure pose convergence.
4) Compute **hinge H-bond consistency** across seeds (requires ProLIF hinge check per pose, or a lightweight hinge detector).

#### Gates (default)

- `pose_convergence`:
  - at least **3/5 runs** have top-pose RMSD < **2.0 Å** to the dominant cluster representative.
- `hinge_robustness`:
  - hinge H-bond appears in **≥ 3/5 runs**.

Interpretation:
- Passing both gates means the molecule’s docking pose and hinge evidence are **not seed-fragile**.
- Failing either gate means: treat as lower confidence; consider re-docking with higher exhaustiveness, alternative protonation state, or deprioritize.

#### Output (per-molecule robustness report)

Write a standardized CSV/TSV with one row per molecule:

- `mol_id`
- `n_seeds`
- `n_converged` (count of runs in dominant cluster within 2.0 Å)
- `hinge_consistency` (fraction in [0,1])
- `robust` (bool; `pose_convergence && hinge_robustness`)

Example schema:

```text
mol_id,n_seeds,n_converged,hinge_consistency,robust
cycle5_hinge_1,5,4,0.8,True
```

#### Practical note on RMSD

- RMSD clustering requires **consistent atom mapping** across poses.
- If your pipeline outputs SDF poses with consistent atom order, you can compute RMSD by RDKit alignment (e.g., `rdMolAlign.AlignMol`).
- If you only have PDBQT, convert to SDF (or RDKit Mol) in a way that preserves atom correspondence.

## Phase 5: Score Interpretation & Analysis

### 5.1 Score Guidelines

```
Vina score interpretation (approximate, target-dependent):
├── < -10.0 kcal/mol: Excellent — strong predicted binder (rare)
├── -8.0 to -10.0: Good — worth pursuing
├── -6.0 to -8.0: Moderate — common for drug-like molecules
├── -4.0 to -6.0: Weak — borderline, likely not active
└── > -4.0: Very weak — unlikely to bind meaningfully

CRITICAL: These are guidelines, not rules. Always:
  1. Compare to a KNOWN ACTIVE ligand docked in the same setup
  2. If known active scores -7.5, candidates < -7.0 are promising
  3. Score differences < 1.0 kcal/mol are within noise
```

### 5.2 Redocking Validation

```python
def validate_docking_setup(receptor_pdbqt, known_ligand_smiles, known_ligand_pdb, box):
    """Validate docking setup by redocking the co-crystallized ligand.

    Success criterion: RMSD < 2.0 Å between docked and crystal pose.
    This is the FIRST thing to do before any virtual screening.
    """
    # Prepare known ligand from SMILES
    known_pdbqt = prepare_ligand(known_ligand_smiles, "known_ligand.pdbqt")
    if not known_pdbqt:
        return {"valid": False, "error": "ligand preparation failed"}

    # Dock
    results, docked_path = dock_single(receptor_pdbqt, known_pdbqt, box,
                                        exhaustiveness=32, n_poses=9)

    # Parse coordinates of docked pose 1
    docked_coords = parse_pdbqt_coords(docked_path, pose=1)
    crystal_coords = parse_pdb_coords(known_ligand_pdb)

    if docked_coords is not None and crystal_coords is not None:
        # Heavy atom RMSD (approximate — atom matching not trivial)
        n = min(len(docked_coords), len(crystal_coords))
        rmsd = np.sqrt(np.mean(np.sum((docked_coords[:n] - crystal_coords[:n])**2, axis=1)))
    else:
        rmsd = None

    validation = {
        "known_ligand_score": results[0]["score_kcal_mol"],
        "rmsd_A": round(rmsd, 2) if rmsd else None,
        "valid": rmsd is not None and rmsd < 2.0,
    }

    print(f"\nRedocking validation:")
    print(f"  Known ligand score: {validation['known_ligand_score']} kcal/mol")
    print(f"  RMSD to crystal: {validation['rmsd_A']} Å")
    print(f"  Valid setup: {'✓' if validation['valid'] else '✗ (RMSD > 2.0 Å — check box/receptor)'}")

    return validation

def parse_pdbqt_coords(pdbqt_path, pose=1):
    """Parse heavy atom coordinates from a PDBQT file."""
    coords = []
    current_pose = 0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_pose += 1
            if current_pose > pose:
                break
            if current_pose == pose and line[:4] in ("ATOM", "HETA"):
                element = line[77:79].strip() if len(line) > 77 else line[12:14].strip()
                if element != "H":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    return np.array(coords) if coords else None

def parse_pdb_coords(pdb_path):
    """Parse heavy atom coordinates from a PDB file."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM"):
                element = line[76:78].strip() if len(line) > 76 else line[12:14].strip()
                if element != "H":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    return np.array(coords) if coords else None
```

### 5.3 VS Results Analysis

```python
import pandas as pd

def analyze_vs_results(results, known_active_score=None):
    """Analyze virtual screening results."""
    df = pd.DataFrame([r for r in results if r.get("best_score") is not None])
    if df.empty:
        print("No valid results to analyze")
        return df

    print(f"\n=== Virtual Screening Analysis ===")
    print(f"  N scored: {len(df)}")
    print(f"  Score range: [{df['best_score'].min():.2f}, {df['best_score'].max():.2f}] kcal/mol")
    print(f"  Mean ± std: {df['best_score'].mean():.2f} ± {df['best_score'].std():.2f}")
    print(f"  Median: {df['best_score'].median():.2f}")

    if known_active_score:
        n_better = (df["best_score"] < known_active_score).sum()
        print(f"  Better than known active ({known_active_score}): {n_better} ({n_better/len(df):.1%})")

    # Score distribution bins
    bins = [(-999, -10), (-10, -8), (-8, -6), (-6, -4), (-4, 999)]
    labels = ["< -10", "-10 to -8", "-8 to -6", "-6 to -4", "> -4"]
    for (lo, hi), label in zip(bins, labels):
        n = ((df["best_score"] >= lo) & (df["best_score"] < hi)).sum()
        print(f"  {label} kcal/mol: {n}")

    return df
```

## Phase 6: Integration with Pipeline

### 6.1 Post-Generation Docking Filter (connects to skill 3 + 6)

```python
def dock_and_filter_candidates(receptor_pdbqt, box, candidate_smiles,
                                score_cutoff=-7.0, admet_filter=True):
    """Full pipeline: prepare → dock → score → ADMET filter → rank.

    This is the main integration point between generative models (skill 3),
    ADMET filtering (skill 6), and docking.
    """
    print(f"=== Docking + Filtering Pipeline ===")
    print(f"  Candidates: {len(candidate_smiles)}")
    print(f"  Score cutoff: {score_cutoff} kcal/mol")

    # Step 1: Prepare ligands
    prepared, failed_prep = prepare_ligands_batch(candidate_smiles)
    print(f"  Prepared: {len(prepared)}, Failed: {len(failed_prep)}")

    # Step 2: Virtual screening
    vs_results = virtual_screen(receptor_pdbqt, prepared, box)

    # Step 3: Score filter
    hits = [r for r in vs_results if r["best_score"] is not None
            and r["best_score"] <= score_cutoff]
    print(f"  Score hits (≤ {score_cutoff}): {len(hits)}")

    # Step 4: ADMET filter (if requested, uses skill 6)
    if admet_filter:
        try:
            # Import from skill 6 code or inline
            from rdkit.Chem import Descriptors, QED

            admet_passed = []
            for hit in hits:
                mol = Chem.MolFromSmiles(hit["smiles"])
                if mol is None:
                    continue
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                qed = QED.qed(mol)

                # Lipinski
                violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
                if violations < 2 and qed > 0.3:
                    hit["mw"] = round(mw, 1)
                    hit["logp"] = round(logp, 2)
                    hit["qed"] = round(qed, 3)
                    hit["lipinski_violations"] = violations
                    admet_passed.append(hit)

            print(f"  ADMET passed: {len(admet_passed)}")
            hits = admet_passed
        except Exception as e:
            print(f"  ADMET filter skipped: {e}")

    # Sort by score
    hits.sort(key=lambda x: x["best_score"])
    print(f"  Final hits: {len(hits)}")

    return hits
```

### 6.2 Full Pipeline Results Table

Always produce this:

```markdown
| Rank | SMILES | Vina Score | QED | MW | LogP | SA Score | Route Steps | Notes |
|------|--------|-----------|-----|-----|------|----------|-------------|-------|
| 1 | ... | -8.5 | 0.82 | 332 | 2.3 | 2.2 | 1 | Top hit |
| 2 | ... | -8.2 | 0.65 | 380 | 3.1 | 4.7 | 2 | |
| 3 | ... | -7.9 | 0.71 | 295 | 1.8 | 3.1 | 1 | Lead-like |
```

### 6.3 Decision Tree: When Docking Adds Value

```
Is docking worth running?
├── You have a crystal structure with co-crystallized ligand
│   ├── Redocking RMSD < 2.0 Å → YES, docking is reliable here
│   └── Redocking RMSD > 2.0 Å → FIX setup first (box? receptor prep?)
├── You have a homology model (no crystal ligand)
│   └── Use with caution. Scores less reliable. Rank only, don't filter.
├── You have > 1000 candidates to screen
│   └── YES — docking is faster than synthesis. Use exhaustiveness=8-16.
├── You have < 10 candidates
│   └── MAYBE — better to just use ADMET + SA + QSAR scores
└── Your target is a GPCR / membrane protein
    └── Extra care needed. Ensure proper receptor preparation.
```

## Checklist Before Reporting

- [ ] **Redocking validated**: Known ligand RMSD < 2.0 Å?
- [ ] **Box defined from reference**: Not just guessed coordinates?
- [ ] **Scores compared to known active**: Relative ranking, not absolute cutoffs?
- [ ] **Score as rank, not energy**: Not claiming Kd from Vina score?
- [ ] **Top poses inspected**: At least described key interactions for top 5?
- [ ] **ADMET filter applied**: Lipinski/QED on hits (skill 6)?
- [ ] **SA Score included**: Synthesizability of hits (skill 4)?
- [ ] **Exhaustiveness stated**: What search depth was used?

## Integration with Knowledge Base

- **Save docking results** to `research/ai4chem/docking/<target>/<date>-vs.md`
- **Save receptor files** to `research/ai4chem/docking/<target>/receptor/`
- **Save top hit poses** to `research/ai4chem/docking/<target>/poses/`
- **Cross-reference** with ADMET profiles (skill 6) and retro routes (skill 4)
- **Git commit**: `cd <OPENCLAW_WORKSPACE> && git add -A && git commit -m "docking: <target> VS results"`
