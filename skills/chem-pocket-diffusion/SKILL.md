---
name: chem-pocket-diffusion
description: Pocket-conditioned 3D ligand generation with diffusion models (DiffSBDD focus) plus an auditable integration path into RDKit→safety→Vina→ProLIF pipelines. Covers environment setup (PyG/equivariance), input preparation (pocket PDB + reference ligand), sampling (timesteps), substructure inpainting, and hard gates (validity).
homepage: https://github.com/arneschneuing/DiffSBDD
metadata: { "openclaw": { "emoji": "🧩", "requires": { "bins": ["python3", "obabel"], "python": ["torch", "pytorch_lightning", "rdkit", "openbabel", "biopython", "biopandas", "pandas", "numpy" ] } } }
---

# Pocket-Conditioned Molecular Diffusion (DiffSBDD) — Practical Playbook

This skill operationalizes **pocket-conditioned 3D ligand generation** using **DiffSBDD** and integrates it into a standard, auditable SBDD loop:

> generate (DiffSBDD) → sanitize (RDKit) → safety gates → docking (Vina) → interaction KPI (ProLIF)

It is written from hands-on deployment experience: dependency pitfalls, required patches, and failure modes.

---

## When to Use

- You have a protein structure and want **de novo ligands** conditioned on the pocket.
- You need **substructure inpainting** to enforce a hinge fragment or a scaffold.
- You want an end-to-end 3D generator but still need downstream **docking + interaction validation**.

## When NOT to Use

- You require guaranteed motif enforcement without post hoc gates (diffusion is probabilistic).
- You cannot tolerate GPU dependencies (DiffSBDD’s practical throughput assumes GPU).

---

## Core philosophy

1. **Diffusion generates candidates, not truth.** You still need RDKit sanitization and docking-based validation.
2. **Hard gates prevent self-deception.** If validity is low, do not proceed to docking.
3. **Inpainting is the right way to “force” chemistry.** Prefer substructure inpainting over soft biases when a fragment is non-negotiable.

---

## Phase 0 — Environment setup (the part that breaks)

DiffSBDD sits on a geometric DL stack. Expect most failures here.

### 0.1 Minimal dependency checklist

You need:
- `torch` with CUDA
- `pytorch-lightning`
- `wandb` (can be disabled, but the import often exists)
- `rdkit`
- `openbabel` (`obabel` CLI)
- `biopandas`
- `biopython`
- PyG stack (`torch-geometric`, `torch-scatter`, etc.)

### 0.2 Install hints (practical)

Use your existing Python (example path):

```bash
PY=/opt/conda/envs/chem/bin/python
PIP=/opt/conda/envs/chem/bin/pip

$PY -c "import torch; print(torch.__version__, torch.cuda.is_available())"

# common python deps
$PIP install -U pytorch-lightning wandb biopandas biopython pyyaml tqdm pandas numpy scipy
```

#### PyTorch 2.6 `weights_only` patch (checkpoint load failure)

PyTorch 2.6 tightened defaults around `torch.load` and may refuse to unpickle some objects (e.g., `argparse.Namespace`) unless explicitly allow-listed.

**Patch (confirmed working):** add the following line at the **top of `generate_ligands.py`** (first line after the shebang, before other imports):

```python
import torch, argparse; torch.serialization.add_safe_globals([argparse.Namespace])
```

If the checkpoint still fails to load, fall back to explicitly disabling `weights_only` in the relevant `torch.load(...)` call:

```python
torch.load(path, map_location=device, weights_only=False)
```

#### BioPython `three_to_one` patch (import/API drift)

Some Biopython versions moved/changed amino-acid utilities (`three_to_one`).

**Patch (confirmed working):** replace:

```python
from Bio.PDB.Polypeptide import three_to_one
```

with:

```python
from Bio.Data.IUPACData import protein_letters_3to1 as _map
three_to_one = lambda x: _map[x.lower()]
```

(Do not silently drop residues: if the pocket becomes empty you will generate garbage.)

### 0.3 Sanity check

```bash
$PY -c "import rdkit; import openbabel; import Bio; import biopandas; import pytorch_lightning; print('deps_ok')"
```

---

## Phase 1 — Inputs (pocket PDB + reference ligand SDF)

DiffSBDD needs a protein structure and a pocket definition.

### 1.1 Protein input (PDB)

- Provide a **PDB** file for the target.
- Ensure chain IDs and residue IDs are stable.
- Hydrogens are optional for generation, but consistency helps downstream steps.

### 1.2 Reference ligand SDF (from PDBQT)

DiffSBDD supports pocket definition via a reference ligand.

If you have a reference ligand as **PDBQT** (e.g., from docking/redocking), convert it to SDF:

```bash
obabel ref_ligand.pdbqt -O ref_ligand.sdf
```

This SDF is used only to define the pocket region; you will still dock the generated ligands later.

---

## Phase 2 — Generation (DiffSBDD)

DiffSBDD exposes a simple generation interface:

```bash
python generate_ligands.py <checkpoint>.ckpt \
  --pdbfile <protein>.pdb \
  --outfile <out>.sdf \
  --ref_ligand <ref>.sdf \
  --n_samples 10 \
  --timesteps 100 \
  --sanitize
```

### 2.1 Walltime benchmark (reference)

Empirical reference (for feasibility planning):
- On **RTX 4080 Laptop (12GB VRAM)**, **10 samples** at `timesteps=100` completed in **~30 seconds**.

(Expect walltime to vary with checkpoint, pocket size, and chosen timesteps.)

### 2.2 Key parameters

- `--n_samples`: number of molecules to generate.
- `--timesteps`: number of denoising steps.
  - **Lower** timesteps → faster, lower quality.
  - **Higher** timesteps → slower, better sample refinement.

**Rule of thumb:**
- feasibility test: `timesteps=50–100`
- higher quality sampling: `timesteps=200+`

### 2.3 Outputs

- Primary output is typically an **SDF** file with 3D coordinates.
- You should extract SMILES and compute validity/uniqueness as immediate QC.

---

## Phase 3 — Post-processing (RDKit)

### 3.1 RDKit sanitization

Even with `--sanitize`, expect failures. A reasonable expectation for 3D generators is ~**80%** RDKit-valid molecules, but this is model/checkpoint dependent.

Extract canonical SMILES from SDF:

```python
from rdkit import Chem

suppl = Chem.SDMolSupplier("generated.sdf", removeHs=False)
smiles=[]
for m in suppl:
    if m is None:
        continue
    smiles.append(Chem.MolToSmiles(m))

n_valid = len(smiles)
n_unique = len(set(smiles))
```

### 3.2 Gate: validity

Validity gates for this skill:
- **validity ≥ 70%** → **PASS** (proceed)
- **60% ≤ validity < 70%** → **WARNING** (proceed only for debugging; expect downstream noise)
- **validity < 60%** → **FAIL** (do not proceed)

If the gate fails:
- reduce timesteps only if you’re debugging speed; otherwise increase timesteps
- verify pocket definition is correct (wrong pocket → nonsense ligands)
- ensure checkpoint loads correctly (no partial-load / wrong device)

---

## Phase 4 — Integrate with our existing pipeline

DiffSBDD output is *not* the final candidate set. Integrate it into the same pipeline used for string generation:

1) **Safety screen** (project denylist / chronic policy)
2) **Docking** with Vina using a consistent box
3) **Interaction validation** (ProLIF): hinge H-bond KPI + interaction rerank

This keeps results comparable across generators.

---

## Phase 5 — Substructure inpainting (fragment conditioning)

DiffSBDD supports **inpainting** to fix parts of a ligand and generate the rest.

> Based on DiffSBDD repo documentation; **not yet validated in our pipeline**. Verify before relying on it.

Typical usage pattern:

```bash
python inpaint.py <checkpoint>.ckpt \
  --pdbfile <protein>.pdb \
  --outfile <out>.sdf \
  --ref_ligand <ref>.sdf \
  --fix_atoms <atom_ids...> \
  --add_n_nodes <k>
```

How to use this for kinase hinge enforcement:
- choose a hinge-binding fragment (from a known binder)
- provide a reference ligand containing that fragment
- fix atoms corresponding to the fragment
- generate the remaining atoms around it

**Expectation:** inpainting increases the probability of retaining the fragment, but still requires docking + hinge H-bond verification.

---

## Failure modes & workarounds

### A) Checkpoint load fails (PyTorch 2.6)
- Symptom: `torch.load` errors about `weights_only`
- Workaround: patch loader to call `weights_only=False`.

### B) Pocket is empty / no atoms selected
- Symptom: generation produces garbage or crashes
- Workaround: verify pocket definition (ref ligand SDF is correct; PDB chain/residue mapping is correct).

### C) Biopython residue conversion errors
- Symptom: import error or key error in residue mapping
- Workaround: implement a local residue mapping fallback (e.g., map `MSE→MET`) and use a supported `three_to_one` source for your Biopython version.

### D) Low RDKit validity
- Symptom: many `None` molecules in SDMolSupplier
- Workarounds:
  - increase timesteps
  - ensure `--sanitize` is used
  - run additional RDKit cleanup (remove fragments, neutralize, canonical tautomer)

### E) SDF contains multiple fragments (`.` in SMILES)
- Workaround: keep largest fragment (RDKit LargestFragmentChooser) before safety/docking.

### F) Downstream mismatch: good diffusion samples but no hinge H-bond
- Workaround: use inpainting with a hinge fragment; enforce hinge KPI as a hard gate after docking.

---

## Deliverables (what this skill should produce)

Minimum artifacts per run:
- Generated SDF (3D)
- Extracted SMILES list + validity/uniqueness summary
- Safety survivors
- Docking scores
- Interaction Top5 report with hinge H-bond evidence

This makes the workflow auditable and comparable across cycles and generators.
