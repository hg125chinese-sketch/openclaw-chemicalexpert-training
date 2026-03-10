---
name: chem-affinity-prediction
description: Protein-ligand binding affinity prediction using Boltz-2. Covers environment isolation (mandatory separate conda env), input preparation (YAML with protein sequence + ligand SMILES + affinity property), inference, output interpretation (IC50 + binder probability), and integration with Vina/ProLIF as an orthogonal scoring signal.
homepage: https://github.com/jwohlwend/boltz
metadata: { "openclaw": { "emoji": "🎯", "requires": { "bins": ["boltz"], "python": ["boltz", "torch", "pytorch_lightning", "rdkit"] } } }
---

# Binding Affinity Prediction (Boltz-2) — Practical Playbook

This skill adds **protein-ligand binding affinity prediction** to the pipeline using Boltz-2, providing IC50 estimates and binder probability as orthogonal scoring signals alongside Vina docking and ProLIF interaction analysis.

> DiffSBDD generates → safety screen → Vina docking → ProLIF interaction → **Boltz-2 affinity** → QE DFT

Written from hands-on deployment experience: environment isolation pitfalls, YAML format discovery, and real results on IPF/ALK5 candidates.

---

## When to Use

- You have DFT-validated candidates and want a second opinion on binding strength beyond Vina scores.
- You want to rank compounds by predicted IC50 rather than docking score.
- You need binder/decoy classification (hit discovery stage).
- You want to compare candidates from different generation methods (e.g., DiffSBDD vs VAE).

## When NOT to Use

- You need binding free energy with kcal/mol precision (use FEP instead).
- Your ligand has >50 heavy atoms (poorly represented in Boltz-2 training data).
- You need protein-protein binding affinity (not yet supported by Boltz-2).

---

## Core Philosophy

1. **Boltz-2 is an orthogonal signal, not a replacement.** Use it alongside Vina and ProLIF, not instead of them. Boltz-2 ranking can disagree with Vina — this is informative, not a bug.
2. **Environment isolation is mandatory.** Boltz-2 requires torch 2.10+ which is incompatible with our chem environment (torch 2.6). Never install Boltz in the chem conda env.
3. **Two outputs for two use cases.** `affinity_probability_binary` for hit discovery (is it a binder?), `affinity_pred_value` for lead optimization (how strong?).

---

## Phase 0 — Environment Setup (CRITICAL: Separate Conda Env)

### 0.1 Why a Separate Environment

Boltz-2 depends on torch 2.10+ and numpy <2.0. Installing it in the chem environment will:
- Upgrade torch from 2.6 to 2.10, breaking RDKit (GLIBCXX version mismatch)
- Downgrade numpy from 2.2 to 1.26, breaking scipy and other packages
- Corrupt the entire chem pipeline

**This was verified the hard way.** Do not repeat this mistake.

### 0.2 Create Isolated Environment

```bash
# Create a fresh conda env (one-time)
/opt/conda/bin/conda create -n boltz python=3.11 -y

# Install Boltz with CUDA support
PYTHONNOUSERSITE=1 /home/node/.conda/envs/boltz/bin/pip install "boltz[cuda]"
```

Note: The environment may be created at `/home/node/.conda/envs/boltz/` if `/opt/conda/envs/` is not writable.

### 0.3 PYTHONNOUSERSITE=1 is Required

The boltz environment shares Python 3.11 with chem. Without `PYTHONNOUSERSITE=1`, it will import packages from `~/.local/lib/python3.11/site-packages/` (chem's user packages), causing version conflicts.

**Every boltz invocation must use:**

```bash
PYTHONNOUSERSITE=1 /home/node/.conda/envs/boltz/bin/boltz predict ...
```

### 0.4 Sanity Check

```bash
PYTHONNOUSERSITE=1 /home/node/.conda/envs/boltz/bin/boltz predict --help
```

Should print usage without import errors.

### 0.5 First Run Downloads (~2GB)

On first run, Boltz-2 downloads:
- CCD data (`~/.boltz/mols.tar`)
- Structure model weights (`~/.boltz/boltz2.ckpt`)
- Affinity model weights (`~/.boltz/boltz2_aff.ckpt`)

This is a one-time cost.

---

## Phase 1 — Input Preparation (YAML Format)

### 1.1 YAML Structure

Boltz-2 takes a YAML file per prediction. Minimum required for affinity:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: <PROTEIN_SEQUENCE>
  - ligand:
      id: B
      smiles: "<SMILES>"
properties:
  - affinity:
      binder: B
```

**Critical details:**
- Entity type for small molecules is `ligand`, NOT `smiles`
- The `properties` section with `affinity.binder` is **required** for affinity prediction — without it, only structure prediction runs
- `binder` value must match the ligand's `id`
- SMILES must be quoted to avoid YAML parsing issues with special characters

### 1.2 Protein Sequence

Use the full kinase domain sequence. Truncated or incomplete sequences produce low confidence scores (we observed confidence ~0.33 with a partial ALK5 sequence). For ALK5/TGFBR1, use the sequence from UniProt P36897.

### 1.3 Batch Preparation

For multiple ligands, create one YAML per ligand in a directory:

```python
import os

def write_boltz_yamls(protein_seq, smiles_dict, output_dir):
    """Write one YAML per ligand for batch Boltz-2 prediction.
    
    smiles_dict: {mol_id: smiles_string}
    """
    os.makedirs(output_dir, exist_ok=True)
    for mol_id, smi in smiles_dict.items():
        yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {protein_seq}
  - ligand:
      id: B
      smiles: "{smi}"
properties:
  - affinity:
      binder: B
"""
        with open(f"{output_dir}/{mol_id}.yaml", "w") as f:
            f.write(yaml_content)
```

Then pass the directory to `boltz predict`.

---

## Phase 2 — Running Predictions

### 2.1 Single Prediction

```bash
PYTHONNOUSERSITE=1 /home/node/.conda/envs/boltz/bin/boltz predict \
  input.yaml \
  --out_dir /path/to/output \
  --use_msa_server \
  --num_workers 0
```

### 2.2 Batch Prediction

```bash
PYTHONNOUSERSITE=1 /home/node/.conda/envs/boltz/bin/boltz predict \
  /path/to/yaml_dir/ \
  --out_dir /path/to/output \
  --use_msa_server \
  --num_workers 0
```

**Important flags:**
- `--use_msa_server`: Required for protein inputs (generates MSA via MMseqs2 server). Without this, prediction fails with "Missing MSA" error.
- `--num_workers 0`: Prevents shared memory allocation failures in Docker containers. Without this, batch jobs may hang at 25%.
- `--override`: Re-run predictions that already exist in the output directory.

### 2.3 Walltime Benchmark

On RTX 4080 Laptop (12GB VRAM):
- Single prediction (structure + affinity): ~35 seconds (MSA ~20s + structure ~12s + affinity ~3s)
- Batch of 4: ~4 minutes total (MSA cached after first protein)

### 2.4 Output Structure

```
out_dir/
└── boltz_results_<dir_name>/
    └── predictions/
        └── <mol_id>/
            ├── <mol_id>_model_0.cif          # 3D structure
            ├── confidence_<mol_id>_model_0.json  # Structure confidence
            └── affinity_<mol_id>.json         # Binding affinity
```

---

## Phase 3 — Interpreting Results

### 3.1 Affinity Output Fields

```json
{
    "affinity_pred_value": -0.62,
    "affinity_probability_binary": 0.43,
    "affinity_pred_value1": -0.93,
    "affinity_probability_binary1": 0.52,
    "affinity_pred_value2": -0.32,
    "affinity_probability_binary2": 0.35
}
```

- **affinity_pred_value**: Mean predicted log10(IC50) in µM. More negative = stronger binding.
  - Convert to IC50: `IC50_uM = 10^(affinity_pred_value)`
  - Example: -0.62 → IC50 ≈ 0.24 µM (submicromolar)
- **affinity_probability_binary**: Mean probability the ligand is a binder (0-1).
  - Use for hit discovery / binder vs decoy classification.
- **Numbered variants** (value1/2, probability1/2): Individual diffusion samples. Spread indicates prediction uncertainty.

### 3.2 Use Cases

| Stage | Use this field | Threshold suggestion |
|-------|---------------|---------------------|
| Hit discovery | affinity_probability_binary | > 0.5 = likely binder |
| Lead optimization | affinity_pred_value | More negative = better; compare relative ranking |

### 3.3 Confidence Context

Structure confidence (from confidence JSON) affects affinity reliability:
- **confidence_score > 0.5**: Affinity prediction is more trustworthy
- **confidence_score < 0.3**: Take affinity with a grain of salt — structure may be unreliable

---

## Phase 4 — Integration with Existing Pipeline

### 4.1 Boltz-2 as Orthogonal Signal

Boltz-2 affinity should be used **alongside** Vina and ProLIF, not as a replacement:

```
DiffSBDD → safety → Vina docking → ProLIF interaction → Boltz-2 affinity → QE DFT
```

### 4.2 Why Rankings Disagree

In our IPF/ALK5 validation, Boltz-2 and Vina produced different rankings:

| Molecule | Vina Score | Boltz-2 IC50 |
|----------|-----------|--------------|
| hinge_1  | -10.01 (best Vina) | 0.47 µM (worst Boltz) |
| hinge_3  | -8.94 (worst Vina) | 0.23 µM (best Boltz) |

This is expected — Vina scores approximate binding free energy from a rigid docking pose, while Boltz-2 predicts IC50 from a co-folded structure with learned affinity patterns. Disagreement is diagnostic information, not error.

### 4.3 Multi-Score Integration

When combining Vina, ProLIF, and Boltz-2, use a panel approach:

1. **Hard gates**: hinge H-bond (ProLIF) = required, MACE prescreen = required
2. **Ranking signals**: Vina score, Boltz-2 affinity_pred_value, Boltz-2 binder probability
3. **Selection**: Prefer candidates that rank well on at least 2 of 3 signals

Do NOT combine into a single composite score without validation — the signals have different scales and calibrations.

---

## Phase 5 — Parsing Results Programmatically

```python
import json
import os
import glob

def parse_boltz_affinity(results_dir):
    """Parse all affinity results from a Boltz-2 batch run."""
    results = []
    pattern = os.path.join(results_dir, "predictions/*/affinity_*.json")
    for path in sorted(glob.glob(pattern)):
        mol_id = os.path.basename(os.path.dirname(path))
        with open(path) as f:
            data = json.load(f)
        ic50_um = 10 ** data["affinity_pred_value"]
        results.append({
            "mol_id": mol_id,
            "affinity_log_ic50": data["affinity_pred_value"],
            "ic50_uM": round(ic50_um, 3),
            "binder_prob": data["affinity_probability_binary"],
        })
    # Sort by affinity (most negative = strongest)
    results.sort(key=lambda x: x["affinity_log_ic50"])
    return results
```

---

## Failure Modes & Workarounds

### A) Import errors / GLIBCXX mismatch
- **Symptom**: `ImportError: libstdc++.so.6: version GLIBCXX_3.4.31 not found`
- **Cause**: Boltz was installed in the chem conda env, upgrading torch
- **Fix**: Uninstall boltz from chem env, restore torch 2.6, use separate boltz env

### B) "Invalid entity type: smiles"
- **Symptom**: YAML parsing error
- **Cause**: Using `smiles:` instead of `ligand:` as entity type
- **Fix**: Use `ligand:` with `smiles:` as a sub-field

### C) "Missing MSA's in input"
- **Symptom**: RuntimeError about missing MSA
- **Fix**: Add `--use_msa_server` flag

### D) No affinity JSON in output
- **Symptom**: Structure prediction runs but no affinity_*.json
- **Cause**: Missing `properties: - affinity: binder: <id>` in YAML
- **Fix**: Add the properties section to YAML

### E) Batch hangs at 25%
- **Symptom**: DataLoader stops progressing, shared memory errors in log
- **Cause**: Docker container shared memory limit
- **Fix**: Use `--num_workers 0`

### F) Low confidence scores (<0.3)
- **Symptom**: Poor structure confidence
- **Possible causes**: Incomplete protein sequence, unusual ligand chemistry
- **Workaround**: Use full-length domain sequence; check if ligand has >50 heavy atoms

---

## Deliverables

Minimum artifacts per Boltz-2 run:
- Affinity results table: mol_id, log10(IC50), IC50_µM, binder_probability
- Comparison with Vina ranking (highlight disagreements)
- Structure confidence scores
- Integration recommendation (which candidates rank well on multiple signals)
