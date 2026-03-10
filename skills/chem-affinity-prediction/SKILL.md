---
name: chem-affinity-prediction
description: Protein-ligand binding affinity prediction using Boltz-2. Covers environment isolation (mandatory separate conda env), input preparation (YAML with protein sequence + ligand SMILES + affinity property), inference, output interpretation (IC50 + binder probability), and integration with Vina/ProLIF as an orthogonal scoring signal.
homepage: https://github.com/jwohlwend/boltz
metadata: { "openclaw": { "emoji": "🎯", "requires": { "bins": ["boltz"], "python": ["boltz", "torch", "pytorch_lightning", "rdkit"] } } }
---

# chem-affinity-prediction — Boltz-2 (affinity) as an orthogonal signal

Docking is great for *poses* and *local interaction evidence* (hinge H-bond, key residues). But it’s also noisy and rigid.

Boltz-2 affinity gives you a **second, orthogonal ranking signal** on a shortlist (ideally **DFT PASS** molecules). The key is doing this **without destroying** the main chemistry environment.

This skill is written from hands-on workflow issues we actually hit:
- boltz install/uninstall can silently break numpy/scipy/torch stacks
- YAML formatting is surprisingly strict
- Docker multiprocessing (`num_workers`) can hang

---

## Workspace variables (use these; do not hardcode)

```bash
# Workspace
OPENCLAW_WORKSPACE="/home/node/.openclaw/workspace-chemicalexpert"

# QMD
QMD="/home/node/.openclaw/.npm-global/bin/qmd"

# Main chemistry env (DO NOT install boltz here)
CHEM_PY="/opt/conda/envs/chem/bin/python"

# Boltz env (MUST be isolated)
BOLTZ_ENV_BIN="$HOME/.conda/envs/boltz/bin"
BOLTZ="$BOLTZ_ENV_BIN/boltz"
BOLTZ_PY="$BOLTZ_ENV_BIN/python"

# Isolation (mandatory)
export PYTHONNOUSERSITE=1
```

---

## When to use (be specific)

Use Boltz-2 affinity when:

1) You have **high-confidence candidates** already:
- **DFT QC PASS** (preferred), or at least
- Vina Top5/Top10 with **hinge H-bond = True** (ProLIF-validated)

2) You need to resolve ranking ambiguity:
- Vina suggests A > B, but interaction evidence is mixed, **or**
- you suspect docking score is overfitting pose artifacts.

3) You want a hit-discovery style signal:
- use `affinity_probability_binary` as binder vs decoy probability
  - **suggested triage**: >0.5 “likely binder”, <0.3 “likely weak/decoy” (heuristic)

Do NOT use Boltz-2 when:
- You need kcal/mol-accurate ΔG (use FEP)
- Ligands are far outside drug-like space (very large, unusual chemistry)

---

## Core rules (non-negotiable)

1) **Environment isolation is mandatory.**
   - Boltz-2 depends on **torch 2.10+** and **numpy < 2.0**.
   - Installing it into the main chem env can break RDKit/scipy/torch.

2) **Always run with `PYTHONNOUSERSITE=1`.**
   - Prevents accidental imports from `~/.local/...` that can poison the boltz env.

3) **Expect disagreement with Vina/ProLIF.**
   - Treat disagreement as information.
   - Use a **panel** strategy (pick a small set that balances signals) rather than trusting a single metric.

---

## Phase 0 — Environment setup (isolated boltz conda env)

### 0.1 Create env

```bash
/opt/conda/bin/conda create -n boltz python=3.11 -y

# Install Boltz with CUDA support
PYTHONNOUSERSITE=1 "$HOME/.conda/envs/boltz/bin/pip" install "boltz[cuda]"
```

Notes:
- The env may live under `$HOME/.conda/envs/boltz/`.
- Do **not** mix this with `/opt/conda/envs/chem`.

### 0.2 Sanity check

```bash
PYTHONNOUSERSITE=1 "$HOME/.conda/envs/boltz/bin/boltz" predict --help
```

### 0.3 First-run downloads (~2GB)

Boltz-2 may download CCD + weights on first run (one-time cost).

---

## Phase 1 — YAML input preparation (the common failure point)

Boltz-2 uses one YAML per prediction. Minimum for affinity:

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

Pitfalls (real):
- small molecule entity type is `ligand` (not `smiles`)
- `properties.affinity.binder` is **required** for affinity; otherwise you only run structure
- `binder` must match ligand id (here `B`)
- **quote SMILES** (YAML parsing breaks on special characters)

### 1.1 Batch writer (script-like example)

```python
#!/opt/conda/envs/chem/bin/python
"""Write one Boltz YAML per ligand.

Run this in the chem env (fast IO), then run boltz in the boltz env.
"""

from pathlib import Path


def write_boltz_yamls(*, protein_seq: str, smiles_by_id: dict[str, str], out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    for mol_id, smi in smiles_by_id.items():
        yml = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {protein_seq}
  - ligand:
      id: B
      smiles: \"{smi}\"
properties:
  - affinity:
      binder: B
"""
        (out / f"{mol_id}.yaml").write_text(yml, encoding="utf-8")
```

Protein sequence note:
- Use the full kinase domain sequence (e.g., ALK5/TGFBR1 from UniProt P36897). Partial sequences reduced confidence in practice.

---

## Phase 2 — Run Boltz predictions (single + batch)

### 2.1 Single

```bash
PYTHONNOUSERSITE=1 "$HOME/.conda/envs/boltz/bin/boltz" predict \
  input.yaml \
  --out_dir <OUT_DIR> \
  --use_msa_server \
  --num_workers 0
```

### 2.2 Batch

```bash
PYTHONNOUSERSITE=1 "$HOME/.conda/envs/boltz/bin/boltz" predict \
  <YAML_DIR>/ \
  --out_dir <OUT_DIR> \
  --use_msa_server \
  --num_workers 0 \
  --override
```

Flags that matter:
- `--use_msa_server`: required for protein inputs (otherwise “Missing MSA”)
- `--num_workers 0`: avoids Docker shared-memory multiprocessing issues (we saw hangs around ~25%)
- `--override`: rerun existing outputs

### 2.3 Walltime benchmark (reference)

On RTX 4080 Laptop (12GB VRAM):
- single (structure + affinity): ~35s (MSA ~20s + structure ~12s + affinity ~3s)
- batch of 4: ~4 min total (MSA cached after first protein)

---

## Phase 3 — Outputs and interpretation

Boltz output layout (typical):

```
<OUT_DIR>/
└── boltz_results_<name>/
    └── predictions/
        └── <mol_id>/
            ├── <mol_id>_model_0.cif
            ├── confidence_<mol_id>_model_0.json
            └── affinity_<mol_id>.json
```

Affinity JSON example:

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

Interpretation:
- `affinity_pred_value`: mean predicted log10(IC50) in µM (more negative = stronger)
  - IC50 (µM) = 10^(affinity_pred_value)
  - example: -0.62 → ~0.24 µM
- `affinity_probability_binary`: mean binder probability (0–1)
- numbered variants indicate uncertainty/spread

Confidence context:
- if structure confidence is low (e.g., <0.3), treat affinity as low-trust.

---

## Phase 4 — Integration with Vina/ProLIF (panel, not single-metric)

### 4.1 Why ranks can disagree

Vina approximates scoring from a rigid docking pose; Boltz-2 is a learned predictor on co-folded structure + patterns.
Disagreement is expected.

Observed example (IPF/ALK5 validation):

| Molecule | Vina Score | Boltz-2 IC50 |
|---|---:|---:|
| hinge_1 | -10.01 (best Vina) | 0.47 µM (worst Boltz) |
| hinge_3 | -8.94 (worst Vina) | 0.23 µM (best Boltz) |

### 4.2 Practical decision rule: panel selection

When the shortlist is small (e.g., DFT PASS set):
- pick a **panel** that balances:
  - docking evidence (hinge H-bond + reasonable key-residue coverage)
  - Vina score
  - Boltz-2 affinity value + binder_prob

Do not “optimize” only one metric unless you have real experimental feedback.
