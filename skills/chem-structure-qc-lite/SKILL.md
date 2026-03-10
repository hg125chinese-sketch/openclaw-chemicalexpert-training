---
name: chem-structure-qc-lite
description: >-
  Lightweight geometry/pose QC for AI-generated 3D ligands and docking poses using PoseBusters (bust) + RDKit sanity checks.
  Two checkpoints: post-generation SDF and post-docking poses.
  Outputs standardized qc.tsv + report with pb_valid, warnings, and failure reasons; includes recovery path (RDKit re-embed + minimize).
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧯", "requires": { "bins": ["bust"], "python": ["pandas", "numpy"] } } }
---

# chem-structure-qc-lite — geometry/pose QC that catches “looks fine but is wrong”

Cycle 5 lesson: AI-generated 3D can pass **docking** and even **interaction checks**, yet be geometrically invalid (wrong chirality, bad bond lengths, atom clashes). PoseBusters (2024) is designed to detect exactly these failure modes.

This skill provides a **lightweight** QC gate at two points:
1) **After generation** (DiffSBDD/VAE outputs SDF)
2) **After docking** (Vina poses)

---

## Workspace variables

```bash
OPENCLAW_WORKSPACE="/home/node/.openclaw/workspace-chemicalexpert"
CHEM_PY="/opt/conda/envs/chem/bin/python"
```

---

## Where this fits in the pipeline (integration points)

### Checkpoint A — after safety, before docking (generated 3D)

**Rationale:** remove geometrically broken 3D before it contaminates docking/pose ranking.

Recommended chain:
`generation → SMILES extraction + sanitize → safety denylist → structure-qc-lite (A) → docking`

### Checkpoint B — after docking, before ProLIF (poses)

**Rationale:** catch pose geometries that are unphysical (clashes / broken stereochem), even if docking score is good.

Recommended chain:
`docking → structure-qc-lite (B) → ProLIF interactions → panel selection`

---

## Hard gate

- **Hard gate:** `pb_valid == True` is required to pass the QC stage.
- **Recovery path (optional):** for failures, attempt `RDKit re-embed + minimize`, then re-run PoseBusters.
  - If it becomes valid, mark as `recovered=True` and keep.
  - If still invalid, reject.

Why hard gate exists:
- Cycle 5: DiffSBDD native 3D produced extreme MACE strain (~450–630 kcal/mol). That is a strong indicator the geometry is not QC-grade.

---

## PoseBusters integration

### A) Generated SDF QC

```bash
bust generated.sdf --outfmt long > qc_generated.tsv
```

### B) Docking pose QC

If you have a multi-pose SDF (or one SDF per ligand):

```bash
bust docked_poses.sdf --outfmt long > qc_docked.tsv
```

**Expected output:** a TSV with PoseBusters checks, per molecule/pose.

Notes:
- The CLI is intentionally simple here; treat PoseBusters as the source of truth for “pb_valid”.
- If your pose files include multiple conformers/poses, ensure you keep `mol_id`/pose identifiers stable.

---

## RDKit lightweight checks (fast + interpretable)

These are *supplementary* diagnostics. They do not replace PoseBusters.

### Checks
- **sanitize warnings / failures** (kekulization, valence)
- **fragment count** (unexpected fragmentation)
- **abnormal bond lengths** (very short/long bonds; crude but catches broken geometries)

### Minimal script skeleton

```python
#!/opt/conda/envs/chem/bin/python
import math
import pandas as pd
from rdkit import Chem


def bond_length(p1, p2):
    dx, dy, dz = p1.x - p2.x, p1.y - p2.y, p1.z - p2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def rdkit_qc_sdf(sdf_path: str, short_A=0.9, long_A=2.0):
    """Return per-mol lightweight QC: sanitize_ok, n_frags, n_bad_bonds.

    Thresholds are intentionally crude.
    - short_A/long_A catch obvious geometry breaks and clashes.
    """
    rows = []
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for i, mol in enumerate(suppl):
        mol_id = None
        if mol is None:
            rows.append({"mol_index": i, "mol_id": None, "sanitize_ok": False, "n_frags": None, "n_bad_bonds": None, "rdkit_reason": "SDMolSupplier returned None"})
            continue

        mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"

        # sanitize
        sanitize_ok = True
        reason = ""
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            sanitize_ok = False
            reason = f"sanitize_fail: {e}"[:200]

        # fragments
        n_frags = len(Chem.GetMolFrags(mol))

        # bond lengths (needs conformer)
        n_bad = 0
        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            for b in mol.GetBonds():
                a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                d = bond_length(conf.GetAtomPosition(a1), conf.GetAtomPosition(a2))
                if d < short_A or d > long_A:
                    n_bad += 1
        else:
            reason = (reason + "; no_conformer").strip("; ")

        rows.append({
            "mol_index": i,
            "mol_id": mol_id,
            "sanitize_ok": sanitize_ok,
            "n_frags": n_frags,
            "n_bad_bonds": n_bad,
            "rdkit_reason": reason,
        })

    return pd.DataFrame(rows)
```

---

## Standardized outputs

### Required columns (minimum)
For each mol/pose:
- `mol_id`
- `pb_valid` (boolean)
- `n_warnings` (PoseBusters warning/error count)
- `failure_reasons` (semicolon-separated)

### Recommended extra columns
- `stage` in {`post_gen`, `post_dock`}
- `recovered` in {`True`, `False`}
- RDKit diagnostics: `sanitize_ok`, `n_frags`, `n_bad_bonds`

---

## Recovery path: re-embed + minimize + re-test

When `pb_valid=False`, attempt recovery **from SMILES**:
- ETKDGv3 embed
- MMFF minimize
- write SDF
- re-run PoseBusters

This mirrors the Cycle 5 fix: DiffSBDD 3D was not QC-grade; RDKit re-embed produced QC-feasible geometries.

---

## Reporting checklist

- [ ] Report both checkpoints (A and/or B) and which one was used.
- [ ] Provide `qc.tsv` (PoseBusters) and a compact summary table.
- [ ] For rejects, list top failure reasons.
- [ ] If recovery was used, report recovered count and rerun results.
