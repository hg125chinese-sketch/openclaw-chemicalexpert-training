# IPF / ALK5 — Cycle 5 (Pocket-conditioned Diffusion) Summary

## Goal

Cycle 5 replaced the string-VAE generator with **DiffSBDD**, a **pocket-conditioned 3D diffusion** model, to test whether structure-based generation can improve downstream quality (safety, binding interactions, and QC pass rate).

## What changed vs Cycles 2–4

- **Generator**: DiffSBDD (first time using pocket-conditioned 3D diffusion in this project)
- **Key workflow change**: treat DiffSBDD primarily as a *molecular design engine*, not as a trusted source of QC-ready 3D geometry.

## Results (end-to-end)

### Step 1 — Generation
- Generated: **98** molecules
- Valid: **89 / 98 = 90.8%**
- Unique: **89**

### Step 2 — Safety (chronic project denylist)
- Hard reject rate: **28%**
  - Reference: VAE-based Cycle 4 hard reject was **53%**

### Step 3 — Docking (batch1)
- Docked: **50** safety survivors (first batch)
- Interaction-aware Top5:
  - **3 / 5** have a **hinge H-bond**

### Step 4 — QC prescreen (MACE gate)
- Attempt A (failed): reuse DiffSBDD-provided 3D coordinates
  - Observed MACE strain (|ΔE| proxy): **~450–630 kcal/mol** for hinge-positive candidates
  - Conclusion: diffusion geometries were not suitable as direct QC inputs.
- Attempt B (successful): **RDKit ETKDGv3 + MMFF re-embed**, then MACE prescreen
  - Outcome: hinge-positive set passed the prescreen and was sent to DFT.

### Step 5 — DFT QC (QE)
- DFT outcome: **3 / 3 PASS = 100%**
  - Reference: Cycle 3/4 QC pass rate was **~50%**

## Best candidate

- **cycle5_top5_hinge_1**
  - Vina: **-10.01**
  - Multi-objective: **score_final = 10.404**

## Key finding

DiffSBDD added clear value at the **molecule design** level (validity, lower safety reject rate, and strong docking/QC outcomes). However, its raw 3D coordinates were not QC-ready in this setup: **QC prescreen required RDKit re-embedding** to produce geometries that are compatible with the MACE strain gate and downstream DFT workflows.

In short:
- **Use DiffSBDD for proposing molecules**.
- **Do not assume diffusion 3D is QC-ready**; re-embed (RDKit) before QC.
