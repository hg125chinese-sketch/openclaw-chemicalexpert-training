# IPF / ALK5 Cycle 6 — Case study (ChemicalExpert)

Status: **Skeleton (waiting for QE DFT results)**

This is the first IPF cycle that exercised the **full upgraded 21-skill pipeline**, adding:
- PoseBusters QC gates (post-gen + post-dock)
- multi-seed docking robustness (3 seeds)
- GNINA CNN rescoring
- PLIF recovery vs co-crystal reference
- panel selection under multi-signal disagreement
- Boltz-2 affinity as a weak orthogonal signal (ALK5-calibrated)

---

## Goal

Run a complete upgraded Cycle 6 loop starting from DiffSBDD generation and ending at a DFT handoff panel.

---

## Step-by-step summary (with key numbers)

### Step 1 — DiffSBDD generation (n≈100 target)

- Generated records: **91**
- RDKit-valid: **85/91 (93.4%)**
- Unique (SMILES): **85**
- Output SDF: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle6_diffsbdd.sdf`

### Step 2 — PoseBusters QC (checkpoint A: post-generation)

- Total: **91**
- PB-valid: **56/91 (61.5%)**
- Improvement vs Cycle 5 PB-valid: **55.1% → 61.5%**

Interpretation:
- Post-gen PoseBusters removed ~**39%** of generated 3D geometries as QC-unsafe.

Outputs:
- Raw PoseBusters report: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_posebusters_postgen.tsv`
- PB survivors: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_postgen_pb_survivors.csv`

### Step 3 — Safety screen (IPF denylist)

- Input: **56** (PB-valid)
- Survivors: **43/56**
- Reject rate: **~23%** (historical low)

Outputs:
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_safety.md`
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_safety_survivors.csv`

### Step 4 — Vina docking (batch)

- Docked: **43** molecules
- Top1 Vina (batch1): **-10.04**

Outputs:
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_docking_batch1_scores.csv`

### Step 5 — ProLIF interactions + PLIF recovery (Top20)

- Applied ProLIF to docking Top20
- Hinge H-bond (True): **7** molecules (Top20 subset)

Outputs:
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_step5_top20_prolif_plif_recovery.csv`

### Step 6 — Multi-seed docking robustness (3 seeds)

- Seeds: 3
- Strict robust gate (pose convergence + hinge 3/3): only **mol_0064** passed

Key observation:
- Multi-seed robustness is extremely selective: **7 hinge=True → 1 robust**.

### Step 7 — GNINA CNN rescoring

- GNINA score-only rescoring completed (ranking signal)

### Step 8 — Multi-signal summary

- Summary table: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_step8_hinge_panel_summary.csv`

Key disagreement example:
- **mol_0064**: Vina **-10.04**, hinge 3/3, but **PLIF recovery = 0.25** (low)

### Step 9 — Panel selection (relaxed hinge robustness)

Hard gate relaxed to:
- `hinge_consistency >= 0.67` (≥2/3 seeds)

Selected panel (K=3):
- `mol_0064` (1.0)
- `mol_0017` (1.0)
- `mol_0007` (0.67)

Report:
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_step9_panel_selection.md`

### Step 10 — Boltz-2 affinity (ALK5)

Boltz-2 remained a **weak signal** on ALK5 (consistent with calibration):
- mol_0064: binder_prob=0.12, IC50~1 µM
- mol_0017: binder_prob=0.06, IC50~1 µM
- mol_0007: binder_prob=0.14, IC50~10 µM

Interpretation:
- `binder_prob` can be used as a weak orthogonal tiebreaker.
- absolute IC50 is not trusted for ALK5.

### Step 11 — QE DFT results (2/2 PASS)

QE completed Cycle 6 DFT QC:
- Results JSON: `/home/node/.openclaw/workspace-quantumexpert/exports/qc_results_cycle6.json`
- Outcome: **2/2 PASS (100%)**

Selected for QE (hinge_consistency=1.0 subset):
- `mol_0064`
- `mol_0017`

Results:
- `mol_0064`: Vina -10.04, gap **2.90 eV**, dipole **3.74 D**, walltime **6.4 h**
- `mol_0017`: Vina -7.25, gap **2.68 eV**, dipole **3.44 D**, walltime **7.2 h**

Multi-objective scoring (decision support):
- `score_final = (-vina) + 0.2*gap - 0.05*dipole`
  - mol_0064: **10.432**
  - mol_0017: **7.617**

Handoff artifacts:
- `exports/qc_handoff_cycle6.json`
- `reports/ipf_cycle6_qc_prescreen.md`
- Integrated table: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle6_qe_qc_multiscore.csv`

---

## New findings (Cycle 6)

1) **PoseBusters as a post-gen gate is high impact**
- It rejects ~39% of generated 3D geometries that look RDKit-valid.

2) **Multi-seed robustness is extremely strict**
- Even among hinge=True, only 1/7 passed strict (3/3) robustness.

3) **Multi-signal disagreement is real and operationally important**
- Vina, GNINA, PLIF recovery, and robustness do not necessarily align.
- A panel strategy is required.

---

## Open questions / next steps

- After QE returns: compare DFT PASS/FAIL against:
  - PLIF recovery
  - multi-seed robustness
  - GNINA CNN scores
  - (weak) Boltz-2 binder_prob

Goal: decide which signals deserve hard gates vs ranking-only status.
