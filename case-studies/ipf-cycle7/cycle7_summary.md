# IPF / ALK5 Cycle 7 — Case study (ChemicalExpert)

Status: **Completed through QE DFT for mol_0021**

Cycle 7 is the first IPF cycle to exercise the **full upgraded 24-skill pipeline**, adding three major infrastructure layers on top of the Cycle 6 workflow:
- **Entity resolution** (`chem-entity-resolver`, skill 24)
- **Evidence schema** (`chem-evidence-schema`, skill 23)
- **Three-layer safety** (`chem-reactivity-safety`, upgraded skill 14)

This cycle also reused the upgraded structure-based stack from Cycle 6:
- PoseBusters QC gates
- multi-seed docking robustness
- GNINA CNN rescoring
- PLIF recovery vs co-crystal reference
- panel selection under multi-signal disagreement
- Boltz-2 as a weak orthogonal signal on ALK5

---

## Goal

Run the first **fully infrastructure-upgraded** ALK5 / IPF cycle:
- canonicalize entities first
- apply stronger post-generation and safety gates
- preserve evidence provenance in a normalized schema
- rank finalists under explicit multi-signal disagreement

---

## Step-by-step summary (with key numbers)

### Step 0 — Entity resolution (new)

Resolved canonical IDs before the rest of the pipeline:
- target: **TGFBR1 / ENSG00000106799 / P36897**
- disease: **idiopathic pulmonary fibrosis / EFO_0000768**

Interpretation:
- This is the first cycle where target / disease IDs were standardized up front instead of implicitly carried through scripts.

### Step 1 — DiffSBDD generation

- Generated records: **98**
- RDKit-valid: **93/98 (94.9%)**
- Unique: **93**

Interpretation:
- This is the best raw generation validity seen so far in the IPF sequence.

### Step 2 — PoseBusters QC (checkpoint A: post-generation)

- Input (RDKit-valid): **93**
- PB-valid: **61/93 (65.6%)**
- Historical comparison: **highest post-gen PB-valid rate so far**

Interpretation:
- Post-gen geometry quality improved again relative to earlier cycles.
- PoseBusters continues to remove a large fraction of geometrically unsafe 3D structures that would otherwise pass RDKit validity.

### Step 3 — Three-layer safety (new)

Three-layer safety stack:
1. SMARTS structural alerts
2. ADMET-AI layer
3. FAERS / real-world evidence layer

Results:
- Input: **61** PB-valid molecules
- Survivors: **51**
- Reject rate: **~16%**

Interpretation:
- This is the **lowest reject rate so far** in the IPF program.
- The generated pool looks cleaner both geometrically and structurally.

### Step 4 — Vina docking (batch)

- Docked: **51** molecules
- Top1 Vina: **-10.27**

Interpretation:
- The dockable pool remains strong even after the stricter upstream filters.

### Step 5 — ProLIF interactions + PLIF recovery (Top20)

- Applied ProLIF + PLIF recovery to docking Top20
- Hinge H-bond positive: **5** molecules

Interpretation:
- Only a minority of strong docking hits recover hinge-like interaction evidence under the reference-guided interaction layer.

### Step 6 — Multi-seed docking robustness (3 seeds)

- Hinge-positive candidates tested with 3 seeds
- Strict robust pass: **only `mol_0021`**
- `mol_0021`: hinge **3/3**, robust=True

Interpretation:
- Multi-seed robustness remains highly selective.
- As in Cycle 6, most docking hits do not survive pose-stability scrutiny.

### Step 7 — GNINA CNN rescoring

- GNINA rescoring completed for the hinge panel
- This was used as an additional ranking signal, not a hard gate

### Step 8 — Multi-signal summary

Signals merged:
- Vina
- hinge consistency / robustness
- PLIF recovery
- GNINA CNN affinity

This preserved disagreement rather than forcing a single score too early.

### Step 9 — Boltz-2 affinity

Panel Boltz-2 result highlight:
- **`mol_0021`: binder_prob = 0.698**

Interpretation:
- This is the **highest Boltz-2 binder probability seen across all ALK5 cycles so far**.
- Compared against the ALK5 calibration set (active mean ≈ **0.268**), **0.698 is substantially above the known-active average**.
- This does **not** make Boltz-2 decisive, but it is a strong positive orthogonal signal relative to prior ALK5 cycles.

### Step 10 — Panel selection with evidence schema (new)

Panel ranking used:
- hard gate: hinge consistency
- ranking signals: Vina + GNINA + PLIF recovery + Boltz-2
- evidence stored in standardized evidence objects (`skill 23`)

Final recommendation:
- **`mol_0021`**

Why:
- only strict robust pass
- hinge 3/3
- acceptable PLIF recovery
- supportive GNINA
- strongest Boltz-2 signal in the entire ALK5 program

### Step 11 — QE DFT handoff and result

Selected for QE:
- **`mol_0021`** only

Handoff preparation used:
- RDKit ETKDGv3 + MMFF re-embed
- MACE prescreen
- geometry source = `rdkit_embed`

Artifacts:
- `exports/qc_handoff_cycle7.json`
- `reports/ipf_cycle7_qc_prescreen.md`
- `/home/node/.openclaw/workspace-quantumexpert/exports/qc_results_cycle7.json`
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle7_qe_qc_multiscore.csv`
- `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle7_qe_qc_summary.md`

QE result for `mol_0021`:
- **QC flag: PASS**
- method: **B3LYP-D3(BJ)/def2-SVP**
- `E_total = -1008.623 Ha`
- `gap = 2.95 eV`
- `dipole = 7.21 D`
- `n_imag = 0`
- walltime: **1.9 h**

Cycle 7 multi-objective score:
- `score_final = (-vina) + 0.2*gap - 0.05*dipole`
- for `mol_0021`: **9.136**

Interpretation:
- Cycle 7 is the **first ALK5 cycle** to combine:
  - canonical entity resolution,
  - three-layer safety,
  - evidence-schema panel adjudication,
  - strict multi-seed robustness,
  - a very strong Boltz-2 signal,
  - and a **DFT PASS** outcome
  in the same lead candidate.
- This is also the first cycle using **true D3 dispersion** (`B3LYP-D3(BJ)`) on the QE side.

---

## Key new features in Cycle 7

1. **First use of entity resolver (skill 24)**
- Canonical target/disease IDs were established before downstream queries.

2. **First use of three-layer safety**
- SMARTS + ADMET-AI + FAERS / real-world context.

3. **First use of evidence schema (skill 23)**
- Panel signals were written in a standardized evidence format with provenance and conflicts.

4. **Best Boltz-2 signal seen on ALK5 so far**
- `mol_0021` binder_prob = **0.698**, well above the calibration active mean (~0.268).

---

## Main takeaways

1. **Upstream infrastructure upgrades are paying off**
- post-gen PB-valid is the highest yet
- safety reject rate is the lowest yet

2. **The selection funnel is getting narrower for the right reasons**
- docking alone still over-selects
- hinge + PLIF + multi-seed robustness continue to do the real filtering

3. **Cycle 7 produced the strongest ALK5 orthogonal ML signal seen so far**
- `mol_0021` is the first candidate where structural robustness and Boltz-2 are both supportive at the same time

4. **Cycle 7 closed the loop successfully with DFT PASS**
- `mol_0021` passed QE at **B3LYP-D3(BJ)/def2-SVP**
- `gap = 2.95 eV`, `dipole = 7.21 D`, `n_imag = 0`
- final triage score = **9.136**

---

## Next steps

- Compare Cycle 7 directly against Cycle 6 in a final cross-cycle ranking table
- Decide whether `mol_0021` should advance as:
  - the primary ALK5 / IPF Cycle 7 lead, or
  - a complementary lead beside Cycle 6 `mol_0064`
- Continue improving the learned relationship between docking / robustness / Boltz signals and QE outcomes
