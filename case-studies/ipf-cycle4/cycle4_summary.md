# IPF / ALK5 (TGFBR1) Cycle 4 — Summary (CE↔QE end-to-end loop)

**Objective:** run a complete CE↔QE collaboration loop for ALK5 and accumulate more data to calibrate the **MACE prescreen** before DFT QC.

---

## 1) Pipeline snapshot (Cycle 4)

**Generation (Strategy E: logit-bias decoding)**
- Base model: SELFIES GRU VAE (healthy checkpoint)
- Targeted generation: **1,000** accepted unique molecules (hinge-motif present)

**Safety (IPF chronic denylist, hard reject)**
- Input: 1,000
- **Survivors: 466**

**Docking (ALK5 1VJY, Vina)**
- Docked batch: first **50** safety survivors
- Scored: **49** (1 docking error)

**Interactions (chem-docking-interactions, ProLIF)**
- From Vina Top10 → interaction rerank → Cycle 4 Top5
- **Top5 hinge H-bond: 3/5**

---

## 2) QC handoff to QuantumExpert (pre-screen + DFT QC)

**QC prescreen inputs:** Cycle 4 Top5, filtered to hinge H-bond = True (max 3)

Prescreen gates:
1. RDKit 3D embedding + MMFF optimization (must converge)
2. Dominant protomer (pH 7.4 proxy) → formal charge
3. MACE-OFF prescreen (must converge; strain threshold = 50 kcal/mol)

**Outcome:**
- Selected (hinge H-bond = True): 3
- Rejected before QE: 1 (RDKit embedding failed)
- Sent to QE: **2**
- QE DFT QC: **1 PASS, 1 OPT_FAIL**

---

## 3) MACE calibration update

New data point confirms the earlier observation:

- The current MACE prescreen **strain** proxy does **not** correlate with DFT OPT_FAIL.
- A recurrent OPT_FAIL structure was observed again (same SMILES as a prior OPT_FAIL), while still passing the MACE prescreen.

**Interpretation:** strain (as currently defined by single-geometry relaxation energy under MACE-OFF) is not sufficient as a DFT-feasibility predictor in this workflow.

---

## 4) Retrosynthesis status

Retrosynthesis + reaction-conditions annotations have been completed for the final **DFT PASS** set (see project conclusion report).
