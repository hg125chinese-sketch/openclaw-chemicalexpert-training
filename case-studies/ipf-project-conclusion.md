# IPF / ALK5 (TGFBR1) — Project Conclusion (CE↔QE)

This report summarizes the full IPF/ALK5 project across **7 cycles**, focusing on what was learned about conditioned generation, docking/interaction KPIs, and the CE↔QE collaboration loop for DFT QC.

---

## 1) Seven-cycle comparison (high-level)

| Cycle | Main focus | Generator strategy | Key KPI outcome (Top5 hinge H-bond) | QE DFT QC status | Safety hard reject |
|---:|---|---|---|---|---|
| 1 | Baseline DMTA | SELFIES GRU VAE (unconditioned / weak control) | **1/5** hinge H-bond | Not used | (not tracked) |
| 2 | Fix generator drift | Strategy D (rejection sampling) | **4/5** hinge H-bond | Not used | (tracked; cycle-specific) |
| 3 | Improve efficiency vs rejection | Strategy E (logit-bias decoding) succeeded; cross-attn conditioning failed | **4/5** hinge H-bond | 4 sent → 2 PASS / 2 OPT_FAIL | (tracked; cycle-specific) |
| 4 | Full CE↔QE loop + prescreen calibration | Strategy E scaled; docking batch1; QC prescreen + QE | **3/5** hinge H-bond (Top5 from docked 50) | 2 sent → 1 PASS / 1 OPT_FAIL | **53%** (VAE-based) |
| 5 | Pocket-conditioned generation + QC stabilization | **DiffSBDD** (pocket-conditioned 3D diffusion) | **3/5** hinge H-bond (Top5 from docked 50) | 3 sent → 3 PASS / 0 OPT_FAIL | **28%** |
| 6 | Full upgraded structure pipeline (PoseBusters + multi-seed + GNINA + PLIF recovery + panel selection) | DiffSBDD + upgraded gates | (panel-driven; hinge robustness + recovery used) | 2 sent → 2 PASS / 0 OPT_FAIL | **23%** (post-PB set) |
| 7 | Full infrastructure-upgraded pipeline (entity resolver + three-layer safety + evidence schema + upgraded structure stack) | DiffSBDD + full 24-skill pipeline | **5 hinge+ Top20; 1 strict robust lead** | 1 sent → 1 PASS / 0 OPT_FAIL | **16%** (post-PB set) |

Notes:
- “Top5 hinge H-bond” is evaluated on docked poses via interaction analysis (ProLIF), not inferred from motifs alone.
- Cycle 7 is the first cycle to standardize entity IDs up front and preserve panel evidence in a normalized schema.

---

## 2) Final DFT PASS set (multi-objective ranking)

DFT PASS molecules discovered in the CE↔QE loop:

| Final rank* | Cycle | ID | Vina | gap (eV) | score_final | Comment |
|---:|---:|---|---:|---:|---:|---|
| 1 | 6 | mol_0064 | **-10.040** | 2.90 | **10.432** | Best PASS so far by the current score (very strong docking; multi-seed hinge robust) |
| 2 | 5 | cycle5_top5_hinge_1 | -10.010 | 2.71 | 10.404 | Prior best PASS (Cycle 5) |
| 3 | 4 | cycle4_top5_hinge_2 | -9.102 | 2.38 | 9.424 | Prior best PASS before Cycle 5 |
| 4 | 5 | cycle5_top5_hinge_2 | -9.134 | 2.38 | 9.315 | PASS; higher dipole penalty |
| 5 | 5 | cycle5_top5_hinge_3 | -8.935 | 2.42 | 9.186 | PASS |
| 6 | 7 | mol_0021 | -8.906 | 2.95 | 9.136 | PASS; strongest Boltz-2 signal so far, first full infrastructure-upgraded winner |
| 7 | 3 | cycle3_top5_hinge_4 | -8.295 | 2.77 | 8.658 | PASS |
| 8 | 3 | cycle3_top5_hinge_2 | -8.003 | 2.09 | 8.354 | PASS |
| 9 | 6 | mol_0017 | -7.253 | 2.68 | 7.617 | PASS (weaker docking but QC PASS; included for completeness) |

\*Ranking basis: a simple, auditable multi-objective score used for handoff triage:
- primary: docking (more negative Vina is better)
- secondary: larger HOMO–LUMO gap (proxy stability)
- tertiary: smaller dipole (proxy polarity)

This is a decision-support ranking, not a claim of true binding affinity.

---

## 3) Key findings

### A) Training-time conditioning can be ignored by string VAEs
- Multiple Cycle 2 “conditioning” variants showed that a SELFIES GRU VAE decoder can produce valid molecules while structurally **ignoring conditioning**.
- This was only visible when evaluated against a **hard domain KPI** (hinge H-bond / hinge-motif coverage), not by generic VAE diagnostics alone.

### B) Inference-time control works: Strategy E (logit-bias decoding)
- Logit-bias decoding reliably steered generation toward hinge-motif chemistry.
- It improved efficiency over rejection sampling while keeping high hinge-H-bond rates in Top5 for Cycle 3.

### C) MACE prescreen calibration: current strain proxy does not predict DFT OPT_FAIL
- With the current definition (single-geometry relaxation energy under MACE-OFF), both PASS and OPT_FAIL can exhibit similar strain values.
- A recurrent OPT_FAIL SMILES was observed again in Cycle 4, strengthening the conclusion that strain alone is insufficient.

### D) DiffSBDD geometry is not QC-ready (in this setup)
- In Cycle 5, directly reusing **DiffSBDD-provided 3D coordinates** produced extremely high MACE relaxation/strain values (**~450–630 kcal/mol**), indicating geometries that are not suitable as direct QC inputs.
- Switching to **RDKit ETKDGv3 + MMFF re-embedding** before MACE prescreen restored a stable QC workflow and achieved **3/3 PASS (0% failure)**.
- Cycle 7 confirmed that this stabilized handoff recipe remains valid in the upgraded pipeline: `mol_0021` passed QE after RDKit re-embed + MACE prescreen.

### E) Infrastructure upgrades improved candidate quality before docking and after docking
- Cycle 7 added **entity resolution**, **three-layer safety**, and **evidence-schema panel adjudication**.
- The cycle achieved the best upstream gate profile so far:
  - **65.6%** post-gen PoseBusters pass rate
  - **16%** safety reject rate
- It also produced the strongest ALK5 Boltz-2 signal seen so far (`mol_0021`, binder_prob **0.698**), though Boltz remains a weak orthogonal signal rather than a primary gate.

---

## 4) CE↔QE collaboration statistics

Across the CE↔QE QC loop:
- Sent to QE (DFT QC): **12 molecules**
- DFT PASS: **9**
- DFT OPT_FAIL: **3**
- Fail rate: **25%**

Breakdown:
- Cycle 3: 4 sent → 2 PASS / 2 OPT_FAIL
- Cycle 4: 2 sent → 1 PASS / 1 OPT_FAIL
- Cycle 5: 3 sent → 3 PASS / 0 OPT_FAIL
- Cycle 6: 2 sent → 2 PASS / 0 OPT_FAIL
- Cycle 7: 1 sent → 1 PASS / 0 OPT_FAIL

---

## 5) Retrosynthesis (close-out)

For all DFT PASS molecules, a retrosynthesis + reaction-conditions annotation has been prepared (manual, heuristic, template-guided).

---

## 6) Limitations

- Docking scores are not binding affinities; pose quality and protein state remain key uncertainties.
- Interaction KPIs depend on docking pose correctness and the chosen fingerprint/interaction definitions.
- The retrosynthesis/conditions section is heuristic without automated route search and purchasable building-block validation.
- MACE prescreen needs richer signals (optimizer diagnostics, multi-conformer screening, or additional geometry sanity checks) to predict DFT feasibility.

---

## 7) Future work

1) Expand the docking batch beyond 50 (Cycle 4) and rerun interaction rerank to reduce sampling noise.
2) Upgrade QC prescreen:
   - multi-seed RDKit embed retries,
   - multi-conformer MACE screening,
   - track MACE optimizer diagnostics (steps, forces),
   - consider a repeated-OPT_FAIL blacklist only if failures persist.
3) Add a more formal multi-objective model (learned from CE↔QE outcomes) once more labeled data accumulates.
4) Convert retrosynthesis to concrete purchasable routes via an external planner + forward validation.
