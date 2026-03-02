# IPF / ALK5 (TGFBR1) Cycle 1 — Summary

**Goal:** establish a complete DMTA baseline (data → QSAR → generation → rules-first filtering → docking → reranking) for ALK5/TGFBR1. For a chronic indication (IPF), we prioritize **safety / developability > activity**.

---

## (1) Key metrics

### Data (ChEMBL)
- Target: **CHEMBL260 (TGFBR1/ALK5)**
- Raw activity records iterated: **16,324**
- Unique molecules (median pChEMBL aggregated): **5,067**
- File: `data/chembl_alk5_activity.csv`

### QSAR / GNN (scaffold split)
- Split: N=5,067 (train=4,053, valid=507, test=507)
- RF (Morgan, 500 trees): RMSE **0.586**, R² **0.567**
- GCN (3-layer, hidden=128): RMSE **0.821**, R² **0.150** (worse than RF)
- File: `reports/model_metrics.md`

### Generation (SELFIES VAE)
- Trained 8 epochs; **KL collapse** observed (AU ≈ 0/32)
- Generated valid unique SMILES: **4,936**
- File: `generated/generated_raw.csv`

### Candidate pool (rules-first)
- Rules-first ADMET (Lipinski + Veber + QED > 0.3) + novelty → pool: **N=2,999**
- With single-step template routes + forward validation: **1,689**
- Files: `reports/candidate_pool.csv`, `reports/routes_5.md`

### Docking QC (chem-docking)
- PDB: **1VJY**; co-crystal ligand resname: **460**
- Redocking: **passed** (best heavy-atom RMSD ≈ **0.05 Å**)
- Co-crystal redock mode1 score ≈ **-10.23 kcal/mol**
- Docking scale-up: Top300 (successful scores for **297** molecules)
  - `reports/top300_dock_scores.csv`
  - reranked scored table: `reports/top300_scored.csv`

---

## (2) Issues discovered & fixes implemented

### Issue A: Generated distribution drifted away from kinase-like chemical space
This is more fatal than downstream filtering.

- Evidence: hinge-binder / kinase-like motif coverage in `candidate_pool` is near-zero, while the training set has substantial coverage.
- Outputs: `reports/scaffold_coverage.md` / `reports/scaffold_coverage.json`

### Issue B: Chronic safety context (IPF) discourages N–N / hydrazone / hydrazine motifs
- Action: promote N–N / hydrazone-like motifs to a **hard denylist** (hard reject, not downweight).

### Issue C: Docking scale-up can fail due to environment/toolchain brittleness
- Action: run Top300 docking in chunks and merge.
- Engineering: docking script includes a fallback path (when meeko CLI is unavailable) using OpenBabel to write ligand PDBQT.

---

## (3) Cycle 2 recommendation (generator-first)

**Key conclusion:** Cycle 2 should primarily fix the generator to pull samples back into kinase/ALK5-relevant chemical space. Otherwise, more docking/filtering just means selecting the “best” molecules from the wrong distribution.

Recommended priority order:
1. **Prior-informed generation:** condition on hinge-binder fragments / scaffold constraints to guarantee candidates contain a plausible hinge binder.
2. **Anti-collapse training:** KL annealing / cyclical annealing / free-bits / word dropout to restore latent usage.
3. **Representation / model family change:** consider graph-based models, JT-VAE, fragment-based models, Transformer + constraints.
4. **Training data alignment:** audit frequent scaffolds/motifs in ALK5 actives and feed them back into the generator.

Minimum verifiable milestones for Cycle 2:
- Hinge-binder motif coverage ratio (candidate vs train) is no longer near-zero.
- Top300 docking best scores move closer to the co-crystal score (smaller Δ), rather than being consistently worse by >2–3.
