# Analog exploration summary — mol_0021 neighborhood (ChemicalExpert)

Status: **Completed through QE DFT**

This case study archives the focused analog campaign launched from the Cycle 7 ALK5 / IPF lead `mol_0021`.

Goal:
- test whether `mol_0021` was a one-off hit or the anchor of a reproducible local chemotype neighborhood
- preserve the Cycle 7 structural logic (hinge, multi-seed, GNINA, Boltz, QE)
- compare analogs against the parent in a conflict-aware way

---

## Starting point

Parent lead:
- `mol_0021`
- key Cycle 7 properties:
  - strict robust
  - hinge-stable
  - Boltz-2 binder_prob **0.698**
  - QE PASS

Main question:
- can close analogs retain or improve this cross-signal coherence?

---

## Analog mini-panel

A 10-member analog panel was designed around the `mol_0021` scaffold across five SAR directions:
1. F-scan
2. hinge-facing aza-aryl tuning
3. bridge analogs / stereochemistry
4. polarity trimming
5. solvent-edge growth

All 10 analogs:
- passed RDKit parsing
- stayed within simple Lipinski-like bounds
- passed PoseBusters QC in the executed analog validation stack

---

## Validation outcomes

### Structural results
- **7/10** analogs retained hinge H-bond (**70%**)
- strongest structure-based analogs:
  - `A5_01`: Vina **-10.12**, hinge=True
  - `A3_02`: Vina **-10.00**, hinge=True
  - `A3_01`: Vina **-9.76**, hinge=True

### Multi-seed robustness
Top 3 analogs all passed strict 3-seed robustness:
- `A5_01`: robust=True, hinge 3/3
- `A3_02`: robust=True, hinge 3/3
- `A3_01`: robust=True, hinge 3/3

### GNINA rescoring
- `A5_01`: strongest CNNscore
- `A3_01`: strongest CNNaffinity
- `A3_02`: weaker GNINA than the other two, but still robust and hinge-stable

### Boltz-2 affinity
- `A5_01`: binder_prob **0.570**
- `A3_02`: binder_prob **0.704**
- `A3_01`: binder_prob **0.441**

Interpretation:
- `A3_02` exceeded the parent on Boltz (`0.704 > 0.698`), becoming the strongest Vina–Boltz consensus analog

### Safety context
- all analogs remained **safety-constrained** because this scaffold family carries `n_n_single_bond` / chronic N–N review flags
- the series therefore remains promising but safety-constrained, not cleanly de-risked

---

## QE DFT outcome

Top 3 analogs were sent for QE DFT and all passed:
- `A5_01`: PASS, gap **2.93 eV**, dipole **8.15 D**, walltime **2.8 h**
- `A3_02`: PASS, gap **3.56 eV**, dipole **6.55 D**, walltime **1.8 h**
- `A3_01`: PASS, gap **4.79 eV**, dipole **1.62 D**, walltime **2.1 h**

Key implication:
- the `mol_0021` analog neighborhood is not a brittle docking artifact
- it survives the full CE↔QE validation stack

---

## Final analog reading

### A5_01
- best docking / CNNscore analog
- likely a successful peripheral-fit optimization of the parent pose core

### A3_02
- best Vina–Boltz consensus analog
- also QE-balanced
- strongest all-around analog lead after full integration

### A3_01
- best QE-style gap/dipole analog
- strongest GNINA-affinity analog
- valuable alternate hypothesis rather than the primary consensus winner

---

## Main lesson

`mol_0021` was **not** a one-off hit.

The analog campaign showed that:
- the neighborhood is real
- the top analogs are hinge-stable
- the top analogs are robust across seeds
- the top analogs survive QE DFT

The limiting factor is no longer whether this chemotype can work at all.
The main remaining limitation is that the whole series remains under a scaffold-level **N–N safety caution**, so further progression should stay explicit about that tradeoff.
