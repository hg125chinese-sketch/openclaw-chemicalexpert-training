# NNF_05 successor campaign conclusion

## Bottom line

**NNF05_S05 is the current campaign winner.**

It preserves the core properties that made the NNF_05 line attractive:
- N-N free
- PoseBusters PASS
- safety PASS
- hinge 3/3
- DFT PASS

and it improves on both the parent continuity lead NNF_05 and the solvent-edge hybrid successor NNF05_S10 in the most decision-relevant structure/physics signals.

## Key comparison: S05 vs S10

| molecule | Vina | Boltz binder_prob | hinge 3/3 | gap (eV) | dipole (D) | DFT |
|---|---:|---:|---:|---:|---:|---:|
| NNF05_S05 | -9.642 | 0.505 | True | 4.97 | 1.54 | True |
| NNF05_S10 | -9.075 | 0.591 | True | 4.76 | 2.67 | True |

Interpretation:
- **S05 wins on docking strength** (`-9.64` vs `-9.07`)
- **S05 wins on QE-style stability proxies** (higher gap, lower dipole)
- **S10 remains Boltz-favored** (`0.591` vs `0.505`), so it is still the best continuity / orthogonal-support backup
- Because both molecules are hinge 3/3 and DFT PASS, the stronger Vina + better gap/dipole package makes **S05 the more advanced successor lead**

## Campaign-level conclusion

This successor campaign validates that the NNF_05 scaffold is not merely preserved but **optimizable**.
The benzofuran-side **diF variant S05** is the first local modification that upgrades the de-risked lead without breaking the validated ALK5 binding-mode logic.

A second important campaign-level lesson is that the **DiffSBDD era is now operating at 14/14 = 100% DFT PASS** for the validated set. That strongly supports the current geometry rebuild + QC handoff policy as a stable production workflow rather than a one-off success.

## Recommended decision

- **Primary successor lead:** `NNF05_S05`
- **Orthogonal / continuity backup:** `NNF05_S10`
- **Project lesson:** local fluorine patterning on the benzofuran side is now a validated optimization direction for the NNF_05 scaffold.
