# ALK5 hinge-binder motif coverage (heuristic)

Train set (ChEMBL260) N_rdkit_ok=5067

Candidate pool N_rdkit_ok=2999

| motif | train_count | train_frac | pool_count | pool_frac |
|---|---:|---:|---:|---:|
| pyrimidine_ring | 29 | 0.006 | 1 | 0.000 |
| imidazole_ring | 0 | 0.000 | 0 | 0.000 |
| imidazole_ring_alt | 622 | 0.123 | 3 | 0.001 |
| pyrazole_ring | 1184 | 0.234 | 18 | 0.006 |
| quinazoline_core | 1 | 0.000 | 0 | 0.000 |
| 2-aminopyrimidine | 32 | 0.006 | 1 | 0.000 |
| 2-aminopyridine | 933 | 0.184 | 27 | 0.009 |
| pyrazolo_pyrimidine | 0 | 0.000 | 0 | 0.000 |

Interpretation: if pool fractions are near-zero for key motifs while train has substantial coverage,
it suggests the generator drifted away from kinase-like chemistry.
