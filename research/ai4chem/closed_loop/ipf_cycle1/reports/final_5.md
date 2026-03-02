# DRD2 closed-loop: 5 candidates

Selection constraints: novelty vs ChEMBL-train; safety/drug-likeness first (Lipinski+Veber pass; QED>0.3); single-step retrosynthesis template with forward validation; SA score available; SA<=6.0; |RF-GCN|<=1.0; filtered motifs: C=C=N, [SH].

| # | SMILES | pred pChEMBL (RF) | pred pChEMBL (GCN) | QED | MW | LogP | TPSA | HBD | HBA | Lipinski | Veber | Traffic (Oral/QED/BBB/LogP) | SA | Retro step |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---|
| 1 | CC(N)[C@@H]1N=C(N)C=CC=CC=C(C=O)C(c2c(F)cccc2Cl)=Cc2cc(F)ccc2N1 | 7.80 | 7.75 | 0.561 | 468.9 | 4.85 | 93.5 | 3 | 5 | pass | pass | green/yellow/yellow/green | 4.75 | reductive_amination |
| 2 | CCOC=CC=NC1C=C2C=NCN=C2CCc2ccc(F)cc2C(=O)N1 | 7.62 | 7.79 | 0.655 | 368.4 | 2.86 | 75.4 | 1 | 5 | pass | pass | green/green/green/green | 4.82 | reductive_amination |
| 3 | O=C(NC1N=CC=CC2OCC=C2N1)c1ccccn1 | 7.16 | 7.05 | 0.825 | 270.3 | 0.61 | 75.6 | 2 | 5 | pass | pass | green/green/green/green | 4.26 | amide_formation |
| 4 | COc1ccccc1C=CC=NC=C1CCCC(=O)NC=C(CNC(=O)c2ccc(F)cc2F)N(C)CC1 | 7.25 | 7.21 | 0.510 | 522.6 | 4.83 | 83.0 | 2 | 5 | pass | pass | green/yellow/green/green | 3.84 | amide_formation |
| 5 | C=C1C(C=CCNCC2(CCBr)C=CC=CO2)=CCc2ccc(F)cc21 | 7.53 | 7.00 | 0.494 | 416.3 | 5.09 | 21.3 | 1 | 2 | pass | pass | green/yellow/green/yellow | 4.46 | reductive_amination |
