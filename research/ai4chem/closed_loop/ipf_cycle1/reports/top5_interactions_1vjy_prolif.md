# IPF Cycle1 Top5 — docking interaction analysis (ProLIF)

Target PDB: 1VJY | redocking RMSD (best): 0.045 Å

Reference ligand: cocrystal_lig_ph7.4.sdf → interaction profile saved to `interaction_reference_1vjy_prolif.json`.

## Top5: hinge H-bond + key residue coverage + interaction rerank

| rank_interaction | rank_vina | smiles | vina | hinge_hbond | key_cov | n_contacts | interaction_score | flags |
|---:|---:|---|---:|---|---:|---:|---:|---|
| 1 | 4 | `Cc1ccc(F)cc1CC=Cc1ccccc1Cl` | -9.591 | True | 0.111 | 2 | 8.83 |  |
| 2 | 3 | `CCNC(=O)c1cccc(Cc2ccccn2)c1` | -8.701 | False | 0.222 | 5 | 3.17 | CRITICAL: No hinge H-bond detected |
| 3 | 5 | `CC(=O)Nc1cccnc1NC=CNc1cc(F)ccc1C` | -8.653 | False | 0.222 | 5 | 3.17 | CRITICAL: No hinge H-bond detected |
| 4 | 2 | `Cc1cccc(Nc2ccccc2Cl)c1` | -8.117 | False | 0.111 | 4 | 2.33 | CRITICAL: No hinge H-bond detected |
| 5 | 1 | `CSC(=O)c1cccc(COc2cc(F)ccc2C)c1` | -8.574 | False | 0.000 | 3 | 1.50 | CRITICAL: No hinge H-bond detected |

Ranking changed vs pure Vina? **True**
