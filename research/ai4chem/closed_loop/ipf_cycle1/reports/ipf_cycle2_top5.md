# IPF Cycle 2 — Top5 (from batch1 Top20) after interaction rerank (ProLIF-lite)

Docking source: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle2_docking_batch1_scores.csv` (n=20) → select Vina Top10 → interaction scoring → Top5.

Target PDB: 1VJY | co-crystal redock RMSD(best): 0.045 Å

| rank | smiles | vina | interaction_score | hinge_hbond | hinge_details | key_cov | flags |
|---:|---|---:|---:|---|---|---:|---|
| 1 | `CC1COCCN1c1ccc(=O)n(CN=CC=CNCC=C(C=CNC(=O)NC2=CC=C[C@H](C)C2)C(C)(C)C)c1N` | -7.921 | 10.167 | True | SER280.A HBAcceptor | 0.222 |  |
| 2 | `CC[C@H](C)Nc1cc(OC=CC(C(N)=O)[C@@H]2CCCOCCN3C(C2)C2=CN(C=CC=N2)C3CO)ccn1` | -8.166 | 10.000 | True | HIS283.A HBAcceptor | 0.333 |  |
| 3 | `[CH]CCOC=NC=CCNC(F)C(=O)CNc1ccccn1` | -8.931 | 8.667 | True | HIS283.A HBAcceptor | 0.222 |  |
| 4 | `CC1(C)C=CC=CN2N=C1Cc1ccc3ncnn3c1N1C=CC21` | -8.304 | 6.333 | True | HIS283.A HBAcceptor | 0.111 | low_key_residue_coverage |
| 5 | `CNC(=O)CN[C@@H](C)NC=NC=CC=CC(=CC=CC=NCNc1cccc(C)n1)c1ccc(F)cc1` | -8.749 | 4.667 | False |  | 0.222 | no_hinge_hbond |
