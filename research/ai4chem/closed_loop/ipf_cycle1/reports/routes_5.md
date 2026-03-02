# Routes for 5 candidates (single-step, template-based)

| # | Target | Template | Precursors | Forward-validated |
|---:|---|---|---|---|
| 1 | CC(N)[C@@H]1N=C(N)C=CC=CC=C(C=O)C(c2c(F)cccc2Cl)=Cc2cc(F)ccc2N1 | reductive_amination | CC(=O)[C@@H]1N=C(N)C=CC=CC=C(C=O)C(c2c(F)cccc2Cl)=Cc2cc(F)ccc2N1 + N | yes |
| 2 | CCOC=CC=NC1C=C2C=NCN=C2CCc2ccc(F)cc2C(=O)N1 | reductive_amination | O=C1C=C2C=NCN=C2CCc2ccc(F)cc2C(=O)N1 + CCOC=CC=N | yes |
| 3 | O=C(NC1N=CC=CC2OCC=C2N1)c1ccccn1 | amide_formation | O=C(O)c1ccccn1 + NC1N=CC=CC2OCC=C2N1 | yes |
| 4 | COc1ccccc1C=CC=NC=C1CCCC(=O)NC=C(CNC(=O)c2ccc(F)cc2F)N(C)CC1 | amide_formation | O=C(O)c1ccc(F)cc1F + COc1ccccc1C=CC=NC=C1CCCC(=O)NC=C(CN)N(C)CC1 | yes |
| 5 | C=C1C(C=CCNCC2(CCBr)C=CC=CO2)=CCc2ccc(F)cc21 | reductive_amination | C=C1C(C=CC=O)=CCc2ccc(F)cc21 + NCC1(CCBr)C=CC=CO1 | yes |
