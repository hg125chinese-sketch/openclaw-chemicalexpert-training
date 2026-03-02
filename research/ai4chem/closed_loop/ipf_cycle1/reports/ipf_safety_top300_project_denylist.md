# IPF Cycle 1 — project safety denylist backtest (Top300)

Context: chronic disease (IPF), so selected alerts are upgraded to HARD REJECT.

Upgraded-to-hard-reject alerts:
- `acrylamide_chronic`
- `aldehyde`
- `chloroacetamide_chronic`
- `michael_acceptor`
- `vinyl_sulfone_chronic`

Top300 (pre-dock) n_valid=300, hard_reject(danger)=88 → **29.3%**

## Most frequent alerts (all severities)

| alert | count | pct |
|---|---:|---:|
| reactive:ipf_hard_michael_acceptor | 52 | 17.3% |
| reactive:michael_acceptor | 52 | 17.3% |
| reactive:n_n_single_bond | 30 | 10.0% |
| chronic:n_n_bond_chronic | 30 | 10.0% |
| chronic:acrylamide_chronic | 28 | 9.3% |
| reactive:ipf_hard_aldehyde | 27 | 9.0% |
| reactive:aldehyde | 27 | 9.0% |
| chronic:high_logp_chronic | 23 | 7.7% |
| reactive:hydrazone | 15 | 5.0% |
| chronic:hydrazone_chronic | 15 | 5.0% |
| chronic:ester_prodrug_risk | 14 | 4.7% |
| reactive:alkyl_halide | 8 | 2.7% |
| reactive:thiol_free | 7 | 2.3% |
| reactive:acyl_halide | 5 | 1.7% |
| reactive:nitroso | 4 | 1.3% |

## Impact on current Top5 (from top300_scored.csv)

hard_reject_in_top5 = **0 / 5**

| # | smiles | severity | upgraded_hits |
|---:|---|---|---|
| 1 | `CSC(=O)c1cccc(COc2cc(F)ccc2C)c1` | clean |  |
| 2 | `Cc1cccc(Nc2ccccc2Cl)c1` | clean |  |
| 3 | `CCNC(=O)c1cccc(Cc2ccccn2)c1` | clean |  |
| 4 | `Cc1ccc(F)cc1CC=Cc1ccccc1Cl` | warning |  |
| 5 | `CC(=O)Nc1cccnc1NC=CNc1cc(F)ccc1C` | clean |  |
