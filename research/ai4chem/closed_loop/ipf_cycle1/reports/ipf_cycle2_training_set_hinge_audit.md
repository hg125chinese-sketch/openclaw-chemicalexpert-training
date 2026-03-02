# IPF Cycle 2 — Training set hinge motif audit

Input: `research/ai4chem/closed_loop/ipf_cycle1/data/chembl_alk5_activity.csv`

Valid molecules: **5067** (invalid SMILES dropped: 0)

## Patterns (SMARTS)

```json
{
  "aminopyrimidine": "c1ncnc(N)c1",
  "aminopyridine": "c1ccnc(N)c1",
  "pyrazole": "[nH]1ccnc1",
  "purine_like": "c1ncnc2[nH]cnc12"
}
```

## Counts & fractions

- aminopyrimidine: 277 (5.47%)
- aminopyridine: 933 (18.41%)
- pyrazole: 622 (12.28%)
- purine_like: 4 (0.08%)

- ANY of above: 1560 (30.79%)

## Notes

- ‘ANY of above’ is the union over the 4 patterns; this is the composition of the hinge-like KPI used in Cycle 2 gates.
- Next Step 7 decision: if aminopyrimidine count is near-zero, conditioning on aminopyrimidine is out-of-distribution for the training set and we should switch fragments or enrich data.
