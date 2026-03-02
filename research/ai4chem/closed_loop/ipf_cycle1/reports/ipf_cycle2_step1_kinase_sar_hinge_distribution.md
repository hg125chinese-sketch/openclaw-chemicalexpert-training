# IPF Cycle 2 — Step 1 (chem-kinase-sar): hinge motif distribution on candidate pool

Input: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_final.csv`

Valid molecules: **500** (invalid dropped: 0)

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

- aminopyrimidine: 45 (9.00%)
- aminopyridine: 352 (70.40%)
- pyrazole: 107 (21.40%)
- purine_like: 0 (0.00%)

- ANY of above: 500 (100.00%)
- multi-match (>=2 patterns): 4 (0.80%)
