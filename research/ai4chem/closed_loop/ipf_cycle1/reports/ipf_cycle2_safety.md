# IPF Cycle 2 — Step 2 (chem-reactivity-safety): IPF chronic project hard reject screen

Input: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_final.csv`

n_input=500 | n_hard_reject=204 | n_survivors=296

## Upgraded alerts → HARD REJECT

```json
[
  "acrylamide_chronic",
  "aldehyde",
  "chloroacetamide_chronic",
  "michael_acceptor",
  "vinyl_sulfone_chronic"
]
```

## Most common flags (top15)

```json
[
  [
    "michael_acceptor",
    118
  ],
  [
    "high_logp",
    89
  ],
  [
    "n_n_single_bond",
    85
  ],
  [
    "n_n_bond_chronic",
    85
  ],
  [
    "hydrazone",
    57
  ],
  [
    "hydrazone_chronic",
    57
  ],
  [
    "acrylamide_chronic",
    55
  ],
  [
    "thiol_free",
    25
  ],
  [
    "aldehyde",
    25
  ],
  [
    "alkyl_halide",
    20
  ],
  [
    "ester_prodrug_risk",
    14
  ],
  [
    "aniline_unsubst",
    13
  ],
  [
    "pains",
    10
  ],
  [
    "hydrazine",
    3
  ],
  [
    "aziridine",
    3
  ]
]
```

## Outputs

- per-molecule: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle2_step2_safety_per_molecule.csv`
- survivors: `research/ai4chem/closed_loop/ipf_cycle1/reports/ipf_cycle2_step2_safety_survivors.csv`
