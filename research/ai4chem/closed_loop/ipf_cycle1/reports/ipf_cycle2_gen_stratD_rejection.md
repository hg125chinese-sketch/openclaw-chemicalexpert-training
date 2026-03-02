# IPF Cycle 2 — Strategy D (rejection sampling for hinge SMARTS)

## QMD evidence

- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Setup

- Step1 checkpoint: `research/ai4chem/closed_loop/ipf_cycle1/models/selfies_vae_cycle2.pt`
- Target: n_target=500 (ANY of required SMARTS)
- max_attempts=500000

## Required SMARTS (ANY-match acceptance)

```json
{
  "aminopyrimidine": "c1ncnc(N)c1",
  "aminopyridine": "c1ccnc(N)c1",
  "pyrazole": "[nH]1ccnc1",
  "purine_like": "c1ncnc2[nH]cnc12"
}
```

## Rejection sampling efficiency

```json
{
  "attempts": 35000,
  "decoded": 32472,
  "valid": 32346,
  "n_accepted": 500,
  "accept_rate_vs_valid": 0.015457861868546342,
  "accept_rate_vs_attempts": 0.014285714285714285
}
```

accept_rate_vs_valid = **1.546%** → **POOR** (>5% good, 0.5–5% poor, <0.5% terrible)

## validate_generation() on accepted pool

```json
{
  "n_generated": 500,
  "n_valid": 500,
  "n_unique": 500,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 75.0,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 1.573,
  "internal_diversity": 0.882,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 1.5729865771812082,
      "threshold": 0.3,
      "pass": true
    }
  },
  "all_gates_passed": true
}
```

## Output

- Pool CSV: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_final.csv`
- Hard gates passed? **True**
