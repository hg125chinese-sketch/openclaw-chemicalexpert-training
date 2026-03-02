# IPF Cycle 2 — Strategy C validation (robustness)

## QMD evidence

- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Setup

- Model: `research/ai4chem/closed_loop/ipf_cycle1/models/selfies_vae_cycle2_stratC.pt`
- n_per_fragment: 20000
- Total attempted samples: 3 × 20000 = 60000

## Experiment A (scale-up, prefix=True)

- unique_smiles_csv: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_valA.csv`

Per-condition unique counts:

```json
{
  "aminopyrimidine": {
    "fragment": "Nc1nccnc1",
    "n_unique": 540
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 568
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 561
  }
}
```

validate_generation():

```json
{
  "n_generated": 1669,
  "n_valid": 1669,
  "n_unique": 1669,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 61.4,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 1.287,
  "internal_diversity": 0.835,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 1.2867890992878426,
      "threshold": 0.3,
      "pass": true
    }
  },
  "all_gates_passed": true
}
```

## Experiment B (scale-up, prefix=False)

- unique_smiles_csv: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_valB.csv`

Per-condition unique counts:

```json
{
  "aminopyrimidine": {
    "fragment": "Nc1nccnc1",
    "n_unique": 15868
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 16006
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 15626
  }
}
```

validate_generation():

```json
{
  "n_generated": 45477,
  "n_valid": 45477,
  "n_unique": 45477,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 0.1,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 0.002,
  "internal_diversity": 0.881,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 0.002259789704154604,
      "threshold": 0.3,
      "pass": false
    }
  },
  "all_gates_passed": false
}
```

## Decision

**A_PASS_B_FAIL: prefix trivial solution → need true latent conditioning (retrain w/ concat(z, frag_embed))**
