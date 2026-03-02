# IPF Cycle 2 — Strategy C aux loss v2 (fixed fragment inputs + SMARTS labeling)

## QMD evidence

- Phase 4.1 architecture: qmd://chem-scaffold-conditioned-gen/skill.md#Phase-4
- cyclical annealing config: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Labeling (molecule SMARTS → class; multi-label expands to multiple pairs)

```json
{
  "aminopyrimidine": {
    "id": 0,
    "smarts": "c1ncnc(N)c1",
    "smiles": "Nc1nccnc1"
  },
  "aminopyridine": {
    "id": 1,
    "smarts": "c1ccnc(N)c1",
    "smiles": "Nc1ccccn1"
  },
  "pyrazole": {
    "id": 2,
    "smarts": "[nH]1ccnc1",
    "smiles": "[nH]1ccnc1"
  }
}
```

## Dataset summary

n_pairs = **1832** (capped if necessary)

class_counts:

```json
{
  "0": 277,
  "1": 933,
  "2": 622
}
```

examples:

```json
[
  {
    "molecule": "[O-][n+]1ccc2c(-c3ccccc3Cl)cc(NC3CCNCC3)nc2c1-c1c(Cl)cccc1Cl",
    "class": 1,
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "[O-][n+]1ccc2c(-c3ccccc3Cl)cc(N3CCNCC3)nc2c1-c1c(Cl)cccc1Cl",
    "class": 1,
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "[O-][n+]1ccc2c(-c3ccc(F)cc3F)cc(NC3CCNCC3)nc2c1-c1c(F)cccc1F",
    "class": 1,
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "Fc1ccc(Nc2ccc3c(-c4ccc(F)cc4F)nncc3n2)c(F)c1",
    "class": 1,
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "Cc1ccc(NC(=O)c2ccnc(N3CCOCC3)c2)cc1Nc1ncnn2cc(C(=O)N[C@@H](C)c3ccccc3)c(C)c12",
    "class": 1,
    "fragment": "Nc1ccccn1"
  }
]
```

## Training

config:

```json
{
  "n_epochs": 60,
  "n_cycles": 4,
  "free_bits": 0.25,
  "word_dropout": 0.1,
  "lambda_aux": 0.5,
  "z_dropout": 0.2
}
```

final_diagnostics:

```json
{
  "kl_total": 4.95093297958374,
  "kl_per_dim_mean": 0.15471665561199188,
  "active_units": 32,
  "total_dims": 32,
  "au_ratio": 1.0,
  "diagnosis": "healthy"
}
```

final_aux_acc (val) = **1.000**

post diagnose_vae:

```json
{
  "kl_total": 4.95093297958374,
  "kl_per_dim_mean": 0.15471665561199188,
  "active_units": 32,
  "total_dims": 32,
  "au_ratio": 1.0,
  "diagnosis": "healthy"
}
```

## Generation (prefix=False; 20000 per fragment)

per-condition unique:

```json
{
  "aminopyrimidine": {
    "fragment": "Nc1nccnc1",
    "n_unique": 19682
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 19661
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 19465
  }
}
```

## validate_generation()

```json
{
  "n_generated": 58352,
  "n_valid": 58352,
  "n_unique": 58352,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 4.5,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 0.095,
  "internal_diversity": 0.869,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 0.09488814010964239,
      "threshold": 0.3,
      "pass": false
    }
  },
  "all_gates_passed": false
}
```

Model saved? **False**
Generated SMILES: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_auxloss_v2_raw.csv`
