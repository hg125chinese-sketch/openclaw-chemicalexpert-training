# IPF Cycle 2 — Strategy C (aux loss forcing fragment retention)

## QMD evidence

- Phase 4.1 architecture: qmd://chem-scaffold-conditioned-gen/skill.md#Phase-4
- cyclical annealing config: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Dataset

n_pairs (classified into 3 fragment classes) = **1272**

skipped_unclassified = **276**

class_counts:

```json
{
  "0": 0,
  "1": 664,
  "2": 608
}
```

example_pairs:

```json
[
  {
    "molecule": "CC(=O)Nc1cc(-c2c(-c3ccc(F)cc3)nc(S/C=C\\C(=O)O)n2CC(O)CO)ccn1",
    "fragment": "Nc1ccccn1",
    "class": 1
  },
  {
    "molecule": "CC(C)C(C)Nc1cc(-c2[nH]c(SCCO)nc2-c2ccc(F)cc2)ccn1",
    "fragment": "c1c[nH]cn1",
    "class": 2
  },
  {
    "molecule": "CC(C)(C)c1cccc(NC(=O)Nc2ccc(CCc3ccnc(NC(=O)C4CC4)c3)cc2F)c1",
    "fragment": "Nc1ccccn1",
    "class": 1
  },
  {
    "molecule": "CC(=O)Nc1cc(CCc2cnc(NC(=O)Nc3cc(C(C)(C)C)nn3-c3cccc(F)c3)s2)ccn1",
    "fragment": "Nc1ccccn1",
    "class": 1
  },
  {
    "molecule": "CC(=O)N1CCCC(C(=O)Nc2cc(-c3ccnc(Nc4ccccc4)c3)ccn2)C1",
    "fragment": "Nc1ccccn1",
    "class": 1
  }
]
```

## Training (recon + beta*KL + lambda_aux*aux)

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
  "kl_total": 4.157121658325195,
  "kl_per_dim_mean": 0.12991005182266235,
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
  "kl_total": 4.157121658325195,
  "kl_per_dim_mean": 0.12991005182266235,
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
    "n_unique": 19598
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 19575
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 19418
  }
}
```

## validate_generation()

```json
{
  "n_generated": 58297,
  "n_valid": 58297,
  "n_unique": 58297,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 1.9,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 0.04,
  "internal_diversity": 0.847,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 0.040185624342279694,
      "threshold": 0.3,
      "pass": false
    }
  },
  "all_gates_passed": false
}
```

Model saved? **False**
Generated SMILES: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_auxloss_raw.csv`
