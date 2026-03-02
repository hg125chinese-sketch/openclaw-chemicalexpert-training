# IPF Cycle 2 — Strategy C (true latent conditioning, no prefix)

## QMD evidence

- Phase 4.1 architecture: qmd://chem-scaffold-conditioned-gen/skill.md#Phase-4
- cyclical annealing config: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Dataset (Step 3 pairs reused)

n_pairs = **1548**

Example pairs:

```json
[
  {
    "molecule": "CSc1nc(-c2ccc(F)cc2)c(-c2ccnc(NC(C)=O)c2)n1CC(=O)N1CCOCC1",
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "CC(=O)Nc1cc(CCc2cnc(NC(=O)Nc3cc(C(C)(C)C)on3)s2)ccn1",
    "fragment": "Nc1ccccn1"
  },
  {
    "molecule": "O=S(=O)(NCCNc1cc(-c2[nH]c(-c3ccccc3)nc2-c2cccc(O)c2)ccn1)c1ccccc1",
    "fragment": "c1c[nH]cn1"
  },
  {
    "molecule": "CCN(C(=O)c1ccccc1)c1ccnc(N[C@@H](C)c2ccccc2)n1",
    "fragment": "Nc1ccncn1"
  },
  {
    "molecule": "CC(=O)Nc1cc(CCc2ccc(NC(=O)Nc3cccc(C(C)(C)C)c3)c(F)c2)ccn1",
    "fragment": "Nc1ccccn1"
  }
]
```

## Training (anti-collapse preserved)

Final diagnostics:

```json
{
  "kl_total": 4.518337249755859,
  "kl_per_dim_mean": 0.1411980390548706,
  "active_units": 32,
  "total_dims": 32,
  "au_ratio": 1.0,
  "diagnosis": "healthy"
}
```

Post-training diagnose_vae (re-collapse check):

```json
{
  "kl_total": 4.518337249755859,
  "kl_per_dim_mean": 0.1411980390548706,
  "active_units": 32,
  "total_dims": 32,
  "au_ratio": 1.0,
  "diagnosis": "healthy"
}
```

## Generation (prefix=False; conditioning via concat(z, frag_embed))

Per-condition unique counts:

```json
{
  "aminopyrimidine": {
    "fragment": "Nc1nccnc1",
    "n_unique": 19668
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 19645
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 19534
  }
}
```

## validate_generation() (hard gates)

```json
{
  "n_generated": 58545,
  "n_valid": 58545,
  "n_unique": 58545,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 2.6,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 0.054,
  "internal_diversity": 0.863,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 0.053843455473074404,
      "threshold": 0.3,
      "pass": false
    }
  },
  "all_gates_passed": false
}
```

Model saved? **False**
Generated SMILES: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_latent_raw.csv`
