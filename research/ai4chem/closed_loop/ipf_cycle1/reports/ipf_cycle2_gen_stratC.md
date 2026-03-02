# IPF Cycle 2 — Strategy C (Fragment-Conditioned VAE)

## QMD evidence

- Phase 4 (FragmentConditionedVAE): qmd://chem-scaffold-conditioned-gen/skill.md#Phase-4
- hinge binder SMARTS (skill 13): qmd://chem-scaffold-conditioned-gen/skill.md#Skill-13
- cyclical annealing config: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Step 1-3) Conditioning fragments + (molecule, fragment) pairs

Paired dataset size (hinge-containing only) = **286**

Example pairs (first 5):

```json
[
  {
    "molecule": "COCC(=O)Nc1cc(Oc2ccc(NC(=O)Nc3cc(C4(C)COC4)nn3-c3ccc(C)cc3)cc2)ccn1",
    "fragment": "CNc1cc(O)ccn1"
  },
  {
    "molecule": "C=CCN(C)c1nccc(Nc2cc(NC(=O)c3ccnc(N4CCOCC4)c3)ccc2C)n1",
    "fragment": "Cc1ccnc(N(C)C)c1"
  },
  {
    "molecule": "Cc1ccc(NC(=O)c2ccnc(N3CCOCC3)c2)cc1-n1c(C)nc2ccc(N3CCN(C)CC3)cc2c1=O",
    "fragment": "Cc1ccnc(N(C)C)c1"
  },
  {
    "molecule": "COC(=O)Nc1cc(Oc2ccc(NC(=O)Nc3cc(C4(C)CSC4)nn3-c3ccc(C)cc3)cc2)ccn1",
    "fragment": "CNc1cc(O)ccn1"
  },
  {
    "molecule": "Cc1ccc(-n2nc(C(C)(C)C)cc2NC(=O)Nc2ncc(CCc3ccnc(NC(=O)C(C)C)c3)s2)cc1",
    "fragment": "CNc1cc(C)ccn1"
  }
]
```

## Step 4) Train FragmentConditionedVAE (anti-collapse kept)

Final diagnostics: {'kl_total': 3.9413070678710938, 'kl_per_dim_mean': 0.12316584587097168, 'active_units': 32, 'total_dims': 32, 'au_ratio': 1.0}

## Step 5) Conditioned generation fragments

```json
{
  "aminopyrimidine": {
    "fragment": "Nc1nccnc1",
    "n_unique": 54
  },
  "aminopyridine": {
    "fragment": "Nc1ccccn1",
    "n_unique": 63
  },
  "pyrazole": {
    "fragment": "[nH]1ccnc1",
    "n_unique": 67
  }
}
```

## Step 6) validate_generation() — hard gates

```json
{
  "n_generated": 184,
  "n_valid": 184,
  "n_unique": 184,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 65.8,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 1.379,
  "internal_diversity": 0.759,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 1.3792128683980158,
      "threshold": 0.3,
      "pass": true
    }
  },
  "all_gates_passed": true
}
```

Model saved? **True**
Model path: `research/ai4chem/closed_loop/ipf_cycle1/models/selfies_vae_cycle2_stratC.pt`
Generated SMILES: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratC_raw.csv`
