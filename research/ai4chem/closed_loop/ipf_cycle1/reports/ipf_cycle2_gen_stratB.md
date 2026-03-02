# IPF Cycle 2 — Strategy B (training set enrichment)

## QMD evidence

- enrich_training_set: qmd://chem-scaffold-conditioned-gen/skill.md#L259-L310
- kinase enrichment data sources: qmd://chem-scaffold-conditioned-gen/skill.md#L327-L340
- cyclical annealing config: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- validate_generation hard gates: qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Step 1) enrich_training_set() — current training set composition

Training hinge-binder fraction (sample 5000, quick patterns) = **0.308**

enrich_training_set output:

```json
{
  "n_added": 0,
  "original_kinase_frac": 0.30787448194197753,
  "final_kinase_frac": 0.30787448194197753,
  "note": "Training set already meets kinase fraction target"
}
```

## Step 2) Kinase ligands retrieved (existing source: PubChem)

Fetched kinase_smiles_n = **1200**

```json
[
  {
    "source": "ChEMBL activity.json",
    "target": {
      "name": "EGFR",
      "target_chembl_id": "CHEMBL203"
    },
    "n_smiles": 357,
    "example": [
      "Cc1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",
      "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",
      "N#CC(C#N)=C(N)/C(C#N)=C/c1ccc(O)cc1"
    ]
  },
  {
    "source": "ChEMBL activity.json",
    "target": {
      "name": "BRAF",
      "target_chembl_id": "CHEMBL5145"
    },
    "n_smiles": 433,
    "example": [
      "CC(C)(C)c1nc(-c2ccc(F)cc2)c(-c2ccncc2)[nH]1",
      "Oc1ccc(-c2[nH]c(-c3ccccc3)nc2-c2ccncc2)cc1",
      "CC(C)(CNC(=O)Nc1ccc(Cl)cc1)c1nc(-c2ccc(Cl)c(O)c2)c(-c2ccncc2)[nH]1"
    ]
  },
  {
    "source": "ChEMBL activity.json",
    "target": {
      "name": "JAK2",
      "target_chembl_id": "CHEMBL2971"
    },
    "n_smiles": 413,
    "example": [
      "CC(C)(C)c1nc2c3ccc(F)cc3c3c(=O)[nH]ccc3c2[nH]1",
      "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ccc(F)cc32)c1C",
      "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4"
    ]
  }
]
```

## Step 3-4) Retrain VAE with cyclical annealing (same anti-collapse config)

Final diagnostics: {'kl_total': 5.65869665145874, 'kl_per_dim_mean': 0.17683427035808563, 'active_units': 32, 'total_dims': 32, 'au_ratio': 1.0, 'recon_accuracy': 0.0, 'diagnosis': 'healthy'}

## Step 5) validate_generation() — hard gates

```json
{
  "n_generated": 4807,
  "n_valid": 4807,
  "n_unique": 4807,
  "validity_pct": 100.0,
  "uniqueness_pct": 100.0,
  "novelty_pct": 100.0,
  "hinge_coverage_pct": 1.8,
  "train_hinge_pct": 47.7,
  "coverage_ratio": 0.038,
  "internal_diversity": 0.921,
  "gate_results": {
    "validity": {
      "value": 100.0,
      "threshold": 85.0,
      "pass": true
    },
    "hinge_coverage_ratio": {
      "value": 0.037522181717657276,
      "threshold": 0.3,
      "pass": false
    }
  },
  "all_gates_passed": false
}
```

Model saved? **False**
Generated SMILES: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_stratB_raw.csv`
