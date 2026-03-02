# IPF Cycle 2 — Generation diagnosis & repair (Step 1)

## QMD evidence (skill references)

- diagnose_vae definition: qmd://chem-scaffold-conditioned-gen/skill.md#L41-L88
- decision based on diagnosis: qmd://chem-scaffold-conditioned-gen/skill.md#L91-L106
- cyclical annealing trainer + free_bits/word_dropout: qmd://chem-scaffold-conditioned-gen/skill.md#L110-L215
- strategy selection guide A/B/C/D: qmd://chem-scaffold-conditioned-gen/skill.md#L651-L671
- validate_generation hard gate (hinge coverage ratio >=0.3, validity>=85): qmd://chem-scaffold-conditioned-gen/skill.md#L503-L586

## Inputs

- Training CSV: `research/ai4chem/closed_loop/ipf_cycle1/data/chembl_alk5_activity.csv`
- Cycle1 VAE ckpt: `research/ai4chem/closed_loop/ipf_cycle1/models/selfies_vae.pt`

## Step 1) diagnose_vae() — Cycle1 checkpoint

diagnosis=full_collapse | KL/dim=0.0007 | AU=0/32 (AU ratio=0.000) | recon_acc=0.000

## Step 2) Strategy selection

Selected strategy: **A**

## Step 3) Strategy A — retrain with cyclical annealing

Training config: n_epochs=30, n_cycles=4, free_bits=0.25, word_dropout=0.1

Final diagnostics (post-retrain):

diagnosis=healthy | KL/dim=0.1681 | AU=32/32 (AU ratio=1.000) | recon_acc=0.000

## Step 4) validate_generation() — hard gates

Generated unique valid SMILES written to: `research/ai4chem/closed_loop/ipf_cycle1/generated/generated_cycle2_raw.csv`

Metrics:

- validity_pct: 100.0
- uniqueness_pct: 100.0
- novelty_pct: 100.0
- hinge_coverage_pct: 1.0 (train_hinge_pct=47.7)
- hinge_binder_coverage_ratio: 0.021 (gate >= 0.3)
- internal_diversity: 0.904

Gate results:

- validity: value=100.000 threshold=85.0 pass=True
- hinge_coverage_ratio: value=0.021 threshold=0.3 pass=False

ALL hard gates passed? **False**
