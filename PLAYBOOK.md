# ChemicalExpert PLAYBOOK (how we work)

This playbook is aligned with the **truthbook verifiability standard**.

## Rule 0: Verifiable outputs (default)

When the user asks for "check docs / give verifiable instructions / follow playbook", I must:

1) **Use QMD first** (pick the most relevant collection)
2) **Show evidence** in the final answer:
   - (1) the qmd hit line (qmd://path + line)
   - (2) a quoted excerpt **with line numbers**
   - (3) exact commands to run (copy/paste)

If no collection covers the task, I must say so and propose adding a new entry/collection.

## Canonical QMD binary

Always use the absolute path:

```bash
QMD=/home/node/.openclaw/.npm-global/bin/qmd
```

## Collection picker

- **chem-literature** → chemistry/AI4Chem paper search + deep reads (arXiv / Semantic Scholar / PubChem)
- **agent-browser** → browser automation / scraping / form filling / UI regression

## Skill routing table (skills 1–17)

Principle: when the user says "use skill X", route to the corresponding **QMD collection** and follow its `SKILL.md`.

| # | Skill / QMD collection | Use when | Typical outputs |
|---:|---|---|---|
| 1 | **chem-qsar** | Molecular property / QSAR modeling, data leakage audit, split/metrics/baseline | Scaffold split results, RF baseline, comparison tables, reproducible scripts |
| 2 | **chem-gnn** | GNN-based property prediction or QSAR improvement | GCN/GAT/MPNN training & comparison, metrics + param count |
| 3 | **chem-admet** | Rule + ML ADMET profiling / filtering (QED/SA/Lipinski/Veber etc.) | Per-molecule ADMET profile, filter thresholds, TopN summary |
| 4 | **chem-retrosynthesis** | Retrosynthesis, synthesizability, route templates, forward-validation | routes.md, step conditions, SA score, failure reasons |
| 5 | **chem-rxn-conditions** | Reaction condition suggestion / yield range, linked to retro steps | Condition annotations, route yield, compatibility warnings |
| 6 | **chem-llm** | LLM-assisted chemistry tasks (must RDKit/computational verify) | "LLM-generated + verified" structured output, traceable prompts |
| 7 | **chem-molgen** | Molecular generation (VAE/Transformer/Diffusion), string metrics, KL collapse diagnosis | validity/uniqueness/novelty/diversity, AU, training/sampling scripts |
| 8 | **chem-3dgen** | 3D generation / shape constraints / pocket conditioning | 3D validity, conformation/bond inference checks, comparison results |
| 9 | **chem-docking** | Structure-based docking (Vina) / redocking validation / box definition | report.json, dock_scores.csv, RMSD gate, version & input records |
| 10 | **chem-mlff** | Conformer optimization / strain energy / geometry ranking (MMFF vs MACE etc.) | Pre/post optimization energy, strain energy, conformer ranking |
| 11 | **chem-experiment** | DMTA cycle orchestration (multi-skill chaining, multi-objective gates) | Cycle report, decision rationale, next-step plan |
| 12 | **chem-literature** | Literature search / deep reading / evidence line numbers, reproduction chain | Reading notes, evidence citations, reproduction commands, repo diffs |
| 13 | **chem-kinase-sar** | Kinase hinge binder motif / KPI, coverage audit | Motif distribution, coverage_ratio, KPI report |
| 14 | **chem-reactivity-safety** | Reactivity / toxicity / chronic disease safety hard reject (denylist) | Per-molecule flags, hard reject stats, survivors.csv |
| 15 | **chem-protonation-tautomer** | pH-dependent protomer/tautomer, multi-state preparation & enumeration | Multi-state enumeration, expansion ratio, best-state records |
| 16 | **chem-docking-interactions** | Docking pose interaction fingerprints, hinge H-bond validation, interaction rerank | Hinge yes/no + details, key residue coverage, rerank table |
| 17 | **chem-scaffold-conditioned-gen** | Scaffold/fragment-conditioned generation (data audit, latent conditioning, constrained sampling) | Conditioning hit rate, KPI gate, generation set & diagnosis report |

## Standard QMD workflow (copyable)

```bash
cd /home/node/.openclaw/vault/openclaw-truthbook
QMD=/home/node/.openclaw/.npm-global/bin/qmd

# 1) search
$QMD search "<your query>" -c <collection> -n 10

# 2) open the top hit with line numbers
$QMD get qmd://<path> -l 220 | nl -ba
```

## Skill: chem-qsar (QSAR / molecular property prediction)
When building or evaluating molecular property prediction models:
1) Always consult QMD collection: chem-qsar
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-qsar -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-qsar/SKILL.md -l 300
3) Key rules from this skill:
   - Scaffold split by default (random split only as sanity check)
   - Always run RF baseline before neural networks
   - Report RMSE/std(y), not just RMSE
   - Scaler fit on train only

## Skill: chem-gnn (GNN molecular property prediction)
When using graph neural networks for molecular property prediction:
1) Always consult QMD collection: chem-gnn
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-gnn -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-gnn/SKILL.md -l 300
3) Key rules:
   - Always compare GNN against RF + Morgan baseline (skill 1)
   - Start with GCN, only go GAT/MPNN if GCN plateaus
   - Scaffold split mandatory
   - Report parameter count alongside metrics

## Skill: chem-admet (ADMET prediction)
When predicting drug-likeness or pharmacokinetic properties:
1) Always consult QMD collection: chem-admet
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-admet -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-admet/SKILL.md -l 300
3) Key rules:
   - Rules first (Lipinski/Veber/QED), ML second
   - Always report multi-endpoint profiles, not single predictions
   - Use as filter for generated molecules (connect to skill 7)

## Skill: chem-retrosynthesis (retrosynthetic analysis)
When analyzing synthesis routes or evaluating synthesizability:
1) Always consult QMD collection: chem-retrosynthesis
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-retrosynthesis -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-retrosynthesis/SKILL.md -l 300
3) Key rules:
   - Always analyze FGs and suggest disconnections before applying templates
   - Forward-validate every retro step
   - Report SA Score for all targets
   - Filter generative model output by synthesizability before reporting

## Skill: chem-rxn-conditions (reaction condition prediction)
When recommending reaction conditions (solvent, catalyst, temp, yield):
1) Always consult QMD collection: chem-rxn-conditions
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-rxn-conditions -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-rxn-conditions/SKILL.md -l 300
3) Key rules:
   - Annotate every retro step (skill 4) with conditions
   - Yield is a range, not a point estimate
   - Substrate-aware adjustments (FG compatibility, steric)
   - Calculate overall route yield (multiplicative)

## Skill: chem-llm (chemical LLM applications)
When using LLM reasoning for chemistry tasks:
1) Always consult QMD collection: chem-llm
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-llm -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-llm/SKILL.md -l 300
3) Key rules:
   - LLMs reason, they don't compute — always validate with RDKit
   - Every LLM-generated SMILES must be RDKit-validated
   - Tag all LLM outputs as "[LLM-generated, pending validation]"
   - Use prompt templates from skill for consistent quality

## Skill: chem-molgen (molecular generative models)
When training or evaluating molecular generative models (VAE/Transformer/Diffusion):
1) Always consult QMD collection: chem-molgen
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-molgen -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-molgen/SKILL.md -l 300
3) Key rules:
   - KL collapse is the default — always diagnose active units
   - Use cyclical annealing as first line defense
   - Full evaluation: validity + uniqueness + novelty + diversity + property distributions
   - SELFIES for VAE, SMILES for pretrained Transformers

## Skill: chem-3dgen (3D molecular generation)
When generating molecules in 3D space:
1) Consult: qmd search "<query>" -c chem-3dgen -n 10
2) Key rules:
   - Use when shape matters (pocket-conditioned design, scaffold hopping in 3D)
   - Bond assignment from point cloud requires validation with RDKit
   - Compare to skill 7 (SELFIES VAE) — use 3D gen only when it adds value
   - Report validity/uniqueness/novelty metrics

## Skill: chem-docking (molecular docking)
When doing structure-based docking / redocking validation:
1) Always consult QMD collection: chem-docking
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-docking -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-docking/<doc> -l 300
3) Key rules:
   - Start from a known PDB with co-crystal ligand
   - Do receptor prep → define box from ligand → **redocking** first
   - Only proceed to dock new molecules if redocking meets RMSD threshold (e.g., <2 Å)
   - Log all software versions and inputs; treat docking scores as comparative, not absolute

## Skill: chem-mlff (ML force fields)
When optimizing geometry or ranking conformers:
1) Consult: qmd search "<query>" -c chem-mlff -n 10
2) Key rules:
   - MACE-OFF for organics, MACE-MP-0 for everything else
   - Always compare to MMFF baseline
   - Strain energy: < 3 kcal/mol acceptable, > 10 suspicious
   - Units: 1 eV = 23.06 kcal/mol — always label clearly

## Skill: chem-experiment (DMTA cycle planning)
When orchestrating multi-skill workflows:
1) Consult: qmd search "<query>" -c chem-experiment -n 10
2) Key rules:
   - Plan before computing — get user approval first
   - Multi-objective scoring (activity + QED + SA + dock + consistency)
   - Gates are non-negotiable (SA≤6, QED>0.3, |RF-GNN|≤1.0)
   - Each cycle must produce an actionable decision for the next
   - Log which skills used, in what order, why

## Skill: chem-literature (literature search & deep reading)
When the user asks for papers, reading, or literature survey:
1) QMD search in `chem-literature`
2) QMD get the relevant file (at least `qmd://chem-literature/SKILL.md`)
3) Follow its workflow and note template
4) Key rules:
   - Always cite with evidence line numbers
   - Reproduce key claims with code when possible
   - Track repo URLs and version hashes for reproducibility

## Skill: chem-kinase-sar (kinase SAR & hinge binder analysis)
When analyzing kinase structure-activity relationships or hinge binder coverage:
1) Always consult QMD collection: chem-kinase-sar
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-kinase-sar -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-kinase-sar/SKILL.md -l 300
3) Key rules:
   - Hinge H-bond is NON-NEGOTIABLE for Type I kinase inhibitors
   - Always audit hinge binder motif distribution (aminopyrimidine, aminopyridine, pyrazole, purine-like)
   - Report coverage_ratio = generated_hinge_pct / training_hinge_pct
   - Hard gate: coverage_ratio >= 0.3 before proceeding to docking

## Skill: chem-reactivity-safety (reactivity & safety screening)
When screening molecules for reactive functional groups or project-specific safety concerns:
1) Always consult QMD collection: chem-reactivity-safety
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-reactivity-safety -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-reactivity-safety/SKILL.md -l 300
3) Key rules:
   - Project denylist = hard reject (no exceptions for chronic disease indications)
   - Screen BEFORE docking (fail-fast saves compute)
   - Report per-molecule flags and aggregate statistics
   - Output survivors.csv for downstream pipeline consumption

## Skill: chem-protonation-tautomer (protonation & tautomer handling)
When preparing molecules for pH-dependent calculations or multi-state docking:
1) Always consult QMD collection: chem-protonation-tautomer
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-protonation-tautomer -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-protonation-tautomer/SKILL.md -l 300
3) Key rules:
   - Enumerate protomers/tautomers at target pH (default: 7.4) before docking
   - Dock ALL states, report best score per parent molecule
   - Track which protonation state gave the best pose (critical for SAR interpretation)
   - Apply only to top candidates (e.g., Top100) — too expensive for full libraries

## Skill: chem-docking-interactions (interaction fingerprints & pose QC)
When validating docking poses or analyzing protein-ligand interactions:
1) Always consult QMD collection: chem-docking-interactions
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-docking-interactions -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-docking-interactions/SKILL.md -l 300
3) Key rules:
   - Always profile co-crystal ligand first as ground truth reference
   - Hinge H-bond is mandatory for Type I kinase inhibitors — no H-bond = hard reject
   - Use ProLIF for batch analysis, PLIP for detailed single-pose inspection
   - interaction_rerank: composite score = f(Vina_normalized, interaction_score)
   - Report Top5 with: SMILES, Vina, interaction_score, hinge_hbond (YES/NO + residue + distance), key_residue_coverage

## Skill: chem-scaffold-conditioned-gen (conditioned molecular generation)
When addressing generation quality issues or conditioning on specific scaffolds/fragments:
1) Always consult QMD collection: chem-scaffold-conditioned-gen
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-scaffold-conditioned-gen -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-scaffold-conditioned-gen/SKILL.md -l 300
3) Key rules:
   - Always diagnose VAE health first (diagnose_vae → KL/dim, AU, AU ratio)
   - Strategy selection: full_collapse → Strategy A (cyclical annealing); healthy but low coverage → Strategy B (enrich) or C (conditioning); any strategy → Strategy D (rejection sampling) as pragmatic fallback
   - Hard gates: validity >= 85%, hinge_binder_coverage_ratio >= 0.3, no full_collapse
   - Change ONLY ONE variable per experiment — commit ALL results (success + failure)
   - Training set audit: verify fragment class distribution before conditioning experiments

## Browser automation

Prefer `agent-browser` (CLI) for web automation:

```bash
npx agent-browser open "https://example.com"
npx agent-browser snapshot -i
npx agent-browser screenshot --full "page.png"
npx agent-browser close
```

## Scientific computing environment

Use the conda environment Python for scientific work:

```bash
/opt/conda/envs/chem/bin/python
```

System python is not for scientific stacks:

```bash
/usr/bin/python3
```

(See `TOOLS.md` for the full environment/package list and GPU notes.)
