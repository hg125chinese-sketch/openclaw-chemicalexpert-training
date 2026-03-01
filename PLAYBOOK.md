# ChemicalExpert PLAYBOOK (how we work)

This playbook is aligned with the **truthbook verifiability standard**.

## Rule 0: Verifiable outputs (default)

When 陛下 asks for “查文档 / 给可验证指令 / 按 playbook 做”, I must:

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

## Standard QMD workflow (copyable)

```bash
cd /home/node/.openclaw/vault/openclaw-truthbook
QMD=/home/node/.openclaw/.npm-global/bin/qmd

# 1) search
$QMD search "<your query>" -c <collection> -n 10

# 2) open the top hit with line numbers
$QMD get qmd://<path> -l 220 | nl -ba
```

## Skill: chem-literature (must use for literature tasks)

When the user asks for papers, reading, or literature survey:

1) QMD search in `chem-literature`
2) QMD get the relevant file (at least `qmd://chem-literature/SKILL.md`)
3) Follow its workflow and note template

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

## Skill: chem-gnn (GNN molecular property prediction)
When using graph neural networks for molecular property prediction:
1) Always consult QMD collection: chem-gnn
2) Workflow:
   - /home/node/.openclaw/.npm-global/bin/qmd search "<query>" -c chem-gnn -n 10
   - /home/node/.openclaw/.npm-global/bin/qmd get qmd://chem-gnn/SKILL.md -l 300
3) Key rules:
   - Always compare GNN against RF + Morgan baseline (skill 2)
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
   - Use as filter for generated molecules (connect to skill 3)

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

## Skill: chem-3dgen (3D molecular generation)
When generating molecules in 3D space:
1) Consult: qmd search "<query>" -c chem-3dgen -n 10
2) Key rules:
   - Use when shape matters (pocket-conditioned design, scaffold hopping in 3D)
   - Bond assignment from point cloud requires validation with RDKit
   - Compare to skill 3 (SELFIES VAE) — use 3D gen only when it adds value
   - Report validity/uniqueness/novelty metrics

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

## Skill: chem-kinase-sar (kinase hinge binder SAR & coverage KPIs)
When targeting any kinase:
1) Consult: qmd search "<query>" -c chem-kinase-sar -n 10
2) Key rules:
   - Run pool coverage analysis BEFORE docking
   - Hard gate: ≥5% hinge binder coverage or stop and fix generator
   - Compare generated pool coverage to training set (ratio ≥0.5)
   - Classify final candidates by kinase inhibitor type (I/II/covalent)
   - This skill is kinase-specific — skip entirely for non-kinase targets

## Skill: chem-reactivity-safety (structural alerts & safety screening)
When screening candidates for safety liabilities:
1) Consult: qmd search "<query>" -c chem-reactivity-safety -n 10
2) Key rules:
   - Run AFTER ADMET (skill 6), BEFORE docking (skill 9)
   - 4 layers: PAINS → reactive metabolites → genotoxicity → chronic disease flags
   - "danger" = hard reject, "warning" = flag for review
   - Set is_chronic=True for chronic disease targets (IPF, diabetes, autoimmune)
   - Every flag must be auditable (SMARTS + matched atoms)
   - If >20% rejected → flag generation quality issue

## Skill: chem-protonation-tautomer (docking input standardization)
Before any docking:
1) Consult: qmd search "<query>" -c chem-protonation-tautomer -n 10
2) Key rules:
   - ALWAYS standardize (remove salts, neutralize, canonicalize) before docking
   - Enumerate protomers at pH 7.4 for molecules with ionizable groups
   - Dock ALL states, report best score + which state produced it
   - Especially critical for kinase targets (imidazole, aminopyrimidine near hinge)
   - Cap enumeration at 3-5 tautomers per protomer to avoid explosion
   - For large-pool screening: enumerate only top candidates after first-pass docking
