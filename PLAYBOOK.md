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

## Skill routing table (skills 1–17)

原则：用户说“用 skill X”时，按下表路由到对应 **QMD collection**，并遵循该 collection 的 `SKILL.md`。

| # | Skill / QMD collection | Use when | Typical outputs |
|---:|---|---|---|
| 1 | **chem-qsar** | 分子性质/QSAR 建模、数据泄漏审计、split/metrics/baseline | scaffold split 结果、RF baseline、对比表、可复现实验脚本 |
| 2 | **chem-gnn** | 用 GNN 做性质预测或 QSAR 提升 | GCN/GAT/MPNN 训练与对照、指标+参数量 |
| 3 | **chem-admet** | 规则+ML 的 ADMET 画像/过滤（QED/SA/Lipinski/Veber 等） | per-molecule ADMET profile、过滤门槛、TopN 汇总 |
| 4 | **chem-retrosynthesis** | 逆合成、可合成性、路线模板、forward-validate | routes.md、每步条件、SA 分数、失败原因 |
| 5 | **chem-rxn-conditions** | 反应条件建议/产率范围、与 retro 步骤联动 | 条件注释、route yield、兼容性警告 |
| 6 | **chem-llm** | LLM 辅助化学任务（必须 RDKit/计算校验） | “LLM 生成+校验”结构化输出、可追溯提示词 |
| 7 | **chem-molgen** | 分子生成（VAE/Transformer/Diffusion）、string metrics、KL collapse 诊断 | validity/uniqueness/novelty/diversity、AU、训练/采样脚本 |
| 8 | **chem-3dgen** | 3D 生成/形状约束/口袋条件化 | 3D validity、构型/键推断校验、对照结果 |
| 9 | **chem-docking** | 结构对接（Vina）/redocking 验证/box 定义 | report.json、dock_scores.csv、RMSD gate、版本与输入记录 |
| 10 | **chem-mlff** | 构象优化/应变能/几何排序（MMFF vs MACE 等） | 优化前后能量、应变能、构象排序 |
| 11 | **chem-experiment** | DMTA 周期编排（多技能串联、multi-objective gates） | cycle report、决策理由、下一步计划 |
| 12 | **chem-literature** | 文献检索/深读/证据行号、复现链路 | reading notes、证据引用、复现命令、repo 差异 |
| 13 | **chem-kinase-sar** | kinase hinge binder motif/KPI、覆盖率审计 | motif 分布、coverage_ratio、KPI 报告 |
| 14 | **chem-reactivity-safety** | 反应性/毒性/慢病安全 hard reject（denylist） | per-molecule flags、hard reject 统计、survivors.csv |
| 15 | **chem-protonation-tautomer** | pH 相关 protomer/tautomer、多状态准备与展开 | multi-state 枚举、展开比例、best-state 记录 |
| 16 | **chem-docking-interactions** | docking pose 的相互作用指纹、hinge H-bond 验证、interaction rerank | hinge yes/no+细节、key residue coverage、rerank 表 |
| 17 | **chem-scaffold-conditioned-gen** | scaffold/fragment 条件生成（含数据审计、latent conditioning、约束采样等） | 条件命中率、KPI gate、生成集与诊断报告 |

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
