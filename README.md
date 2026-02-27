# Teaching AI Agents to Discover Drugs: A Systematic Training Methodology

**How I trained an LLM agent to autonomously run drug discovery workflows — and what it taught itself along the way.**

---

## TL;DR

I trained an AI agent ("ChemicalExpert") with 12 specialized chemistry skills using a systematic 7-step methodology. The agent learned to autonomously plan and execute Design-Make-Test-Analyze (DMTA) cycles for drug discovery. When given an open-ended task — "find anti-fibrosis drug candidates" — it independently selected a target, ran the full pipeline, and then *diagnosed its own output quality problems* without being asked.

This post describes the methodology, the skills, and what happened when the agent was tested on a target it had never seen before.

---

## The Problem

Large language models know a lot about chemistry. But knowing facts and *doing chemistry* are different things. An LLM can explain what QSAR is, but can it train a random forest on ChEMBL data, validate with scaffold splits, and tell you when the model is unreliable?

I wanted an agent that could:
- Run computational chemistry workflows end-to-end (data → models → molecules → evaluation)
- Know *when* to use each tool and *when not to*
- Catch its own mistakes before I catch them
- Work on new targets without being retrained

The platform is [OpenClaw](https://github.com/openclaw-ai/openclaw), an open-source agent framework. The underlying model is GPT-5.2 via OpenAI's Codex API. The agent runs in a Docker container with GPU support on NixOS/WSL2, with access to conda environments, RDKit, PyTorch, and standard scientific Python.

## The 7-Step Training Methodology

Through trial and error training two agents (a debugging agent and ChemicalExpert), I converged on a 7-step process:

### Step 1: Pre-Assessment
Before writing any skill, ask the agent diagnostic questions across four dimensions: current knowledge, methodology, tool capabilities, and learning preferences. Let it tell you what it knows — don't assume.

### Step 2: Deep Probing
Design pointed questions targeting suspected gaps. The goal is to expose the difference between "knows about" and "can do."

### Step 3: Gap Analysis
Compare pre-assessment against probing results. Only fill gaps — don't reteach what the agent already knows.

### Step 4: Skill Design
Each skill is a structured document (~200-500 lines) with decision boundaries, runnable code examples, philosophy sections, and failure modes.

**Key insight**: Skills should be structured around "when to do what" rather than "how to do everything."

### Step 5: Guided Practice
Give the agent a concrete task that requires the new skill. Watch it work. Intervene only when necessary.

### Step 6: Behavioral Correction
When the agent makes mistakes, challenge it to find the error itself rather than giving the answer directly.

### Step 7: Self-Review
Ask the agent to write its own skill review. Compare against your observations. Discrepancies reveal shallow understanding.

## The 12 Skills

| # | Skill | Core Capability |
|---|-------|----------------|
| 1 | Literature Review | Search strategies, source verification, evidence tables |
| 2 | QSAR Modeling | Scaffold splits, RF baseline, descriptor selection |
| 3 | Molecular Generation | SELFIES VAE, KL collapse detection, novelty metrics |
| 4 | Retrosynthetic Analysis | Template matching, forward validation, route scoring |
| 5 | Graph Neural Networks | GCN for molecular properties, scaffold generalization vs RF |
| 6 | ADMET Prediction | Rules-first (Lipinski/Veber), QED gates, chronic disease flags |
| 7 | Reaction Conditions | Yield prediction, reagent selection, total yield calculation |
| 8 | Chemical LLM | When to use LLM reasoning vs computation |
| 9 | Molecular Docking | Receptor prep, redocking validation, virtual screening |
| 10 | 3D Generation | E(3)-equivariant diffusion, point cloud to molecule |
| 11 | ML Force Fields | MACE-OFF geometry optimization, conformer ranking |
| 12 | Experiment Planning | DMTA cycle orchestration, multi-objective scoring |

Each skill document is in `skills/*/SKILL.md`.

## Skill Design Principles

1. **Philosophy over procedure.** "Always compare MACE to MMFF baseline" beats 50 lines of API docs.
2. **Decision trees over checklists.** "If X, do A; if Y, do B" matches real drug discovery decisions.
3. **Failure modes are non-negotiable.** Every skill includes what failure looks like.
4. **Gates, not guidelines.** Hard thresholds (SA ≤ 6, QED > 0.3) that must be passed.
5. **One sentence per skill.** If you can't state the core rule in one sentence, the skill is too broad.

## The Migration Test: IPF/ALK5

After validating on DRD2, I gave the agent an open-ended task with no guidance:

> "Find anti-fibrosis (IPF) drug candidates. Pirfenidone and nintedanib have limited efficacy."

### What the agent did autonomously
- **Selected ALK5/TGFBR1** as the target (justified by data availability and pathway relevance)
- **Chose safety-first weighting** (chronic disease → prioritize QED/SA over raw activity)
- **Requested approval on 3 decisions** before computing anything (plan-before-compute)
- **Ran the full pipeline**: ChEMBL data → QSAR → VAE generation → ADMET gates → docking → retrosynthesis → scoring

### What the agent found wrong with its own work
Without being prompted, the agent identified 4 problems:

1. **Docking scores weren't competitive** (-7.58 vs co-crystal -10.23)
2. **Missed a dangerous motif** (hydrazone N-N bond, inappropriate for chronic dosing)
3. **Insufficient docking coverage** (only 5 of 3000 candidates docked)
4. **Root cause**: VAE KL collapse → molecules didn't resemble kinase inhibitors

The agent's diagnosis was confirmed by scaffold coverage analysis: <1% kinase hinge-binder motifs in the candidate pool vs 12-23% in training data.

After remediation (N-N deny list, Top 300 docking, scaffold analysis), best Vina score improved from -7.58 to -9.591.

See `case-studies/ipf-cycle1/` for the full report.

## Practical Lessons

**For agent trainers:**
- Don't teach what the agent already knows. Start with assessment.
- Skills encode decision boundaries, not API documentation.
- Challenge the agent to find its own mistakes. Hand-fed corrections don't stick.
- Keep skills <500 lines. Split references into separate files.

**For AI + drug discovery:**
- Self-diagnosis > raw capability. An agent that flags mediocre candidates as mediocre is more useful than one that calls them excellent.
- The generation problem remains hard. No amount of post-hoc filtering fixes a generator exploring the wrong chemical space.

## Repository Structure
```
├── README.md                            # This file
├── OpenClaw-Agent-Training-Guide.md     # Full methodology (Sections 1-12.9)
├── PLAYBOOK.md                          # Agent's operational playbook
├── skill-review-v3.md                   # Agent's self-review (3rd iteration)
├── skills/
│   ├── chem-literature/SKILL.md
│   ├── chem-qsar/SKILL.md
│   ├── chem-molgen/SKILL.md
│   ├── chem-retrosynthesis/SKILL.md
│   ├── chem-gnn/SKILL.md
│   ├── chem-admet/SKILL.md
│   ├── chem-rxn-conditions/SKILL.md
│   ├── chem-llm/SKILL.md
│   ├── chem-docking/SKILL.md
│   ├── chem-3dgen/SKILL.md
│   ├── chem-mlff/SKILL.md
│   └── chem-experiment/SKILL.md
└── case-studies/
    └── ipf-cycle1/
        ├── cycle1_summary.md
        └── scaffold_coverage.md
```

## Technical Stack

- **Platform**: [OpenClaw](https://github.com/openclaw-ai/openclaw)
- **Model**: GPT-5.2 via OpenAI Codex OAuth
- **Infrastructure**: NixOS on WSL2, Docker with GPU passthrough
- **Chemistry**: RDKit, AutoDock Vina, MACE-OFF, PyTorch, e3nn, ASE

## License

MIT License. Chemical data from ChEMBL is under CC BY-SA 3.0.

## Citation
```bibtex
@misc{chemicalexpert2026,
  title={Teaching AI Agents to Discover Drugs: A Systematic Training Methodology},
  author={Heng Gao},
  year={2026},
  url={https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training}
}
```
