# Teaching AI Agents to Discover Drugs: A Systematic Training Methodology

How I trained an LLM agent to autonomously run drug discovery workflows — and what it taught itself along the way.

## TL;DR

I trained an AI agent ("ChemicalExpert") with **21 specialized chemistry skills** using a systematic 7-step methodology. The agent learned to autonomously plan and execute Design-Make-Test-Analyze (DMTA) cycles for drug discovery. Over **6 iterative cycles** on the same target, the agent:

- Diagnosed its own output quality problems (VAE KL collapse, missing hinge binders)
- Ran 6 controlled experiments to fix molecular generation, documenting every failure
- Improved Top5 candidate quality from **20% to 80%** passing the critical hinge H-bond gate
- Discovered that SELFIES GRU VAE decoders structurally ignore conditioning signals — a finding validated across 6 independent experiments
- Upgraded to pocket-conditioned 3D diffusion (DiffSBDD), achieving **100% DFT pass rate** across two consecutive cycles
- Built a full multi-signal pipeline with PoseBusters geometry QC, GNINA CNN rescoring, multi-seed pose robustness, Boltz-2 affinity prediction, and conflict-aware panel selection

## The Problem

Large language models know a lot about chemistry. But knowing facts and *doing* chemistry are different things. An LLM can explain what QSAR is, but can it train a random forest on ChEMBL data, validate with scaffold splits, and tell you when the model is unreliable?

I wanted an agent that could:

- Run computational chemistry workflows end-to-end (data → models → molecules → evaluation)
- Know when to use each tool and when not to
- Catch its own mistakes before I catch them
- Work on new targets without being retrained

The platform is [OpenClaw](https://github.com/openclaw/openclaw), an open-source agent framework. The underlying model is GPT-5.2 via OpenAI's Codex API. The agent runs in a Docker container (GPU optional), with access to conda environments, RDKit, PyTorch, and standard scientific Python.

## The 7-Step Training Methodology

Through trial and error training two agents (a debugging agent and ChemicalExpert), I converged on a 7-step process:

**Step 1: Pre-Assessment** — Before writing any skill, ask the agent diagnostic questions across four dimensions: current knowledge, methodology, tool capabilities, and learning preferences. Let it tell you what it knows — don't assume.

**Step 2: Deep Probing** — Design pointed questions targeting suspected gaps. The goal is to expose the difference between "knows about" and "can do."

**Step 3: Gap Analysis** — Compare pre-assessment against probing results. Only fill gaps — don't reteach what the agent already knows.

**Step 4: Skill Design** — Each skill is a structured document (~200-700 lines) with decision boundaries, runnable code examples, philosophy sections, and failure modes. Key insight: Skills should be structured around "when to do what" rather than "how to do everything."

**Step 5: Guided Practice** — Give the agent a concrete task that requires the new skill. Watch it work. Intervene only when necessary.

**Step 6: Behavioral Correction** — When the agent makes mistakes, challenge it to find the error itself rather than giving the answer directly.

**Step 7: Self-Review** — Ask the agent to write its own skill review. Compare against your observations. Discrepancies reveal shallow understanding.

## The 21 Skills

| # | Skill | Core Capability |
|---:|---|---|
| 1 | QSAR Modeling | Scaffold splits, RF baseline, descriptor selection |
| 2 | Graph Neural Networks | GCN for molecular properties, scaffold generalization vs RF |
| 3 | ADMET Prediction | Rules-first (Lipinski/Veber), QED gates, chronic disease flags |
| 4 | Retrosynthetic Analysis | Template matching, forward validation, route scoring |
| 5 | Reaction Conditions | Yield prediction, reagent selection, total yield calculation |
| 6 | Chemical LLM | When to use LLM reasoning vs computation |
| 7 | Molecular Generation | SELFIES VAE, KL collapse detection, novelty metrics |
| 8 | 3D Generation | E(3)-equivariant diffusion, point cloud to molecule |
| 9 | Molecular Docking | Receptor prep, virtual screening, batch robustness, multi-seed pose verification, GNINA CNN rescoring |
| 10 | ML Force Fields | MACE-OFF geometry optimization, conformer ranking, QC prescreen with calibration |
| 11 | Experiment Planning | DMTA cycle orchestration, multi-objective scoring |
| 12 | Literature Review | Search strategies, source verification, evidence tables |
| 13 | Kinase SAR | Hinge binder motif analysis, coverage KPI, structure-activity |
| 14 | Reactivity Safety | Reactive group flagging, project denylists, chronic disease hard rejects |
| 15 | Protonation & Tautomers | pH-dependent multi-state enumeration, confidence scoring, best-state docking |
| 16 | Docking Interactions | ProLIF interaction fingerprints, hinge H-bond validation, PLIF recovery quantification |
| 17 | Scaffold-Conditioned Generation | VAE diagnostics, anti-collapse strategies, constrained decoding (Strategy E) |
| 18 | Pocket-Conditioned Diffusion | DiffSBDD 3D ligand generation, GPU inference, substructure inpainting |
| 19 | Binding Affinity Prediction | Boltz-2 affinity (isolated env), IC50 + binder probability, target calibration |
| 20 | Panel Selection | Multi-signal conflict-aware selection, 2x2 disagreement grid, auditable rationale |
| 21 | Structure QC Lite | PoseBusters geometry/pose plausibility, post-gen + post-dock checkpoints |

Each skill document is in `skills/*/SKILL.md`.

### Skill Design Principles

- **Philosophy over procedure.** "Always compare MACE to MMFF baseline" beats 50 lines of API docs.
- **Decision trees over checklists.** "If X, do A; if Y, do B" matches real drug discovery decisions.
- **Failure modes are non-negotiable.** Every skill includes what failure looks like.
- **Gates, not guidelines.** Hard thresholds (SA ≤ 6, QED > 0.3, hinge H-bond = required) that must be passed.
- **One sentence per skill.** If you can't state the core rule in one sentence, the skill is too broad.

## The Migration Test: IPF/ALK5 (6 Cycles)

After validating on DRD2, I gave the agent an open-ended task with no guidance:

> "Find anti-fibrosis (IPF) drug candidates. Pirfenidone and nintedanib have limited efficacy."

### Cycle 1: Autonomous Pipeline + Self-Diagnosis

**What the agent did autonomously:**
- Selected ALK5/TGFBR1 as the target (justified by data availability and pathway relevance)
- Chose safety-first weighting (chronic disease → prioritize QED/SA over raw activity)
- Requested approval on 3 decisions before computing anything (plan-before-compute)
- Ran the full pipeline: ChEMBL data → QSAR → VAE generation → ADMET gates → docking → retrosynthesis → scoring

**What the agent found wrong with its own work (unprompted):**
1. Docking scores weren't competitive (-7.58 vs co-crystal -10.23)
2. Missed a dangerous motif (hydrazone N-N bond, inappropriate for chronic dosing)
3. Insufficient docking coverage (only 5 of 3000 candidates docked)
4. Root cause: VAE KL collapse → molecules didn't resemble kinase inhibitors

After remediation, best Vina score improved from -7.58 to -9.591.

**Interaction analysis (skill 16) revealed:** 4/5 Top5 candidates lacked the fundamental hinge H-bond required for Type I kinase inhibition. Pure Vina ranking completely obscured this failure.

### Cycle 2: Fixing the Generator (6 Controlled Experiments)

The agent systematically attempted to fix the generation problem:

| Step | Strategy | Result | Coverage Ratio |
|------|----------|--------|---------------|
| 1 | Cyclical annealing (fix KL collapse) | Latent fixed (AU 0→32), but wrong chemical space | 0.021 |
| 2 | Training set enrichment | Already at target (47.7%), no change | 0.038 |
| 3 | Fragment-conditioned VAE (prefix) | Passed — but trivial solution (string concatenation) | 1.379 |
| 4 | Prefix=False validation | Confirmed: conditioning was fake | 0.002 |
| 5 | True latent conditioning (concat z, frag_embed) | Model ignores frag_embed | 0.054 |
| 6 | Auxiliary fragment classifier loss | Class imbalance bug found, fixed, still insufficient | 0.095 |
| **D** | **Rejection sampling (pragmatic fallback)** | **Passed all gates** | **1.573** |

**Key finding:** SELFIES GRU VAE decoders structurally ignore external conditioning signals. The autoregressive decoder learns to rely on its own hidden state history, bypassing any conditioning input (concat, auxiliary loss, cross-attention — all failed).

**Cycle 2 result:** 4/5 Top5 with hinge H-bond (vs 1/5 in Cycle 1).

### Cycle 3: Logit Bias Decoding

Built on the Cycle 2 finding, tested inference-time interventions:

- **Logit bias (+2.0 on fragment tokens):** Coverage ratio 1.393, efficiency 2.84% (2x better than rejection sampling)
- Cross-attention conditioning also attempted — failed (coverage 0.049), confirming the architectural limitation

**Cycle 3 result:** 4/5 Top5 with hinge H-bond, best Vina -8.927.

### Cycle 4: Full CE↔QE Collaboration

The complete multi-agent loop with QuantumExpert (QE) for DFT verification:

- Generated 1000 molecules (logit bias) → safety screen (466 survivors) → docked 50 → Top5 (3/5 hinge H-bond)
- QC pre-processing (RDKit + protomer + MACE prescreen) → 2 molecules handed to QE
- QE batch DFT (B97-D/def2-SVP) → 1 PASS, 1 OPT_FAIL (50% fail rate)
- Retrosynthesis completed for all PASS molecules

### Cycle 5: Pocket-Conditioned Diffusion (DiffSBDD)

Replaced the SELFIES VAE with DiffSBDD, a pocket-conditioned 3D diffusion model:

- **Generation:** 98 molecules from ALK5 pocket, 89 valid (90.8%), no logit bias needed
- **Safety:** 28% reject rate (vs VAE's 53%) — DiffSBDD produces cleaner molecules
- **Docking:** Best Vina -10.270, surpassing the co-crystal ligand (-10.23)
- **Key discovery:** DiffSBDD's 3D coordinates are optimized for docking pose, not molecular stability (MACE strain 450-630 kcal/mol). Solution: use DiffSBDD for molecular design, RDKit re-embed for QC geometry.
- **DFT:** 3/3 PASS (100% vs Cycle 3-4's 50%) — best candidate Vina -10.01, score_final 10.404

### Cycle 6: Full 21-Skill Pipeline

First cycle using all new tools (PoseBusters, multi-seed docking, GNINA rescoring, Boltz-2 affinity, panel selection):

- **Generation:** 91 molecules, 85 RDKit-valid (93.4%)
- **PoseBusters QC (new):** 56/91 PB-valid (61.5%) — caught 35% more geometry issues than RDKit alone
- **Safety:** 43 survivors (23% reject — lowest ever)
- **Docking:** Best Vina -10.04
- **Multi-seed robustness (new):** 7 hinge-positive molecules → only 1 survived 3-seed pose convergence + hinge consistency check
- **GNINA rescoring (new):** CNN-based scoring as orthogonal ranking signal
- **Boltz-2 affinity (new):** Weak signal on ALK5 (calibration: active mean binder_prob 0.27 vs decoy 0.10)
- **Panel selection (new):** Conflict-aware 2x2 grid, selected 3 candidates with relaxed hinge gate (≥0.67)
- **DFT:** 2/2 PASS (100%) — best candidate mol_0064, score_final **10.432** (new record)

### Progress Across Cycles

| Metric | Cycle 1 | Cycle 2 | Cycle 3 | Cycle 4 | Cycle 5 | Cycle 6 |
|--------|---------|---------|---------|---------|---------|---------|
| Generation method | Unconditional VAE | Rejection sampling | Logit bias | Logit bias | DiffSBDD | DiffSBDD |
| PB-valid rate | — | — | — | — | 55.1% | **61.5%** |
| Safety reject rate | 29% | 41% | 54% | 53% | 28% | **23%** |
| Best Vina score | -9.591 | -9.158 | -8.927 | -9.997 | **-10.270** | -10.04 |
| DFT PASS rate | — | — | 50% | 50% | 100% | **100%** |
| Best score_final | — | — | — | — | 10.404 | **10.432** |

**CE↔QE collaboration statistics:** 11 molecules submitted for DFT across 4 cycles. 8 PASS, 3 OPT_FAIL (27% overall fail rate). DiffSBDD era (Cycles 5-6): **5/5 = 100% PASS**.

## Architecture Lessons Learned

**1. VAE conditioning failure.** Six experiments across two cycles demonstrated that SELFIES GRU VAE decoders structurally ignore external conditioning signals. Inference-time intervention (logit bias, rejection sampling) is the validated solution. Documented as Strategy E in skill 17.

**2. Pocket-conditioned diffusion changes everything.** DiffSBDD produces better molecules (lower safety reject, better Vina scores, higher DFT pass rate) without needing any conditioning hacks. But its 3D coordinates are not suitable for QC — use the SMILES, not the geometry.

**3. Interaction analysis is mandatory.** Vina scores alone mask critical binding mode failures. In Cycle 1, 80% of top candidates lacked hinge H-bonds despite acceptable scores.

**4. Multi-signal validation catches silent failures.** Cycle 6 demonstrated that 7 molecules passing single-run hinge H-bond checks collapsed to just 1 when multi-seed robustness was applied. PoseBusters caught 35% more geometry issues than RDKit sanitization alone.

**5. Boltz-2 requires per-target calibration.** On ALK5, Boltz-2 binder probability shows weak but directional signal (active mean 0.27 vs decoy 0.10), while absolute IC50 predictions are unreliable (2 orders of magnitude off). Use as tiebreaker, not primary signal.

## Practical Lessons

**For agent trainers:**
- Don't teach what the agent already knows. Start with assessment.
- Skills encode decision boundaries, not API documentation.
- Challenge the agent to find its own mistakes. Hand-fed corrections don't stick.
- Keep skills < 700 lines. Split references into separate files.
- Every experiment gets committed — failures are as valuable as successes.

**For AI + drug discovery:**
- Self-diagnosis > raw capability. An agent that flags mediocre candidates as mediocre is more useful than one that calls them excellent.
- The generation problem remains hard. No amount of post-hoc filtering fixes a generator exploring the wrong chemical space.
- Pocket-conditioned generation (DiffSBDD) is a step change over string-based VAEs for structure-based drug design.
- Multi-agent collaboration works. CE→QE handoff with typed JSON contracts and explicit gates closed the loop from generation to quantum-chemical verification.
- Multi-signal validation is worth the compute cost. Single-metric ranking produces false positives that only surface in expensive downstream experiments.

## Repository Structure

```
├── README.md                            # This file
├── OpenClaw-Agent-Training-Guide.md     # Full methodology (Sections 1-12.9)
├── PLAYBOOK.md                          # Agent's operational playbook (21 skills)
├── SECURITY.md                          # Security & privacy guidelines
├── skills/
│   ├── chem-qsar/SKILL.md              # 1: QSAR modeling
│   ├── chem-gnn/SKILL.md               # 2: Graph neural networks
│   ├── chem-admet/SKILL.md             # 3: ADMET prediction
│   ├── chem-retrosynthesis/SKILL.md    # 4: Retrosynthetic analysis
│   ├── chem-rxn-conditions/SKILL.md    # 5: Reaction conditions
│   ├── chem-llm/SKILL.md              # 6: Chemical LLM
│   ├── chem-molgen/SKILL.md           # 7: Molecular generation
│   ├── chem-3dgen/SKILL.md            # 8: 3D generation
│   ├── chem-docking/SKILL.md          # 9: Molecular docking
│   ├── chem-mlff/SKILL.md             # 10: ML force fields
│   ├── chem-experiment/SKILL.md       # 11: Experiment planning
│   ├── chem-literature/SKILL.md       # 12: Literature review
│   ├── chem-kinase-sar/SKILL.md       # 13: Kinase SAR
│   ├── chem-reactivity-safety/SKILL.md # 14: Reactivity safety
│   ├── chem-protonation-tautomer/SKILL.md # 15: Protonation & tautomers
│   ├── chem-docking-interactions/SKILL.md # 16: Docking interactions
│   ├── chem-scaffold-conditioned-gen/SKILL.md # 17: Conditioned generation
│   ├── chem-pocket-diffusion/SKILL.md # 18: Pocket-conditioned diffusion
│   ├── chem-affinity-prediction/SKILL.md # 19: Binding affinity prediction
│   ├── chem-panel-selection/SKILL.md  # 20: Panel selection
│   └── chem-structure-qc-lite/SKILL.md # 21: Structure QC lite
├── case-studies/
│   ├── ipf-cycle1/                    # Cycle 1: autonomous pipeline + self-diagnosis
│   ├── ipf-cycle2/                    # Cycle 2: 6 conditioning experiments
│   ├── ipf-cycle3/                    # Cycle 3: logit bias decoding
│   ├── ipf-cycle4/                    # Cycle 4: full CE↔QE collaboration
│   ├── ipf-cycle5/                    # Cycle 5: DiffSBDD + 100% DFT pass
│   ├── ipf-cycle6/                    # Cycle 6: full 21-skill pipeline
│   └── ipf-project-conclusion.md      # Final results across all cycles
└── .gitignore
```

## Technical Stack

- **Platform:** [OpenClaw](https://github.com/openclaw/openclaw)
- **Model:** GPT-5.2 via OpenAI Codex API
- **Infrastructure:** Linux + Docker (GPU optional; configuration-dependent)
- **Chemistry:** RDKit, AutoDock Vina, GNINA 1.3, ProLIF, MACE-OFF, DiffSBDD, PoseBusters, PyTorch, e3nn, ASE
- **Affinity prediction:** Boltz-2 (isolated conda environment)
- **Quantum Chemistry (via QuantumExpert):** PySCF (B97-D/def2-SVP)

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
