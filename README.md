# Multi-Agent AI for Computational Science: A Systematic Training and Collaboration Methodology

**One person + OpenClaw + 7-step training = 3 cooperating AI agents for end-to-end drug discovery**

## Abstract

We trained three specialized AI agents — ChemicalExpert (28 skills), ProteinEngineer (7 skills), and QuantumExpert (4 skills) — using a systematic 7-step methodology on the OpenClaw platform. The agents share a vault-based skill system and collaborate through typed JSON handoff protocols. Validated on an IPF/ALK5 drug discovery campaign over seven DMTA cycles plus an analog exploration campaign, the system achieved 100% DFT pass rate across three consecutive cycles and the analog campaign (9/9 in the DiffSBDD era) using pocket-conditioned 3D diffusion (DiffSBDD) and a full multi-signal validation pipeline (PoseBusters geometry QC, GNINA CNN rescoring, multi-seed pose robustness, Boltz-2 affinity prediction, conflict-aware panel selection, ToolUniverse-powered target validation, and standardized evidence objects). All training artifacts are open-source.

---

## 1. Introduction

Large language models know chemistry, protein engineering, and quantum mechanics — but knowing and doing are different. An LLM can explain QSAR, but can it train a random forest on ChEMBL data with scaffold splits and flag when its own model is unreliable? Can it diagnose why its molecular generator produces candidates that dock well but lack the pharmacophoric interactions required for activity?

Single-agent approaches hit a depth ceiling: one agent cannot maintain expert-level skill across medicinal chemistry, protein design, and quantum chemistry simultaneously. A multi-agent system with specialized roles and structured collaboration protocols can.

**What we built:**

- **ChemicalExpert (CE)**: 28 skills across four layers — computational chemistry (1-21), scientific infrastructure (22-24), and cognitive capabilities (25-28) — covering molecular generation (VAE + diffusion), QSAR, docking (Vina + GNINA), interaction analysis (ProLIF + PLIF recovery), three-layer safety screening, geometry QC (PoseBusters), binding affinity prediction (Boltz-2), multi-signal panel selection with hinge-aware pre-ranking, synthesis planning, target validation (ToolUniverse), entity resolution, evidence standardization, self-diagnosis, scientific reasoning, cross-cycle learning, and supervised autonomous cycle planning
- **ProteinEngineer (PE)**: 7 skills for protein interface redesign, stability assessment, and developability screening
- **QuantumExpert (QE)**: 4 skills for DFT calculations, excited-state methods, xTB prescreening, and cross-agent workflow orchestration

All three agents share the same 7-step training methodology, the same vault-based skill system, and can exchange data through typed JSON contracts.

---

## 2. Platform and Architecture

### 2.1 OpenClaw

OpenClaw is an open-source agent framework that provides the runtime for all three agents. Each agent runs in a Docker container with access to domain-specific conda environments and tools.

### 2.2 System Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    Shared Vault                          │
│  (openclaw-truthbook)                                   │
│  ┌──────────────┐ ┌──────────────┐ ┌─────────────────┐  │
│  │ CE Skills    │ │ PE Skills    │ │ QE Skills       │  │
│  │ (28 QMD      │ │ (7 QMD       │ │ (4 QMD          │  │
│  │ collections) │ │ collections) │ │ collections)    │  │
│  └──────────────┘ └──────────────┘ └─────────────────┘  │
│  ┌──────────────┐ ┌──────────────┐                      │
│  │ Debug Skills │ │ ToolUniverse │                      │
│  │ (2 QMD)      │ │ (1996 tools) │                      │
│  └──────────────┘ └──────────────┘                      │
└────────────┬───────────────┬───────────────┬────────────┘
             │               │               │
     ┌───────▼────────┐ ┌───▼─────┐ ┌───────▼────────┐
     │ ChemicalExpert │ │ Protein │ │ QuantumExpert  │
     │   Workspace    │ │Engineer │ │   Workspace    │
     │                │ │Workspace│ │                │
     │ exports/ ──────┼─┼─────────┼─▶ imports/       │
     │ imports/ ◀─────┼─┼─────────┼── exports/       │
     └────────────────┘ └─────────┘ └────────────────┘
             │               │               │
     ┌───────▼────────┐ ┌───▼─────┐ ┌───────▼────────┐
     │  conda: chem   │ │  prot   │ │  chem (PySCF)  │
     │  RDKit, Vina   │ │ ESM-2   │ │  DFT, TD-DFT   │
     │  ProLIF, MACE  │ │ESMFold  │ │  CASSCF/NEVPT2 │
     │  DiffSBDD      │ │         │ │  xTB 6.7.1     │
     │  GNINA, PB     │ │         │ │                │
     │  ToolUniverse  │ │         │ │                │
     ├────────────────┤ └─────────┘ └────────────────┘
     │ conda: boltz   │
     │  Boltz-2       │
     │  (isolated)    │
     └────────────────┘
```

### 2.3 Shared Vault Architecture

The vault is the central knowledge store that makes multi-agent collaboration possible. It contains:

- **QMD collections**: searchable, versioned skill documents that agents query at runtime via `qmd search` and `qmd get` with line-number evidence
- **Skills**: modular capability definitions (SKILL.md files, 200-700 lines each) with decision boundaries, failure modes, and runnable examples
- **Cross-references**: any agent can cite specific lines from any skill document, enabling cross-domain auditability

Each agent reads from the shared vault but writes only to its own workspace. This prevents cross-contamination while enabling knowledge sharing.

### 2.4 Agent Dispatch

Each agent has its own workspace, identity (IDENTITY.md), behavioral principles (SOUL.md), and operational playbook (PLAYBOOK.md). Agents are invoked via command-line dispatch (`ce "..."`, `qu "..."`, `pr "..."`) and execute tasks within their workspace.

### 2.5 Cross-Agent Communication

Agents communicate through structured JSON files placed in `exports/` and `imports/` directories. The handoff protocol is defined by explicit schemas — the sending agent packages data with all required fields, and the receiving agent validates the schema before processing. This decoupled design means agents can be developed and tested independently.

---

## 3. The 7-Step Training Methodology

Through iterative refinement across all three agents, we converged on a 7-step training process that proved transferable across domains.

### Step 1: Pre-Assessment

Before writing any skill, ask the agent diagnostic questions across four dimensions: current knowledge, methodology, tool capabilities, and learning preferences. Let it reveal what it already knows.

### Step 2: Deep Probing

Design pointed questions targeting suspected gaps. The goal is to expose the difference between "knows about" and "can do." For example, an agent might explain KL collapse in VAEs but fail to diagnose it from training logs.

### Step 3: Gap Analysis

Compare pre-assessment against probing results. Only fill gaps — don't reteach what the agent already knows. This keeps skills lean and focused.

### Step 4: Skill Design

Each skill is a structured document (200-700 lines) with:

- **Decision boundaries**: "if X, do A; if Y, do B" — matching real domain decisions
- **Hard gates**: non-negotiable thresholds that must pass before proceeding
- **Failure modes**: what failure looks like and what to do about it
- **Runnable examples**: concrete code and commands, not abstract descriptions
- **Philosophy sections**: one-sentence core principles that guide behavior

### Step 5: Guided Practice

Give the agent a concrete task requiring the new skill. Watch it work. Intervene only when necessary.

### Step 6: Behavioral Correction

When the agent makes mistakes, challenge it to find the error itself rather than giving the answer. Hand-fed corrections don't stick; self-diagnosed corrections become permanent behavioral changes.

### Step 7: Self-Review

Ask the agent to write its own skill review. Compare against your observations. Discrepancies reveal shallow understanding that needs reinforcement.

### Cross-Domain Transferability

This methodology was applied identically across chemistry (CE), protein engineering (PE), and quantum chemistry (QE). The same structure — pre-assess, probe, fill gaps, design skills with gates, practice, correct, review — produced effective agents in all three domains.

---

## 4. Agent Profiles

### 4.1 ChemicalExpert (CE)

**Mission.** CE is a pragmatic "lab-mate" agent for digital chemistry: molecular representations, property modeling, molecule generation, safety filtering, docking/interaction analysis, reproducible experiment pipelines, and autonomous cycle planning.

**Skill set (28 skills in 4 layers):**

| Layer | Skills |
|-------|--------|
| **Computational chemistry (1-21)** | chem-qsar, chem-gnn, chem-admet, chem-retrosynthesis, chem-rxn-conditions, chem-llm, chem-molgen, chem-3dgen, chem-docking (+ multi-seed + GNINA), chem-mlff, chem-experiment, chem-literature, chem-kinase-sar, chem-reactivity-safety (3-layer), chem-protonation-tautomer, chem-docking-interactions (+ PLIF recovery), chem-scaffold-conditioned-gen, chem-pocket-diffusion, chem-affinity-prediction (Boltz-2), chem-panel-selection (+ hinge-aware pre-ranking), chem-structure-qc-lite (PoseBusters) |
| **Scientific infrastructure (22-24)** | chem-target-validation (ToolUniverse 10-phase), chem-evidence-schema (T1-T4 grading + conflict detection), chem-entity-resolver (multi-DB canonicalization) |
| **Cognitive capabilities (25-28)** | chem-self-diagnosis (failure classification + auto-recovery), chem-reasoning (anomaly/conflict/trend/hypothesis), chem-cycle-learning (cross-cycle lessons), chem-autonomous-cycle (supervised next-step proposals) |

**Design principles:**

1. Gating over guessing: mandatory KPIs (e.g., hinge H-bond) are hard gates, not soft scores
2. Separation of concerns: generation, safety, docking, interpretation, and QC are explicit pipeline stages
3. Failure must be informative: every failure produces a traceable artifact
4. Cross-agent contracts: handoffs use strict schemas with minimal assumptions
5. Multi-signal validation: no single metric is trusted alone; panel selection resolves disagreements
6. Cognitive loop: after each cycle, CE can autonomously propose next steps but requires human approval before execution

**Capability boundary.** CE is not authoritative for final synthetic feasibility without external validation, DFT-level electronic properties (delegated to QE), or claims of true binding affinity from docking scores.

### 4.2 ProteinEngineer (PE)

**Mission.** PE is an autonomous agent for protein engineering workflows, emphasizing interface redesign campaigns that must survive multi-tool scrutiny.

**Skill set (7 skills):**

| Skill | Core Function |
|-------|---------------|
| prot-seqdesign | Constrained sequence design (LigandMPNN + ESM-IF cross-validation) |
| prot-structure | Structure evaluation and fold-back validation (ESMFold) with WT sanity check |
| prot-stability | ESM-2 log-likelihood ratio scoring with interface-mutation mode |
| prot-interface | Interface extraction, contact typing, hotspot proxy, burial-aware analysis |
| prot-msa | Pseudo-conservation and consensus guidance |
| prot-developability | Sequence-based red-flag screening with WT-relative gating |
| prot-dbtl | Design-Build-Test-Learn cycle orchestration |

**Key training discoveries:**

- **1YCR (short peptide binder)**: Exposed metric pathologies — pLDDT scaling artifacts, clash counting without bonded exclusions. Led to peptide-aware fold-back rules.
- **1BRS (standard protein interface)**: Revealed tool disagreement as diagnostic signal — fold-back PASS + interface PASS but stability FAIL. Led to interface-mutation mode with burial-aware strictness.
- **1Z92 (cytokine interface, full DBTL cycle)**: Uncovered system-level blockers — ESMFold cannot reliably fold disulfide-rich proteins (IL-2, pLDDT≈45). Led to WT sanity check mechanism (HARD_GATE/DIAGNOSTIC/SKIP modes) and WT-relative developability screening.

**Capability boundary.** PE uses evolutionary and structural proxies, not physics-based energy calculations.

### 4.3 QuantumExpert (QE)

**Mission.** QE provides quantum-chemistry rigor inside the multi-agent workflow, turning candidate molecules into verifiable electronic-structure outputs with explicit QC gates.

**Skill set (4 skills):**

| Skill | Core Function |
|-------|---------------|
| qchem-dft | DFT calculations, SCF troubleshooting, geometry optimization, frequency analysis |
| qchem-excited-state | TD-DFT, multi-reference escalation (SA-CASSCF + NEVPT2), active space validation |
| qchem-workflow | Orchestration layer with batch screening, checkpointing, and failure classification |
| qchem-prescreen-xtb | Fast xTB 6.7.1 prescreening for geometry sanity before DFT |

**Key training discoveries:**

- Gate-first execution: if SCF doesn't converge, don't report energies; if geometry isn't converged, don't proceed to frequencies
- Method labeling as a contract: every number is tied to functional/basis/software version
- Batch robustness: per-molecule independence with checkpointing prevents one failure from stopping the batch
- Dispersion correction vigilance: QE initially ran B97-D without verifying dispersion was active — a silent failure mode. Permanently encoded in PLAYBOOK: "don't lower standards, change route" (switch to B97-D or ωB97M-V rather than silently dropping D3)

**Capability boundary.** QE works within PySCF's capabilities. Methods requiring unavailable corrections are flagged, not silently downgraded. GPU4PySCF is verified usable but marked "DF compatibility pending" due to CuPy 14 incompatibility.

---

## 5. Cross-Agent Collaboration Protocol

### 5.1 Handoff Schema

CE → QE handoff uses a strict JSON contract:

```json
{
  "source_agent": "ChemicalExpert",
  "task": "qc_screening",
  "project": "IPF_ALK5_cycle7",
  "molecules": [
    {
      "mol_id": "mol_0021",
      "smiles": "...",
      "charge": 0,
      "spin": 0,
      "geometry_source": "rdkit_embed",
      "geometry_file": null
    }
  ],
  "requested_properties": ["E_total", "HOMO", "LUMO", "gap", "dipole"]
}
```

QE returns results with `qc_flag` per molecule (PASS / OPT_FAIL / SCF_FAIL), enabling CE to make gate decisions without interpreting raw QC data.

### 5.2 Gate Contracts

Each handoff includes pre-processing gates on the sending side:

- **CE pre-processing before QE handoff**: PoseBusters QC → safety screen → Vina docking → hinge-aware pre-ranking → multi-seed robustness → GNINA rescoring → ProLIF + PLIF recovery → panel selection → RDKit 3D embedding + MMFF convergence → MACE geometry prescreen → JSON packaging
- **QE gates during DFT**: SCF convergence → geometry optimization convergence → frequency analysis (n_imag = 0) → property extraction

### 5.3 Failure Handling

The protocol is designed so failures at any stage produce actionable diagnostics rather than silent drops:

- OPT_FAIL molecules are flagged but retained in the audit trail
- Repeated OPT_FAIL on the same SMILES across cycles triggers a structural blacklist
- MACE prescreen rejections are logged with strain values for calibration
- PoseBusters failures are logged with specific failure reasons

### 5.4 Potential Cross-Agent Interfaces (Not Yet Tested)

- **PE↔CE (protein-ligand co-design)**: CE proposes ligands; PE redesigns the protein interface to improve complementarity.
- **PE↔QE (mechanism-driven enzyme engineering)**: QE provides energetic insights on catalytic steps; PE translates those into sequence changes.

---

## 6. End-to-End Case Study: IPF/ALK5

### 6.1 Problem Statement

Idiopathic pulmonary fibrosis (IPF) has limited treatment options. We tasked the multi-agent system with finding novel ALK5/TGFBR1 inhibitor candidates, starting from an open-ended prompt with no target guidance.

### 6.2 Target Validation (Skill 22)

Before any DMTA cycle, CE performed systematic target validation using the ToolUniverse 10-phase framework, querying Open Targets, ChEMBL, UniProt, Pharos, PDB, and PubMed. Result: **76/100 (T1 Strong GO)** — strong biology/chemistry validation but safety-constrained (TGF-β pathway systemic risk).

### 6.3 Cycle 1: Autonomous Pipeline + Self-Diagnosis

CE autonomously selected ALK5, chose safety-first weighting (chronic disease), and ran the full pipeline. Self-diagnosis (unprompted): 4/5 Top5 candidates lacked hinge H-bond. Root cause: VAE KL collapse.

**Result:** Top5 hinge H-bond rate: 1/5 (20%)

### 6.4 Cycle 2: Fixing the Generator (6 Controlled Experiments)

CE systematically tested 6 conditioning strategies, concluding that SELFIES GRU VAE decoders structurally ignore external conditioning signals. Adopted rejection sampling as pragmatic fallback.

**Result:** Top5 hinge H-bond rate: 4/5 (80%)

### 6.5 Cycle 3: Strategy E (Logit Bias Decoding)

Logit bias (+2.0 on fragment tokens) achieved 2x better efficiency than rejection sampling.

**Result:** Top5 hinge H-bond rate: 4/5 (80%)

### 6.6 Cycle 4: Full CE↔QE Collaboration

Complete multi-agent loop. QE: 1 PASS, 1 OPT_FAIL (50% fail rate).

### 6.7 Cycle 5: Pocket-Conditioned Diffusion (DiffSBDD)

Replaced SELFIES VAE with DiffSBDD. Key discovery: DiffSBDD's 3D coordinates are for molecular design, not QC geometry (MACE strain 450-630 kcal/mol). Best Vina -10.270, surpassing the co-crystal ligand. DFT: 3/3 PASS (100%).

### 6.8 Cycle 6: Full Multi-Signal Validation Pipeline

First cycle using PoseBusters, multi-seed docking, GNINA, Boltz-2, and panel selection. PoseBusters caught 35% more geometry issues than RDKit alone. Multi-seed reduced 7 hinge-positive candidates to 1 robust survivor. DFT: 2/2 PASS (100%). Best candidate mol_0064, score_final **10.432**.

### 6.9 Cycle 7: Full 28-Skill Pipeline with Cognitive Capabilities

First cycle using entity resolution, three-layer safety (SMARTS + ADMET-AI + openFDA/FAERS), evidence schema, and ToolUniverse integration:

- **Entity resolution (new):** TGFBR1/ENSG00000106799/P36897 + EFO_0000768 canonicalized before any computation
- **PoseBusters QC:** 61/93 PB-valid (65.6%) — new record
- **Three-layer safety:** 51 survivors (16% reject — all-time lowest)
- **Boltz-2 affinity:** mol_0021 binder_prob = **0.698** — highest ever on ALK5 (calibration active mean: 0.268)
- **DFT:** 1/1 PASS (100%, B3LYP-D3(BJ)/def2-SVP) — mol_0021

### 6.10 Analog Exploration: mol_0021 Neighborhood Validation

A 10-member analog panel around mol_0021 validated that the lead is not a one-off hit:

- 7/10 analogs retained hinge H-bond (70%)
- Top 3 (A5_01, A3_02, A3_01) all passed strict 3-seed robustness + QE DFT (3/3 PASS)
- **A3_02**: strongest Vina–Boltz consensus (Vina -10.00, Boltz binder_prob 0.704 — exceeding parent)
- **A5_01**: strongest docking + GNINA (Vina -10.12, CNNscore 0.703)
- **A3_01**: strongest QE profile (gap 4.79 eV, dipole 1.62 D, score_final 10.635 — project-wide #1)
- All analogs carry N–N safety caution (scaffold-level, not compound-specific)

The mol_0021 chemotype neighborhood is real and survives the full CE↔QE validation stack.

### 6.11 Cognitive Capabilities in Practice

After Cycle 7, CE used its cognitive skill stack for the first time:

- **Skill 27 (cycle-learning):** Aggregated 7 cycles into `cycle_history.json` + 7 validated lessons with confidence grading
- **Skill 28 (autonomous-cycle):** Identified hinge-compatible generation as the primary bottleneck, proposed 3 next-step options, recommended hinge-aware generation benchmark
- **Skill 20 update (from P001 backtest):** Hinge-aware pre-ranking improved Top20 hinge rate from 25% to 65% on Cycle 7 data — adopted as default

### 6.12 Cumulative Results

**Progress across cycles:**

| Metric | Cycle 1 | Cycle 2 | Cycle 3 | Cycle 4 | Cycle 5 | Cycle 6 | Cycle 7 |
|--------|---------|---------|---------|---------|---------|---------|---------|
| Generation method | Unconditional VAE | Rejection sampling | Logit bias | Logit bias | DiffSBDD | DiffSBDD | DiffSBDD |
| PB-valid rate | — | — | — | — | 55.1% | 61.5% | **65.6%** |
| Safety reject rate | 29% | 41% | 54% | 53% | 28% | 23% | **16%** |
| Best Vina score | -9.591 | -9.158 | -8.927 | -9.997 | **-10.270** | -10.04 | -10.27 |
| DFT PASS rate | — | — | 50% | 50% | 100% | 100% | **100%** |
| Best score_final | — | — | — | — | 10.404 | **10.432** | 9.136 |
| Boltz-2 best binder_prob | — | — | — | — | — | 0.12 | **0.698** |

**CE↔QE collaboration statistics:** 15 molecules submitted for DFT across 5 cycles + analog campaign. 12 PASS, 3 OPT_FAIL (20% overall fail rate). DiffSBDD era (Cycles 5-7 + analogs): **9/9 = 100% PASS**.

---

## 7. Key Findings

### 7.1 SELFIES GRU VAE Conditioning Failure

Six controlled experiments demonstrated that SELFIES GRU VAE decoders structurally ignore external conditioning signals. Inference-time intervention (logit bias, rejection sampling) is the validated solution.

### 7.2 Pocket-Conditioned Diffusion Changes Everything

DiffSBDD produced better molecules across every metric and achieved 100% DFT pass rate across three consecutive cycles plus an analog campaign (9/9). However, its 3D coordinates are optimized for docking pose, not molecular stability. Use DiffSBDD for molecular design, RDKit for geometry preparation.

### 7.3 Multi-Signal Validation Catches Silent Failures

PoseBusters caught 35% more geometry issues than RDKit alone; multi-seed docking reduced 7 hinge-positive candidates to just 1 robust survivor; GNINA and Boltz-2 provided orthogonal ranking signals whose disagreements are informative.

### 7.4 Interaction Analysis is Mandatory

Pure docking scores mask critical binding mode failures. In Cycle 1, 80% of top candidates lacked the fundamental hinge H-bond despite acceptable Vina scores.

### 7.5 Boltz-2 Requires Per-Target Calibration

On ALK5, Boltz-2 binder probability shows weak but directional signal (active mean 0.27 vs decoy 0.10). But mol_0021's 0.698 proved to be a genuine signal — the analog campaign confirmed that its neighborhood retains strong Boltz-2 values (A3_02: 0.704).

### 7.6 Three-Layer Safety Outperforms Rules Alone

SMARTS structural alerts + ADMET-AI ML predictions + real-world evidence (openFDA/FAERS) reduced safety reject rate from 23% to 16% while maintaining protective gates.

### 7.7 Entity Resolution Prevents Silent Evidence Corruption

Standardizing target/molecule/disease IDs before computation prevents errors from different tools using different identifiers for the same entity.

### 7.8 Self-Diagnosis Outweighs Raw Capability

An agent that flags its own mediocre output is more useful than one that presents flawed results confidently. CE's unprompted identification of KL collapse, QE's detection of silent dispersion failures, PE's diagnosis of fold-back rejecting valid designs — in each case, self-diagnosis led to permanent behavioral improvement.

### 7.9 Cognitive Capabilities Close the Planning Loop

Skills 25-28 allow CE to move beyond executing instructions to proposing next steps. The hinge-aware pre-ranking update (25%→65% improvement) was proposed, backtested, and adopted through the autonomous cycle planning workflow — not by human command.

### 7.10 Analog Campaigns Validate Lead Neighborhoods

mol_0021 was initially a single robust hit. The analog campaign proved the neighborhood is real: 7/10 hinge retention, 3/3 DFT PASS, and A3_01 achieving the project-wide highest score_final (10.635). This transforms a single lead into a validated chemotype series.

---

## 8. Limitations and Future Work

### 8.0 Time and Compute Summary

| Item | Approximate Time |
|------|-----------------|
| CE training (28 skills + 7 cycles + analogs) | ~5 weeks of iterative sessions |
| PE training (7 skills + 3 practice runs) | ~1 week |
| QE training (4 skills + validation) | ~4 days |
| One CE↔QE collaboration round (2-3 molecules) | ~6-10 hours DFT wall time |
| Full Cycle 7 (generation → all gates → DFT → scoring) | ~20 hours total |
| Analog campaign (10 analogs → validation → DFT) | ~24 hours total |

All computation ran on a single machine (Docker on Linux/WSL2, RTX 4080 Laptop GPU).

### 8.1 Current Limitations

- **MACE prescreen**: Current strain proxy does not correlate with DFT OPT_FAIL. Needs richer signals or more calibration data.
- **DiffSBDD 3D coordinates**: Not suitable for direct QC input (strain 450-630 kcal/mol). RDKit re-embed is still required.
- **Boltz-2 on ALK5**: Weak signal; per-target calibration essential. Analog campaign showed it can identify genuine binders but should not be the sole gate.
- **N–N safety constraint**: The mol_0021 series carries a scaffold-level N–N bond safety flag. Further progression requires explicit medchem assessment of this tradeoff.
- **No experimental validation**: All results are computational.
- **GPU4PySCF**: Verified usable but density fitting compatibility pending (CuPy 14 issue).

### 8.2 Future Directions

- **Hinge-aware generation**: Pharmacophore-constrained diffusion or DiffSBDD substructure inpainting to directly enforce hinge motifs (P001 identified this as the primary remaining bottleneck)
- **PE integration**: Use PE to engineer ALK5 stability variants for assay development
- **Higher-level QC**: TZVP refinement or DLPNO-CCSD(T) for final candidate ranking
- **Three-agent loop**: CE generates ligands, PE optimizes the protein target, QE validates binding energetics
- **Expanded ToolUniverse integration**: More skills using the 1996-tool SDK for automated evidence gathering

---

## 9. How to Reproduce

To train a new agent using this methodology:

1. **Set up OpenClaw** with a Docker container and the relevant conda environment for your domain.
2. **Pre-assess** the base model by asking diagnostic questions across your domain's four dimensions.
3. **Design skills iteratively** — start with 2-3 core skills, validate on a practice task, then expand. Each skill should fit in one SKILL.md file (200-700 lines) with decision boundaries and hard gates.
4. **Use a real project as the training ground** — abstract exercises don't produce robust agents. Our agents became effective through IPF drug discovery (CE), protein interface redesign (PE), and DFT validation (QE).
5. **For multi-agent collaboration**, define the JSON handoff schema before the first data exchange. Include `mol_id`, all required fields, and a `flag` field for gate status.

Each of the three repositories below contains the complete skill set, training guide, and case studies needed to reproduce or extend the corresponding agent.

---

## 10. Repository Links

| Repository | Content | Skills |
|------------|---------|--------|
| [openclaw-chemicalexpert-training](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) | CE training, 28 skills, IPF case study (7 cycles + analog exploration) | 28 |
| [openclaw-proteinengineer-training](https://github.com/hg125chinese-sketch/openclaw-proteinengineer-training) | PE training, 7 skills, 3 protein engineering cases | 7 |
| [openclaw-quantumexpert-training](https://github.com/hg125chinese-sketch/openclaw-quantumexpert-training) | QE training, 4 skills, CE↔QE collaboration test | 4 |

---

## 11. Citation

```bibtex
@misc{multiagent-compchem-2026,
  title={Multi-Agent AI for Computational Science: A Systematic Training and Collaboration Methodology},
  author={Heng Gao},
  year={2026},
  url={https://github.com/hg125chinese-sketch}
}
```
