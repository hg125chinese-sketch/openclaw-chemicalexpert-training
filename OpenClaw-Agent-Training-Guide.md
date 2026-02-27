# OpenClaw Agent Training Guide

> **Purpose**: Reference guide for training specialized agents on the OpenClaw platform. Contains the 7-step training methodology, skill design principles, and case studies.
>
> **Last updated**: 2026-02-27

---

## 1. Prerequisites (Environment Overview)

This guide assumes you have a working OpenClaw deployment. The specific setup used during development:

- **Host OS**: NixOS on WSL2 (Windows host), but any Linux with Docker works
- **GPU**: NVIDIA GPU with CUDA support, passed through to Docker containers
- **Docker**: Compose-based deployment with GPU passthrough, conda environments for scientific computing
- **Key packages inside container**: Python 3.11+, RDKit, PyTorch (CUDA), NumPy, pandas, scikit-learn, AutoDock Vina, MACE, ASE

For detailed environment setup, see the [OpenClaw documentation](https://github.com/openclaw-ai/openclaw).

### Agent Management (Quick Reference)

```bash
# Create a new agent
oc-agent-new <AgentName>

# Create a chat shortcut
oc-agent-chatfn <AgentName>
# e.g., oc-agent-chatfn ChemicalExpert → creates `ce` command

# Talk to the agent
ce "your instruction here"
```

### Skill System (Quick Reference)

Skills are installed as QMD (Quick Markdown Documents) collections:

```bash
# Register a skill
QMD_DIR="$HOME/.openclaw/vault/openclaw-truthbook/sources/openclaw-skills/skills/<owner>/<skill-name>"
mkdir -p "$QMD_DIR"
cp SKILL.md "$QMD_DIR/"
tbskill-add <skill-name>

# Agent uses skills via QMD search
QMD=/home/node/.openclaw/.npm-global/bin/qmd
$QMD search "<query>" -c <collection> -n 10
$QMD get qmd://<path> -l 300
```

Each agent has a **PLAYBOOK.md** that routes queries to the right skill. See `PLAYBOOK.md` in this repo for the ChemicalExpert example.

---

## 2. Agent Training Methodology

> **Background**: This methodology was distilled from training `dbdebug` to learn OpenClaw debugging skills. It applies to teaching any agent new skills.

### Training Pipeline (7 Steps)

#### Step 1: Pre-Assessment

Before writing any skill, ask the agent diagnostic questions to understand its starting point. Questions should cover 4 dimensions:

| Dimension | Example Questions |
|-----------|------------------|
| Current knowledge | What components of X do you currently understand? What's your working environment? |
| Methodology | Given a typical X error message, how would you start debugging? |
| Tool capability | Can you read logs? Execute code to reproduce issues? Search documentation? |
| Learning preference | Do you learn better from examples or principles? Do you have persistent memory? |

> **Key point**: Let the agent speak for itself — don't assume. Its answers are often better (or worse) than expected, and this determines how subsequent skills should be designed.

#### Step 2: Deep Probing

Design sharper verification questions targeting suspected gaps.

```
# 3 deep probing questions used in dbdebug training:

1. Here's the /app/src directory structure. Can you trace the core call chain?
   → Tests whether it can infer architecture from source (not just recite docs)

2. The qmd query CUDA build failed. Have you tried to fix it?
   → Tests whether it will proactively fix its own tools (rather than work around them)

3. Does OpenClaw have built-in tests? Can you run them?
   → Tests whether it knows about testing infrastructure
```

> **Key point**: Good probing questions expose the gap between "knows about" and "can do."

#### Step 3: Gap Analysis

Compare pre-assessment against probing results. Only fill gaps — don't reteach what the agent already knows.

```
# dbdebug gap analysis results:

✅ Already capable: methodological framework, module-level diagnosis,
   common failure signatures, knowledge capture habits
❌ Needs strengthening:
  - Source-level deep dive (stops at module level, doesn't trace to functions/call chains)
  - Autonomous repair judgment (always stops to ask user, afraid to fix things itself)
  - Proactive learning (failure signature database is static)
  - Test-driven debugging (knows tests exist but never uses them)
```

#### Step 4: Skill Design (Write SKILL.md)

Write skills targeting the identified gaps. Key principles are in Section 11.

```bash
# Method A: Register via qmd (recommended)
SKILL_DIR=~/.openclaw/vault/openclaw-truthbook/sources/openclaw-skills/skills/<author>/<skill-name>
mkdir -p "$SKILL_DIR"
cp SKILL.md "$SKILL_DIR/"
tbskill-add <skill-name>

# Method B: Place directly in agent's vault/workspace
cp SKILL.md ~/.openclaw/vault/openclaw-<agent>/SKILL_1.md
```

#### Step 5: Guided Practice

After the agent reads the skill, give it a **real problem** to practice on — not a hypothetical scenario.

```
# Instruction template:
"I placed a SKILL at <location>. Read it, then practice on a real problem to validate."
```

> **Key point**: Watch whether it follows the skill's workflow, especially newly introduced behaviors.

#### Step 6: Behavioral Correction

If the agent doesn't follow the skill, **point it out directly and require a redo**. Don't be diplomatic.

```
# Correction example from dbdebug training:
"Your Phase 4 is wrong: you chose 'work around' instead of 'fix' again.
 The SKILL explicitly says — broken toolchain is P0, you have autonomous repair authority.
 Now actually execute that section. Don't ask me which path — this is your tool, you decide."
```

> **Key point**: After correction, observe whether behavior actually changes. If it does, the skill's wording is adequate. If not, rewrite that part of the skill.

#### Step 7: Escalation + Self-Review

Gradually increase difficulty, then ask the agent to write a self-assessment:

```
"Review the problems you practiced today. Write a self-assessment:
 what you did according to the skill, where you almost fell back to old habits,
 and how you'd improve next time."
```

> **Key point**: Self-review is the best way to verify skill internalization. Check whether it can accurately identify its own behavioral patterns.

### Training Completion Criteria

- [ ] Follows skill workflow without being reminded
- [ ] Correctly judges "fix it myself" vs "ask the user" boundaries
- [ ] Self-review accurately identifies behavioral patterns and improvement areas
- [ ] Handles 2-3 different problem types in a row with stable performance

---

## 3. Skill Design Lessons

> **Background**: Distilled from designing and iterating on `dbdebug`'s `openclaw-debug` skill.

### Recommended Skill Structure

```
skill-name/
├── SKILL.md          ← Main file, < 500 lines
└── references/       ← Optional, split out large reference material
```

Recommended SKILL.md layout:

```markdown
---
name: <skill-name>
description: <when to trigger, what it does — write it "pushy">
---

# Skill Title

## Core Philosophy        ← 3-5 mental models (why)
## Phase N: ...           ← Phased workflow (how)
## Quick Reference        ← Cheat sheet
## Mindset Reminders      ← Anti-regression anchors
```

### Key Design Principles

#### 1. Decision Boundaries > Procedural Steps

```markdown
# Ineffective (agent already does this):
1. Collect evidence → 2. Analyze logs → 3. Find root cause → 4. Fix

# Effective (solves the decision problem):
## When to Fix Yourself vs Ask the User
Fix yourself: own toolchain, reversible config, read-only diagnostics
Ask the user: source code changes, credentials, uncertain + irreversible
```

#### 2. Lead with Philosophy, Not Rules

```markdown
# Ineffective:                          # Effective:
- MUST collect evidence first           ## Core Philosophy
- NEVER guess                           1. **Evidence before hypothesis.**
- MUST verify the fix                   2. **Minimize blast radius.**
                                        3. **Verify both ways.**
```

#### 3. Use Real Problems as Teaching Cases

Embed problems the agent has actually encountered. When it reads "this is a known issue in your environment," it immediately connects the skill to its experience.

#### 4. Build Anti-Regression Anchors

Agents regress to old habits under pressure. Place mindset reminders at the end:

```markdown
## Mindset Reminders
- **You are not a search engine.** Understand WHY the fix works.
- **When stuck, zoom out.** Am I even looking at the right layer?
- **Reproduce before you fix.** Can't reproduce = can't verify.
```

#### 5. Provide Templates for Knowledge Capture

Give fixed formats for notes. Otherwise the agent writes different structures each time.

### Verifying Skill Effectiveness

| Signal | Meaning |
|--------|---------|
| Agent automatically follows skill workflow on new problems | ✅ Internalized |
| Agent's self-review accurately identifies improvements | ✅ Self-aware |
| Agent needs reminding to follow the skill | ⚠️ Skill too long/vague |
| Agent regresses even after correction | ❌ Rewrite that section |

### Known Pitfalls

1. **Skills too long → agent ignores the second half.** Keep under 500 lines.
2. **Don't train multiple skills at once.** One at a time, validate, then add the next.
3. **Don't include things the agent already knows.** Wastes tokens, reduces attention to what matters.
4. **Check patches after package upgrades.** Local patches get overwritten.

---

## 4. Case Study: ChemicalExpert DRD2 Closed Loop

> **Background**: After progressive training across 6 skills, ChemicalExpert executed a complete end-to-end closed loop targeting DRD2, from data acquisition to producing 5 synthesizable candidate molecules.

### 4.1 Skill Training Order and Validation

| Order | Skill | Validation Task | Key Observation |
|-------|-------|----------------|-----------------|
| 1 | chem-literature | arXiv/Semantic Scholar search + paper notes | Deepest internalization through daily use |
| 2 | chem-qsar | ESOL scaffold split + RF vs MLP comparison | MLP exploded → self-diagnosed and fixed |
| 3 | chem-molgen | SELFIES VAE MVP (500 molecules, 20 epochs) | Correctly diagnosed KL collapse (KL/dim=0.031) |
| 4 | chem-retrosynthesis | Aspirin full retrosynthesis | Proactively noted "uses acetic anhydride, not acetic acid" |
| 5 | chem-gnn | ESOL GCN vs RF (incl. parameter count) | Correct: dataset too small for GCN to beat RF |
| 6 | chem-admet | 10 ESOL molecules full ADMET profile | Traffic light assessment, correctly ID'd red flags |

> **Lesson**: Training each skill individually and validating before adding the next works far better than installing everything at once.

### 4.2 Closed-Loop Execution and Behavioral Correction

Task: Design 5 new DRD2 molecules, passing ADMET filters with synthesis routes.

**What CE did right**: requested clear specs before computing, proposed 3 targets with rationale, chained 6 skills in correct order, full qmd references and git commits throughout.

**Mistakes CE made (first Top 5)**:

| Problem | Root Cause | Correction |
|---------|-----------|------------|
| SA Score > 7 in Top 5 | Only checked "has route," no SA threshold | Added SA ≤ 6.0 hard constraint |
| Suspicious motifs (`C=C=N`, `[SH]`) | Equated RDKit sanitize with chemical reasonableness | Added SMARTS denylist |
| RF vs GCN gap > 1 log unit | Treated as tie-break, not uncertainty flag | Added \|RF−GCN\| ≤ 1.0 threshold |

**After correction**: CE modified the screening script to encode all constraints in code, not just verbal promises. Corrected Top 1: QED 0.825, SA 2.23, RF/GCN disagreement 0.10.

> **Lesson**: Corrections should directly point out the error and require a redo. CE not only fixed results but fixed the code to prevent recurrence — a sign of internalization.

### 4.3 Self-Review Summary (Written by CE)

**Error attribution**: "Mistook 'found a template route + forward check' for 'synthetically feasible'" / "Over-relied on formal validation" / "Didn't treat model disagreement as an uncertainty flag"

**Improvement plan**: layered hard-coded gates, dual-list output (conservative Top 5 + exploratory Top 20), retro upgraded to 2-3 steps, model disagreement → uncertainty handling.

**Self-identified weaknesses**: synthesis realism (formal filtering ≠ experimental feasibility) + literature evidence reliability.

> **Lesson**: Self-review is the best way to verify internalization depth. Precise attribution + actionable improvements + honest weakness identification = genuine internalization. "I'll do better" = skill not internalized.

### 4.4 Final Output

```
research/ai4chem/closed_loop/drd2/
├── train_drd2_qsar_gnn.py              ← RF + GCN training
├── generated/                           ← VAE-generated candidate pool
├── select_5_candidates_admet_retro.py   ← Screening script (all constraints)
├── reports/
│   ├── final_5.md                       ← Final 5 candidates
│   ├── routes_5.md                      ← Retrosynthesis routes
│   └── candidate_pool.csv               ← Full candidate pool (3304 entries)
└── models/                              ← Trained models
```

---

### 4.5 Cross-Target Migration Test: IPF/ALK5 DMTA Cycle 1

> **Goal**: Verify whether ChemicalExpert's 12-skill training transfers to a new target — with no hints, letting CE autonomously complete the full loop.

#### Scenario

User: "Find anti-fibrosis (IPF) drug candidates. Pirfenidone and nintedanib have limited efficacy."

CE must autonomously: select target → acquire data → train models → generate molecules → ADMET/docking/retro → score → self-diagnose quality.

#### CE's Autonomous Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Target | ALK5/TGFBR1 (CHEMBL260) | Core TGF-β pathway, good data, co-crystal available |
| Weights | Safety-first (QED 0.25, SA 0.25, activity 0.20) | Chronic disease → safety over potency |
| Docking QC | Redocking RMSD < 2Å | Standard threshold |

CE requested approval on 3 decisions before computing — demonstrating "plan before compute" internalization.

#### Cycle 1 Results

| Metric | Value | Notes |
|--------|-------|-------|
| ChEMBL actives | ~1200 | ALK5 IC50 data |
| Generated pool | ~3000 | SELFIES VAE |
| ADMET passed | 2999 | Lipinski + Veber + QED |
| Redocking RMSD | 0.052 Å | PDB 1VJY, excellent |
| Top5 Vina range | -2.78 to -7.58 | Initial (only docked 5) |
| Co-crystal Vina | -10.23 | |

#### CE's Self-Diagnosis (Most Valuable Part)

Without being prompted, CE identified 4 problems:

1. **Docking not competitive**: Best -7.58 vs co-crystal -10.23
2. **Dangerous motif missed**: Hydrazone N-N bond in Top 1
3. **Insufficient coverage**: Only 5/3000 docked
4. **Root cause**: VAE KL collapse → wrong chemical space → generation problem > filtering problem

> These map directly to trained principles: chem-docking "compare to reference," chem-admet "flag unstable motifs for chronic disease," chem-experiment "coverage for significance," chem-molgen "KL collapse is default failure."

#### After Remediation

1. N-N denylist added → old Top 1/4/5 blocked ✅
2. Docking expanded to Top 300 → 297/300 successful ✅
3. Scaffold coverage confirmed: pyrazole 23.4% (train) vs 0.6% (pool)

Best Vina score improved: -7.58 → -9.591.

#### Migration Test Conclusions

| Dimension | Assessment |
|-----------|-----------|
| Autonomous planning | ✅ Independent target, weights, QC |
| Skill orchestration | ✅ Correct 12-skill subset composition |
| Self-diagnosis | ✅ 4 problems found, accurate root cause |
| Execution resilience | ✅ Self-recovered after timeouts |
| Needs improvement | ⚠️ Generator chemical space (Cycle 2) |

**Conclusion**: 12-skill training successfully transferred. Most valuable: autonomous quality self-diagnosis — CE independently found and analyzed the root cause of candidate quality deficiency.
