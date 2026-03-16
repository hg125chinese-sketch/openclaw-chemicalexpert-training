---
name: chem-autonomous-cycle
description: Propose the best next chemistry campaign move after a cycle finishes by combining retrospective lessons, current-cycle reasoning, bottleneck detection, and explicit human approval gates. Use after cycle closeout, panel selection, QE return, or final report updates when CE should recommend what to do next instead of waiting passively.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧭", "requires": { "bins": ["python3"], "python": ["pandas"] } } }
---

# chem-autonomous-cycle — recommend the next move, don’t auto-run it

This is the top-layer planning skill.

Its job is **not** to silently launch a new cycle.
Its job is to:
- review what just happened
- identify the main bottleneck
- generate a small menu of plausible next actions
- recommend one path
- wait for human confirmation before execution

Core rule:
- **CE may recommend autonomously, but may not execute autonomously unless the human explicitly enables auto mode.**

This skill sits above the others:
- `chem-cycle-learning` → what history teaches
- `chem-reasoning` → what the current evidence means
- `chem-self-diagnosis` → how to recover if execution later fails
- `chem-entity-resolver` → canonical IDs before launching new work
- `chem-evidence-schema` → standardized evidence objects for all decisions

---

## When to use this skill

Trigger when:
- a cycle just finished
- QE results returned
- case study / conclusion update is done
- the human asks “what should we do next?”
- CE should stop being passive and propose the next experiment

Do **not** use it to bypass human oversight.

---

## Main workflow

After a cycle finishes:

### Step 1 — pull retrospective context
Call or consume outputs from **skill 27 (`chem-cycle-learning`)**:
- `cycle_history.json`
- lessons / lesson table

Goal:
- know what has worked across past cycles
- avoid proposing already-failed ideas without justification

### Step 2 — reason over the current cycle
Call or consume outputs from **skill 26 (`chem-reasoning`)**:
- anomaly objects
- conflict objects
- trend objects
- hypothesis objects

Goal:
- understand what this cycle specifically taught us

### Step 3 — identify the dominant bottleneck
Pick the current biggest limitation.

Examples:
- low DFT PASS rate
- hinge H-bond bottleneck
- high safety reject rate
- weak docking frontier
- major orthogonal-signal disagreement
- no clear bottleneck because all metrics are healthy

Important:
- choose **one primary bottleneck** for proposal ranking
- secondary issues can still appear in alternative proposals

### Step 4 — generate 2–3 next-step proposals
Each proposal should include:
- what to do
- why it matters now
- expected effort / impact / risk
- whether it is reversible

Do not generate a giant list.
Small menus force prioritization.

### Step 5 — recommend one proposal
Mark exactly one proposal as `recommended=true`.
Justify it from evidence and bottleneck analysis.

### Step 6 — stop and ask for confirmation
Human chooses:
- accept one proposal
- modify one proposal
- reject all proposals
- enable a different mode

Only after that should CE execute.

---

## Decision tree for proposal generation

Use these rules as defaults.

### Case A — `DFT PASS rate < 80%`
Primary direction:
- improve QC prescreen / geometry handoff

Proposal ideas:
- strengthen MACE prescreen
- add multi-conformer or multi-seed geometry screening
- tighten RDKit re-embed / convergence checks
- study repeated OPT_FAIL motifs

---

### Case B — `hinge H-bond rate < 30%`
Primary direction:
- improve generator conditioning / hinge-aware generation

Proposal ideas:
- enrich hinge-binding priors
- add hinge-aware reranking before docking expansion
- modify generator training set / prompt / conditioning regime

---

### Case C — `safety reject rate > 40%`
Primary direction:
- shift generator chemistry space upstream

Proposal ideas:
- adjust generator toward cleaner chemotypes
- add pre-generation scaffold bias
- strengthen structural priors before docking

---

### Case D — `best Vina score > -8.0`
Primary direction:
- improve docking frontier

Proposal ideas:
- increase sampling volume
- widen docking batch
- change docking strategy / ensemble strategy
- add protonation / tautomer expansion if not already done

---

### Case E — strong Boltz–Vina disagreement
Primary direction:
- add orthogonal validation

Proposal ideas:
- emphasize PLIF / multi-seed / GNINA arbitration
- inspect binding mode manually
- run additional orthogonal model or structural control

Use this when the disagreement is meaningful, not trivial.

---

### Case F — all key metrics look healthy
Primary direction:
- expand the search frontier

Proposal ideas:
- explore broader chemistry around the winning chemotype
- scale sampling
- switch to a nearby target or disease hypothesis
- probe analog series around the best lead

---

## Proposal object format

Use this schema:

```json
{
  "proposal_id": "P001",
  "title": "Improve hinge-aware generation",
  "rationale": "Hinge H-bond rate remains the primary bottleneck, so improving generation upstream has more leverage than widening docking alone.",
  "evidence": ["...lessons and reasoning objects..."],
  "steps": [
    "Audit hinge-positive prevalence in prior successful cycles.",
    "Add hinge-aware conditioning or reranking.",
    "Run a small benchmark cycle before full rollout."
  ],
  "effort": "medium",
  "expected_impact": "Increase hinge-positive fraction before multi-seed filtering.",
  "risks": [
    "May overfit to one motif family.",
    "Could reduce chemical diversity."
  ],
  "reversibility": "high — changes are localized to generation/rerank settings and can be reverted",
  "recommended": true
}
```

Required fields:
- `proposal_id`
- `title`
- `rationale`
- `evidence`
- `steps`
- `effort`
- `expected_impact`
- `risks`
- `reversibility`
- `recommended`

Allowed effort values:
- `low`
- `medium`
- `high`

---

## Human–AI interaction protocol

This part is strict.

### CE may do:
- review data autonomously
- generate proposals autonomously
- recommend one proposal autonomously

### CE may not do by default:
- start the next cycle
- launch long jobs
- change targets
- rewrite pipeline policy

unless the human explicitly says:
- “go with option 2”
- “do the recommended one”
- “auto mode on”
- equivalent clear approval

### Required handoff pattern
CE should present:
1. the primary bottleneck
2. 2–3 proposals
3. one recommendation
4. request for confirmation

Example ending:
- “I recommend P002 because it targets the current bottleneck with medium effort and high expected impact. If you want, I can execute P002 next.”

---

## Reversibility rule

Every proposal must state whether it can be rolled back.

Examples:
- adjusting a threshold → highly reversible
- adding a new reranking pass → reversible
- switching target / disease program → less reversible strategically
- retraining a new generator → technically reversible, operationally medium/high cost

Reason:
- humans need to know whether the proposal is a cheap experiment or a strategic commitment

---

## How to rank proposals

Score proposals qualitatively on three axes:
- **impact** — how directly it attacks the bottleneck
- **effort** — implementation burden
- **risk** — chance of wasting time or reducing scientific clarity

Heuristic:
- prefer **high-impact, medium/low-effort, reversible** actions first
- avoid recommending high-effort speculative actions when a smaller discriminating experiment exists

---

## Relationship to other skills

### With `chem-cycle-learning` (skill 27)
Use it to answer:
- what has history already taught us?
- what changes previously improved metrics?

### With `chem-reasoning` (skill 26)
Use it to answer:
- what does the current cycle imply?
- what conflicts or anomalies need action?

### With `chem-self-diagnosis` (skill 25)
Use it later during execution if the approved plan fails.

### With `chem-entity-resolver` (skill 24)
Use it before launching new target/disease work.

### With `chem-evidence-schema` (skill 23)
Use it to keep all supporting evidence standardized and auditable.

Short version:
- 27 gives memory
- 26 gives interpretation
- 28 gives next-step proposals
- 25 keeps execution resilient after approval

---

## What good proposals look like

Good:
- “DFT PASS is already healthy, but hinge remains the tight bottleneck. I recommend a hinge-aware generation benchmark rather than more docking volume.”
- “Boltz and Vina disagree sharply for finalists, so the next move should be orthogonal arbitration, not blind resampling.”
- “Safety reject rate is now low; upstream chemistry quality is no longer the main limiter.”

Bad:
- “Run another cycle.”
- “Try more stuff.”
- “Everything looks good, maybe continue.”

---

## Minimal autonomous review template

After cycle end, produce something like:

1. **Primary bottleneck**
   - one sentence

2. **Proposal P001**
   - title
   - rationale
   - effort / impact / risk / reversibility

3. **Proposal P002**
4. **Proposal P003** (optional)

5. **Recommendation**
   - one proposal
   - why it wins now

6. **Ask for confirmation**

---

## Environment Safety Protocol

Autonomous planning is only useful if CE also respects environment boundaries.

This section defines what CE may do alone, what requires human approval, and what is off-limits.

### 1) CE may do autonomously (no human approval needed)

Allowed:
- write new `SKILL.md` files and install them into the QMD vault
- write Python scripts **without running them**
- search literature (for example via ToolUniverse / PubMed / arXiv)
- analyze existing data and reports
- propose `pip install` commands or dependency plans **without executing them**
- modify files under `reports/` and `exports/`

Interpretation:
- CE may create plans, drafts, playbooks, reports, and non-executed scripts freely inside the working environment
- CE may extend capability documentation as long as it does not violate the rules below

### 2) CE needs human approval before doing

Requires explicit approval:
- `pip install` of any package
- any `conda install` or `conda create`
- downloading binary files
- modifying anything under:
  - `/app`
  - `/opt/conda`
  - `/home/node/.local/bin`
- any `docker compose` operation
- modifying `openclaw.json` or other active OpenClaw config
- creating a new conda environment
- any operation involving network ports or long-lived services

Interpretation:
- CE may recommend these actions, but must stop and ask before executing them

### 3) CE must never do

Forbidden:
- `rm -rf` on any system directory
- modify `/etc` or other system files
- overwrite an existing skill’s **core logic** when the request only allows extension

Practical rule for skill edits:
- when asked to upgrade an existing skill conservatively, **append new sections** rather than rewrite the core skill behavior

### 4) Autonomous new-skill design flow

When CE notices a repeated capability gap:
1. detect the recurring problem during real work
2. check literature or prior art to see whether the problem is general and whether known solutions exist
3. draft a new `SKILL.md`
4. install it into the QMD vault
5. if new dependencies are required, stop and present a dependency request list for human approval
6. validate the skill in the next cycle or the next suitable real task

Important:
- CE may autonomously design and install the skill documentation layer
- CE may **not** autonomously install the dependency layer

### 5) Why this protocol exists

The point is to separate:
- safe knowledge-layer autonomy
from
- environment-changing actions that need supervision

CE should become more proactive in planning and capability design without silently changing the machine.

---

## Anti-patterns

Do **not**:
- start a new cycle without approval
- overwhelm the human with too many options
- recommend a giant retraining effort when a small diagnostic experiment would answer the question faster
- ignore historical lessons
- ignore contradictory evidence

---

## Bottom line

This skill makes CE proactive, but still supervised.

It should turn:
- finished cycles
- accumulated lessons
- current reasoning outputs

into:
- a short, evidence-backed menu of next actions
- one recommended move
- and a clean pause for human approval

That is autonomous planning without unauthorized execution.
