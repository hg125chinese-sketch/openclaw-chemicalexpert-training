---
name: chem-cycle-learning
description: Learn reusable lessons from multiple chemistry campaign cycles by aggregating metrics, spotting recurring bottlenecks, comparing methods, and extracting validated lessons for future decision-making. Use after several DMTA cycles, case-study updates, or cross-cycle reviews.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "📚", "requires": { "bins": ["python3"], "python": ["pandas"] } } }
---

# chem-cycle-learning — turn history into strategy

Running more cycles is not the same thing as learning from them.

This skill teaches CE to look backward across many cycles and extract:
- what actually improved
- what stayed broken
- which filters matter most
- which tools earned trust
- which workflow changes should become policy

Core rule:
- **do not just archive history; compress it into reusable lessons**

This skill is retrospective.
- **chem-cycle-learning** = learn from past cycles
- **chem-reasoning** = interpret current evidence and propose next actions

Use cycle-learning to build the memory that reasoning can later reuse.

---

## Main outputs

### 1) `cycle_history.json`
A normalized per-cycle summary.

Recommended fields:

```json
{
  "cycle": 7,
  "generation_method": "DiffSBDD",
  "n_generated": 98,
  "rdkit_valid_rate": 0.949,
  "pb_valid_rate": 0.656,
  "safety_reject_rate": 0.16,
  "best_vina": -10.27,
  "hinge_rate": 0.25,
  "dft_pass_rate": 1.0,
  "best_score_final": 9.136
}
```

### 2) `lessons.json` or equivalent lesson table
A list of normalized lesson objects.

Schema:

```json
{
  "lesson_id": "L001",
  "source_cycles": [5, 6, 7],
  "category": "generation",
  "lesson": "DiffSBDD-era generation is cleaner than earlier VAE-era generation in this project.",
  "evidence": "Safety reject rates and PB pass rates improved materially after the switch.",
  "confidence": 0.82,
  "actionable": "Prefer DiffSBDD-style generation for future ALK5 cycles unless a target-specific reason suggests otherwise.",
  "status": "validated"
}
```

---

## Data aggregation workflow

### Step 1 — collect quantitative metrics
From each cycle’s `reports/` directory, extract or estimate:
- `generation_method`
- `n_generated`
- `rdkit_valid_rate`
- `pb_valid_rate`
- `safety_reject_rate`
- `best_vina`
- `hinge_rate`
- `dft_pass_rate`
- `best_score_final`

Use exact numbers when available.
If a metric is missing for older cycles:
- set it to `null`
- do not invent values
- note the missingness explicitly

### Step 2 — collect qualitative lessons
Read the relevant `case-studies/` summaries and project conclusion docs.
Extract:
- method changes
- failure analyses
- workflow upgrades
- explicit decisions that became policy

Examples:
- “raw DiffSBDD 3D is not QC-ready”
- “RDKit re-embed + MACE prescreen stabilized QE handoff”
- “Boltz-2 is a weak orthogonal signal on ALK5”

### Step 3 — normalize into `cycle_history.json`
Each cycle should become one record.
Keep the schema stable across cycles.

### Step 4 — derive lessons
Use the history table + case-study notes to generate lesson objects.
Each lesson should be:
- compact
- evidence-backed
- actionable
- tagged with confidence + status

---

## What to look for automatically

### A) Generation quality trend
Question:
- which generation regime produces the cleanest molecules?

Signals:
- RDKit-valid rate
- PoseBusters-valid rate
- safety reject rate

Typical example:
- compare VAE-era cycles vs DiffSBDD-era cycles

Output style:
- “DiffSBDD-era cycles produce cleaner molecules than earlier VAE-based cycles in this project.”

---

### B) Bottleneck identification
Question:
- which step removes the most candidates?

Signals:
- PB drop
- safety drop
- hinge filter
- multi-seed drop
- QE drop

Important:
- identify both the **largest numerical filter** and the **most decision-critical filter**

These are not always the same.

Example:
- PoseBusters may remove many candidates early
- multi-seed may remove fewer in absolute count but determine final lead quality

---

### C) Tool reliability
Question:
- which tools correlate best with downstream success?

Signals:
- Vina vs QE PASS
- GNINA vs PLIF recovery
- Boltz vs final panel outcome
- MACE prescreen vs DFT feasibility

Do not overclaim correlation from tiny sample sizes.
Use language like:
- “suggestive”
- “consistent with”
- “not yet enough to validate”

---

### D) Failure mode induction
Question:
- what kinds of candidates fail repeatedly?

Look for repeated patterns:
- OPT_FAIL chemotypes
- hinge instability patterns
- molecules that dock well but fail PLIF recovery
- raw 3D geometry failures

Goal:
- convert repeated pain into either:
  - a filter,
  - a policy,
  - or a hypothesis to test.

---

### E) Method-improvement tracking
Question:
- which skill upgrades changed measurable outcomes?

Examples:
- PoseBusters introduced as post-gen gate
- multi-seed docking added
- evidence schema added
- three-layer safety added
- entity resolver added
- RDKit re-embed + MACE prescreen adopted

For each, ask:
1. what changed?
2. what metric moved afterward?
3. how confident are we that the change mattered?

---

## Lesson object rules

Required fields:
- `lesson_id`
- `source_cycles`
- `category`
- `lesson`
- `evidence`
- `confidence`
- `actionable`
- `status`

Allowed categories:
- `generation`
- `filtering`
- `validation`
- `collaboration`
- `method`

Allowed status values:
- `validated`
- `preliminary`
- `hypothesis`

### Confidence guidance
- `0.2–0.4`: weak clue
- `0.4–0.7`: plausible pattern
- `0.7–0.9`: repeated and well-supported
- `>0.9`: rare; reserve for robust repeated observations

---

## Minimum lessons to extract from 7-cycle IPF history

At minimum, produce lessons covering these topics.

### L1 — DiffSBDD vs VAE generation quality
Theme:
- DiffSBDD-era generation appears cleaner than earlier VAE-era generation

Evidence to look for:
- PB-valid rates
- safety reject rates
- post-switch improvements in downstream survivorship

Likely category:
- `generation`

---

### L2 — PoseBusters as a post-gen gate
Theme:
- PoseBusters removes a large fraction of technically valid but poor-quality generated molecules

Evidence to look for:
- RDKit-valid vs PB-valid drop
- examples where PB prevented poor structures from propagating downstream

Likely category:
- `filtering`

---

### L3 — Multi-seed robustness is strict but valuable
Theme:
- multi-seed robustness acts as a hard narrowing stage for plausible binders

Evidence to look for:
- candidate count collapse in the hinge-positive panel (for example, `7 → 1` type filtering)
- robustness separating finalists better than raw docking alone

Likely category:
- `validation`

---

### L4 — Boltz-2 calibration behavior on ALK5
Theme:
- Boltz-2 is informative on ALK5, but should be treated as a weak orthogonal signal rather than a primary gate

Evidence to look for:
- calibration active mean
- disagreement cases with Vina / GNINA / PLIF
- standout signal such as Cycle 7 `mol_0021`

Likely category:
- `method`

---

### L5 — RDKit re-embed is critical for QE success
Theme:
- raw generator geometry is unsafe for QC handoff; RDKit re-embed + MACE prescreen improved DFT success materially

Evidence to look for:
- raw DiffSBDD strain failures
- later DFT PASS-rate stabilization after workflow change

Likely category:
- `validation`

---

### L6 — B3LYP-D3(BJ) vs B97-D efficiency / behavior
Theme:
- compare QE method choices as more cycles accumulate

Evidence to look for:
- walltime
- convergence behavior
- any systematic changes in returned properties or stability of optimization outcomes

Important:
- this lesson may begin as `preliminary` if sample size is still tiny

Likely category:
- `method`

---

### L7 — Three-layer safety changed the post-gen funnel
Theme:
- SMARTS + ADMET-AI + FAERS adds richer safety context without breaking backward compatibility

Evidence to look for:
- reject-rate changes
- survivor quality changes
- whether the softer layers added useful review flags instead of noisy hard rejects

Likely category:
- `filtering`

---

## How to write good lessons

Good lesson:
- specific
- tied to real cycles
- supported by both numbers and workflow context
- ends with a concrete future recommendation

Bad lesson:
- “Cycle 7 was good.”
- “Vina worked well.”
- “We should keep improving.”

Better:
- “Across Cycles 5–7, the combination of DiffSBDD generation and downstream structural gates produced cleaner candidates than earlier VAE-based cycles; future ALK5 work should keep DiffSBDD as default unless target-specific benchmarking shows otherwise.”

---

## Retrospective procedure

### Step 1 — build the history table
For each cycle, make one row.

### Step 2 — compare eras, not just neighbors
Examples:
- VAE era vs DiffSBDD era
- pre-PoseBusters vs post-PoseBusters
- pre-RDKit-reembed vs post-RDKit-reembed

### Step 3 — identify the top bottleneck per era
Different eras can fail for different reasons.

### Step 4 — extract lessons only where evidence exists
If the data is thin, emit a `preliminary` lesson or a `hypothesis`, not a `validated` one.

### Step 5 — connect lesson to future behavior
Every lesson should change something:
- tool choice
- gate policy
- ranking policy
- diagnostic habit
- experiment design

---

## Relationship to chem-reasoning

Keep the split clean:

### `chem-cycle-learning`
- retrospective
- aggregates many cycles
- extracts reusable lessons
- updates institutional memory

### `chem-reasoning`
- prospective / situational
- explains current evidence
- handles anomalies, conflicts, and live decisions
- may consume lessons from cycle-learning as prior knowledge

Short version:
- **cycle-learning writes the playbook**
- **reasoning uses the playbook**

---

## Minimal output set

When this skill runs well, it should usually produce:
1. `cycle_history.json`
2. `lessons.json` or lesson table
3. a short retrospective summary in markdown

Even if the user only asked for analysis, structure the results so later skills can reuse them.

---

## Bottom line

A multi-cycle program should get smarter over time.

This skill exists to make sure CE can say:
- what changed
- what worked
- what failed repeatedly
- what became policy
- and what future cycles should do differently

That is how past experiments become future advantage.
