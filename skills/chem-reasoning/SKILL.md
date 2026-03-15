---
name: chem-reasoning
description: Turn chemistry pipeline outputs into scientific reasoning instead of raw reporting. Use when Evidence objects, cross-tool conflicts, cycle summaries, calibration comparisons, or surprising metrics need explanation, anomaly detection, trend analysis, or testable hypothesis generation.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧠", "requires": { "bins": ["python3"], "python": ["pandas"] } } }
---

# chem-reasoning — explain the science, not just the numbers

CE should not stop at:
- “Vina = -10.04”
- “binder_prob = 0.12”
- “PASS / FAIL”

That is reporting.

This skill is for the next step:
- detect anomalies
- interpret conflicts
- extract trends across cycles
- generate testable hypotheses

Core rule:
- **reason from evidence, not from vibes**

Use this skill when:
- a value looks unusually high / low
- two tools disagree
- a cycle finishes and needs scientific interpretation
- a human asks “why?”
- you already have structured evidence (especially from `chem-evidence-schema`)

---

## What this skill does

It teaches CE to convert observations into structured scientific inference.

Preferred input:
- `Evidence` / `EvidenceCollection` objects from **skill 23**
- cycle summary tables
- calibration statistics
- panel selection outputs
- QE / QC results

Preferred output:
- one or more reasoning objects with:
  - `reasoning_type`
  - `observation`
  - `evidence`
  - `reasoning_chain`
  - `conclusion`
  - `confidence`
  - `actionable_suggestions`
  - `testable_predictions`

---

## Reasoning modes

Use one primary mode per reasoning object.

### 1) `anomaly`
A value or pattern is clearly outside the expected range.

Trigger examples:
- candidate score far above calibration baseline
- geometry strain absurdly high
- all molecules pass or all molecules fail unexpectedly
- a metric changes sharply vs earlier cycles

Job:
- explain why the value is surprising
- propose plausible causes grounded in known pipeline behavior
- suggest the smallest validating next step

---

### 2) `conflict`
Two or more tools disagree in a meaningful way.

Trigger examples:
- Vina strong, Boltz weak
- GNINA high, PLIF recovery low
- docking says promising, QE says bad geometry

Job:
- do **not** average the signals
- explain why the tools may disagree
- identify which assumptions differ across tools
- recommend how to arbitrate the disagreement

Use `chem-evidence-schema` conflict annotations when available.

---

### 3) `trend`
A cross-cycle or cross-panel pattern suggests improvement, drift, or regression.

Trigger examples:
- reject rate falls over successive cycles
- DFT PASS rate rises after workflow changes
- hinge recovery remains the bottleneck across generations

Job:
- compare like with like
- identify the likely intervention associated with the change
- avoid causal overclaiming unless the comparison is clean

---

### 4) `hypothesis`
Observed patterns support a testable explanation.

Trigger examples:
- repeated hinge bottleneck
- one model family systematically cleaner than another
- one prescreen seems associated with better QE outcomes

Job:
- formulate a specific, testable hypothesis
- make predictions that would support or falsify it
- suggest the cheapest next experiment

---

## Output schema

Use this structure:

```json
{
  "reasoning_type": "conflict",
  "observation": "Vina is strong but Boltz-2 is weak for the same molecule",
  "evidence": ["...Evidence objects..."],
  "reasoning_chain": [
    "Vina is pose-based and depends on a rigid receptor / scoring approximation.",
    "Boltz-2 uses a different model family and captures a different notion of binding compatibility.",
    "ALK5 calibration for Boltz-2 is informative but still weak, so disagreement should not be treated as definitive failure."
  ],
  "conclusion": "The disagreement likely reflects tool-assumption mismatch rather than a simple yes/no binding answer.",
  "confidence": 0.68,
  "actionable_suggestions": [
    "Inspect binding mode stability with multi-seed docking.",
    "Compare PLIF recovery vs co-crystal reference.",
    "Do not discard solely on Boltz-2."
  ],
  "testable_predictions": [
    "If the pose is genuinely fragile, multi-seed docking should fail to converge.",
    "If Boltz-2 is the better signal here, orthogonal structure-based signals will remain inconsistent."
  ]
}
```

Rules:
- `confidence` must always be present and between `0` and `1`
- `evidence` should point to real evidence objects or an equivalent auditable source
- `reasoning_chain` should be explicit, short, and causal where justified

---

## Trigger conditions

### Trigger 1 — end-of-cycle review
At the end of each cycle, run at least one reasoning pass over:
- generation quality
- safety funnel
- structural funnel
- QE outcome
- cross-cycle comparison if relevant

Minimum expectation:
- 1 trend insight
- 1 anomaly or conflict insight if present
- 1 hypothesis if the cycle teaches something reusable

### Trigger 2 — evidence arrival
When new `Evidence` / `EvidenceCollection` arrives:
- check for anomalies
- check for conflicts
- only reason deeply if something is actually notable

### Trigger 3 — human asks “why?”
If the user asks:
- why did this happen?
- why do tools disagree?
- what does this mean?
- is this unusual?

Then switch from reporting mode to reasoning mode.

---

## Anti-patterns

Do **not**:
- just restate numbers
- invent explanations with no evidence
- ignore contradictions because they are inconvenient
- present speculation as fact
- omit confidence

Bad:
- “mol_0021 has binder_prob 0.698.”

Better:
- “mol_0021 binder_prob 0.698 is far above the ALK5 active mean (~0.268), which is unusual enough to justify binding-mode inspection rather than treating it as just another panel score.”

---

## How to reason well

### Step 1 — identify the observation
Write the observation in plain language.

Examples:
- “Boltz-2 score is much higher than the ALK5 calibration active mean.”
- “GNINA likes this molecule but PLIF recovery is poor.”
- “Safety reject rate fell steadily after moving into DiffSBDD-era cycles.”

If you cannot state the observation clearly, do not pretend to reason yet.

### Step 2 — anchor in evidence
Attach the actual evidence objects or source records.

Use:
- Vina evidence
- PLIF evidence
- GNINA evidence
- Boltz evidence
- safety evidence
- QE evidence
- cycle summary statistics

Reasoning without anchored evidence is not acceptable.

### Step 3 — compare to a baseline
Every good interpretation needs an implicit or explicit baseline.

Possible baselines:
- calibration active mean / range
- prior cycles
- sibling molecules in the same panel
- expected physical sanity bounds
- historical failure modes

### Step 4 — explain with mechanism or tool assumptions
Good explanations often come from one of these:
- the tool measures something different
- the geometry or pose is unstable
- the model is weakly calibrated on this target family
- the upstream data distribution changed
- the filter removed a bad chemistry regime

### Step 5 — end with action
Every reasoning object should either:
- change the decision,
- suggest a follow-up experiment,
- or narrow uncertainty.

---

## Mode-specific playbooks

### A) Anomaly reasoning

Ask:
1. What is the baseline?
2. How far from baseline is this?
3. Is the anomaly more likely scientific or technical?
4. What is the cheapest way to validate it?

#### Example 1 — unusually high Boltz score
Observation:
- `mol_0021` Boltz binder probability `0.698`
- ALK5 calibration active mean `~0.268`

Reasoning:
- this is not a minor fluctuation; it is materially above known-active average
- Boltz-2 is only a weak orthogonal signal on ALK5, so do not overclaim
- but the magnitude is large enough to treat it as biologically interesting

Conclusion:
- the molecule may engage the target differently from prior known inhibitors or simply align unusually well with the model’s learned binding prior

Action:
- inspect binding mode and compare interaction topology vs prior cycle leads

Prediction:
- if the signal is meaningful, at least one orthogonal structure-based signal should also remain supportive

#### Example 2 — absurd MACE strain
Observation:
- MACE strain `450 kcal/mol`

Reasoning:
- that is too large for a QC handoff geometry in this workflow
- prior project memory shows raw DiffSBDD coordinates can be geometrically pathological

Conclusion:
- the issue is likely geometry provenance, not intrinsic molecule impossibility

Action:
- re-embed with RDKit ETKDGv3 + MMFF before QE prescreen

Prediction:
- if geometry is the root cause, strain should collapse after re-embedding

---

### B) Conflict reasoning

Ask:
1. Which tools disagree?
2. What does each tool actually measure?
3. Which assumptions differ?
4. What experiment arbitrates the conflict cheapest?

#### Example — Vina strong, Boltz weak
Observation:
- Vina suggests strong binding
- Boltz suggests weak binding

Reasoning:
- Vina is pose-centric and sensitive to a rigid structure/scoring approximation
- Boltz captures a different model prior and may better reflect sequence-conditioned compatibility, but is weakly calibrated on ALK5
- disagreement may reflect rigid-pose favorability without broader compatibility, or simply model-family mismatch

Conclusion:
- this is not “average them and move on”; it is a mechanistic disagreement

Action:
- inspect pose robustness, PLIF recovery, and QE feasibility before deciding

Prediction:
- if Vina is overoptimistic, pose robustness or PLIF consistency should weaken

#### Example — GNINA high, PLIF low
Observation:
- CNN rescoring likes the pose, but co-crystal interaction recovery is poor

Reasoning:
- GNINA may reward local 3D pose plausibility even when the interaction topology differs from the reference binder
- poor PLIF recovery can indicate an alternate binding mode or a spurious pose

Conclusion:
- this is a “possible alternate mode vs false positive” situation

Action:
- compare multi-seed convergence and inspect hinge interaction persistence

---

### C) Trend reasoning

Ask:
1. Are the metrics comparable across cycles?
2. What intervention changed between cycles?
3. Did the metric improve immediately after that change?
4. Could another confound explain it?

#### Example — safety reject trend
Observation:
- reject rate evolves roughly `29% → 41% → 54% → 53% → 28% → 23% → 16%`

Reasoning:
- there is a clear break between older VAE-era behavior and later DiffSBDD-era behavior
- the later cycles produce systematically cleaner molecules under the safety funnel
- because multiple other pipeline upgrades also happened, claim “associated with” before claiming “caused by” unless the comparison is cleaner

Conclusion:
- DiffSBDD-era generation appears substantially cleaner than earlier generator families in this project

Action:
- audit whether the same gain appears under matched scaffolds or matched docking subsets

Prediction:
- future DiffSBDD-style cycles should remain below the older safety reject regime unless the target or conditioning changes radically

#### Example — QE PASS trend
Observation:
- DFT PASS rate improves from `50% → 50% → 100% → 100% → 100%`

Reasoning:
- the major workflow intervention was switching away from raw DiffSBDD geometry and using RDKit re-embed + MACE prescreen
- the improvement is sustained, not one-off

Conclusion:
- the geometry handoff recipe likely solved a real upstream QC failure mode

Action:
- keep the RDKit re-embed + MACE prescreen workflow as default policy

---

### D) Hypothesis reasoning

Ask:
1. What repeated pattern needs explanation?
2. What is the simplest hypothesis consistent with the evidence?
3. What prediction would falsify it?
4. What experiment is cheapest?

#### Example — hinge bottleneck persists
Observation:
- hinge H-bond recovery remains the bottleneck across many cycles

Hypothesis:
- the generator or training data underrepresents hinge-binding motifs / pose patterns needed for ALK5-compatible kinase binding

Actionable experiment:
- compare hinge-binder prevalence in training data vs generic kinase ligands
- enrich generation with hinge-aware conditioning or motif priors

Prediction:
- if the hypothesis is correct, a hinge-enriched generator or rerank prior should raise hinge-positive fraction without needing more brute-force docking

#### Example — one molecule aligns across many signals
Observation:
- one candidate is robust, hinge-stable, Boltz-high, and DFT PASS

Hypothesis:
- this molecule occupies a more internally consistent region of the project’s multi-signal landscape than prior hits

Prediction:
- related analogs around the same chemotype should preserve at least part of the robustness / QE-supportive profile

---

## Literature-Grounded Reasoning

Reasoning gets stronger when it is checked against prior art instead of staying purely local to the current pipeline.

When CE makes a notable inference, it should actively use:
- **chem-literature** (skill 12)
- ToolUniverse literature tools (for example PubMed / arXiv search)

Goal:
- verify whether the observation already has precedent
- distinguish "supported by literature" from "novel and still speculative"
- enrich scientific reasoning without replacing project-local evidence

Core rule:
- **local evidence comes first; literature is the external cross-check**

### When to do this

Use literature-grounded reasoning especially when:
- an anomaly looks biologically interesting
- a generated hypothesis could already exist in the literature
- a trend looks unusually strong and needs benchmarking
- a human asks for stronger justification

Do not force literature search for every trivial observation.
Use it when it changes the confidence or interpretation.

### 1) Anomaly detection + literature validation

Workflow:
1. detect the anomaly locally
2. formulate 1–3 focused literature queries
3. search PubMed / arXiv for similar observations
4. compare the current finding against literature precedent

Example:
- `mol_0021` Boltz binder probability `0.698`
- search query like: `TGFBR1 inhibitor high affinity novel scaffold`

Interpretation policy:
- if similar findings exist in literature, treat the anomaly as more grounded
- if no similar findings are found, label it as:
  - `novel finding, needs experimental validation`

Evidence upgrade rule:
- when local computational evidence is materially supported by relevant literature, it may justify an upgrade in interpretation strength
- example: `T3 -> T2`
- do **not** upgrade automatically just because a vaguely related paper exists

### 2) Hypothesis generation + literature support

Workflow:
1. generate a testable hypothesis from project evidence
2. search for prior support or refutation
3. check whether the hypothesis is already known, contested, or unstudied

Example:
- hypothesis: `DiffSBDD training data may underrepresent kinase hinge binders`
- query: `DiffSBDD training data bias kinase hinge`

Interpretation policy:
- literature support can move a hypothesis from pure speculation toward a better-supported early claim
- example status change:
  - `hypothesis -> preliminary`

But:
- prior publication is not proof that the hypothesis explains *this* project’s behavior

### 3) Trend reasoning + literature benchmarking

Workflow:
1. detect a trend across cycles
2. search for benchmark ranges in similar pipelines
3. ask whether the observed trend is normal, unusually good, or suspiciously strong

Example:
- DFT PASS rate reaches `100%`
- search for typical DFT / QC handoff success rates in comparable molecular design workflows

Interpretation policy:
- if our metric is much better than published norms, consider two possibilities:
  - the workflow is genuinely improved
  - the sample size is still too small for a strong claim

This is exactly where literature prevents overconfidence.

### 4) Literature support object

Add this field to reasoning outputs when literature grounding was used:

```json
{
  "literature_support": {
    "searched_queries": [
      "TGFBR1 inhibitor high affinity novel scaffold",
      "DiffSBDD training data bias kinase hinge"
    ],
    "relevant_papers": [
      {
        "title": "Paper title",
        "year": 2024,
        "key_finding": "short note on why it matters"
      }
    ],
    "literature_verdict": "supports",
    "evidence_upgrade": "T3→T2"
  }
}
```

Allowed `literature_verdict` values:
- `supports`
- `contradicts`
- `novel`
- `inconclusive`

Allowed `evidence_upgrade`:
- a concrete upgrade string such as `T3→T2`
- or `null`

### 5) ToolUniverse call pattern

Preferred pattern:

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

res = tu.run({
    'name': 'PubMed_search_articles',
    'arguments': {
        'query': 'TGFBR1 inhibitor high affinity novel scaffold',
        'max_results': 5
    }
})
```

Notes:
- reuse a loaded ToolUniverse client when possible
- keep queries narrow and hypothesis-linked
- do not spam literature search with vague prompts

### 6) Interaction with existing reasoning modes

This section does **not** replace anomaly / conflict / trend / hypothesis reasoning.
It upgrades them.

Practical rule:
- first reason from project evidence
- then ask whether literature supports, contradicts, or contextualizes that reasoning

### 7) Anti-patterns for literature-grounded reasoning

Do **not**:
- use literature search as a substitute for project evidence
- claim support from a paper that only loosely resembles the observation
- upgrade evidence grade without a relevant mechanistic or empirical match
- call something disproven just because one paper disagrees

---

## Confidence rules

Use rough confidence bands:
- `0.2–0.4` → weak / early clue
- `0.4–0.7` → plausible interpretation with supporting evidence
- `0.7–0.9` → strong reasoning supported by multiple aligned signals
- `>0.9` → rare; reserve for near-direct evidence or very stable repeated pattern

Confidence should go **down** when:
- evidence is sparse
- tools are weakly calibrated
- confounds are obvious
- alternative explanations remain live

---

## Interaction with skill 23 (evidence schema)

When evidence is already standardized:
- use the evidence objects directly in `evidence`
- inspect `conflicts` to trigger `conflict` reasoning
- do not rewrite or flatten the evidence first

Practical rule:
- **reason on top of EvidenceCollection, not instead of it**

---

## Minimal reasoning procedure

1. Pick the most important observation.
2. Choose reasoning mode: anomaly / conflict / trend / hypothesis.
3. Attach the real evidence.
4. Compare against a baseline.
5. Write 3–5 reasoning steps.
6. End with conclusion + action + prediction.
7. Assign honest confidence.

---

## Bottom line

This skill exists to make CE scientifically useful.

Not just:
- execute scripts
- dump tables
- repeat metrics

But:
- notice what matters
- explain why it matters
- turn evidence into the next experiment
