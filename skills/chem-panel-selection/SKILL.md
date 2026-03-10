---
name: chem-panel-selection
description: Standardized panel selection protocol for small molecule shortlist (3–10) under conflicting signals (Vina, ProLIF hinge H-bond, Boltz-2 affinity, QE DFT). Defines hard gates vs ranking signals, disagreement analysis (2x2 grid), and auditable decision recording.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧭", "requires": { "bins": ["python3"], "python": ["pandas", "numpy"] } } }
---

# chem-panel-selection — conflict-aware selection (panel, not a single winner)

At the end of a DMTA cycle you usually have a shortlist (often **3–10** molecules) and **multiple signals that disagree**:

- **Docking** (Vina): fast, pose-dependent, brittle
- **Pose evidence** (ProLIF): hinge H-bond and key residues, but depends on docking pose quality
- **Affinity prediction** (Boltz‑2): orthogonal ML signal, can disagree with docking
- **QC** (QE DFT): feasibility gate, but expensive and sparse

If you pick “the best” by one metric, you will overfit noise.

This skill gives a **standard protocol** to select a **panel** that is:
- auditable (explicit rationale)
- disagreement-aware (tests competing hypotheses)
- robust (doesn’t collapse to one metric)

---

## Workspace variables

```bash
OPENCLAW_WORKSPACE="/home/node/.openclaw/workspace-chemicalexpert"
QMD="/home/node/.openclaw/.npm-global/bin/qmd"
CHEM_PY="/opt/conda/envs/chem/bin/python"
```

---

## When to use

Use this when:
- You have a shortlist **N=3–10** (Top5/Top10) and need to decide what to send to:
  - QE DFT QC, or
  - retrosynthesis / synthesis planning, or
  - experimental validation.
- Vina and Boltz‑2 disagree (rank correlation is weak / flips occur).
- You need a **documented decision** that can be audited later.

Do not use this for large libraries (N>100). Use hard gates + ranking + diversity selection upstream.

---

## Signals: hard gates vs ranking signals

### Hard gates (binary; non-negotiable)

Typical hard gates (project-dependent):
- **Hinge H-bond = True** (ProLIF-validated) for Type I kinase inhibitors
- **DFT QC PASS** (QE) if the next stage requires QC-feasible chemistry
- Project safety denylist: **PASS** (no hard reject)

Rule:
- If a molecule fails a hard gate, it **cannot** be selected for the next stage that assumes that gate.

### Ranking signals (continuous; noisy)

Use these to rank **within** the gate-passing set:
- Vina score (more negative is better)
- Boltz‑2 affinity (more favorable log10(IC50 µM) is better) + binder_prob
- Interaction score / key residue coverage
- Simple multi-objective score (auditable, not “truth”)

Rule:
- Ranking signals are for **comparative decisions**, not claims of absolute affinity.

---

## Panel selection protocol (default)

### Step 0 — Define the decision context
Write down:
- What is the next action? (QE DFT / synthesis planning / wet lab)
- What are the hard gates for that action?
- How many molecules can you afford to take forward? (panel size K)

### Step 1 — Apply hard gates
Filter the shortlist.

### Step 2 — Normalize and label ranks (no over-engineering)
Within the gate-passing set:
- compute rank by Vina (best = 1)
- compute rank by Boltz‑2 affinity (best = 1)

### Step 3 — Disagreement analysis (2x2 grid)
Split by “high/low” on each ranking signal (median split is fine):

- **High Vina / High Boltz**: consensus winners
- **High Vina / Low Boltz**: docking-favored hypothesis
- **Low Vina / High Boltz**: ML-favored hypothesis
- **Low Vina / Low Boltz**: deprioritize unless diversity/novelty is the goal

### Step 4 — Choose a panel (K molecules)
Default selection heuristic:
- 1–2 from **High/High** (if exists)
- 1 from **High Vina / Low Boltz**
- 1 from **Low Vina / High Boltz**
- optional 1 “wildcard” for diversity/novelty or mechanistic hypothesis

### Step 5 — Record an auditable rationale
For each selected molecule:
- list gate status
- list key signal values
- explain why it was chosen (which hypothesis it tests)

Avoid:
- “picked because best score” with no context
- changing gates after seeing results (p-hacking)

---

## Report template (copy/paste)

Use this structure in your cycle report:

```markdown
## Panel selection (conflict-aware)

### Context
- next action: <QE DFT / synthesis / etc>
- shortlist size: N=<>
- panel size: K=<>

### Hard gates
- hinge_hbond required: <yes/no>
- DFT PASS required: <yes/no>
- safety hard reject: <yes/no>

### Disagreement grid (Vina vs Boltz-2)
- split rule: median
- summary counts per quadrant

### Selected panel (K)
For each:
- gates: hinge=<>, dft=<>, safety=<>
- signals: vina=, interaction_score=, boltz_affinity=, binder_prob=
- rationale: <one sentence>

### What we learn
- If High Vina / Low Boltz wins experimentally → docking signal stronger in this regime
- If Low Vina / High Boltz wins → affinity predictor is adding value
```

---

## Minimal helper (script-like) to build the 2x2 grid

```python
#!/opt/conda/envs/chem/bin/python
"""Make a Vina vs Boltz-2 2x2 grid for a shortlist.

Expected columns (subset is OK):
- mol_id
- vina_score
- boltz_affinity  (log10(IC50 uM); more negative = better)
- binder_prob
- has_hinge_hbond (bool)
- qc_flag (PASS/OPT_FAIL)
"""

import pandas as pd


def add_quadrants(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    # Higher-is-better transforms
    out["vina_rank_signal"] = -out["vina_score"].astype(float)
    out["boltz_rank_signal"] = -out["boltz_affinity"].astype(float)

    v_med = out["vina_rank_signal"].median()
    b_med = out["boltz_rank_signal"].median()

    out["vina_hi"] = out["vina_rank_signal"] >= v_med
    out["boltz_hi"] = out["boltz_rank_signal"] >= b_med

    def _q(r):
        return ("HighVina" if r.vina_hi else "LowVina") + "/" + ("HighBoltz" if r.boltz_hi else "LowBoltz")

    out["quadrant"] = out.apply(_q, axis=1)
    return out


def quadrant_counts(df: pd.DataFrame) -> dict:
    return df["quadrant"].value_counts().to_dict()
```

---

## Anti-overfitting rules

- Do not tune weights/thresholds on a single cycle’s outcome.
- Keep gates stable across cycles unless a new failure mode forces a change.
- Prefer **panels** over “one best molecule” until you have real experimental feedback.
