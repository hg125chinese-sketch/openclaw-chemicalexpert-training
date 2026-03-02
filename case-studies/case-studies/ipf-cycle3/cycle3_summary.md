# IPF / ALK5 (TGFBR1) Cycle 3 — Summary (Efficiency vs Rejection Sampling)

**Objective:** improve **efficiency** over Cycle 2’s rejection sampling while preserving the same hard KPI: produce candidates that yield a **hinge H-bond** in docked poses.

---

## Core question
Can we do better than rejection sampling?

- Rejection sampling reliably enforces constraints, but can be **sample-inefficient**.
- Cycle 3 asked whether we can keep (or improve) hinge-H-bond satisfaction while **accepting fewer wasted samples**.

---

## Attempt 1: Cross-attention conditioning (failed)
**Hypothesis:** a cross-attention conditioned VAE decoder would attend to the conditioning signal and generate hinge-binder-like chemistry more directly.

**What happened:**
- The model produced valid strings, but **conditioning did not translate into the required hinge-binding motifs** at the KPI level.
- In short: it looked fine by generic VAE diagnostics, yet still failed the domain KPI.

**Conclusion:** training-time architectural changes alone do not guarantee constraint satisfaction for this setting.

---

## Attempt 2: Strategy E — Logit bias decoding (succeeded)
**Idea:** keep the base generator, but enforce preferences at **inference time**.

Mechanism:
- During decoding, add a positive bias to the logits of tokens associated with the desired fragment/motif.
- This is a soft constraint, but applied at every decoding step; combined with validity checks, it becomes an effective practical constraint.

**Result:**
- Achieved ~**2× higher efficiency** than rejection sampling (higher accept rate for constraint-satisfying samples).
- Preserved the KPI outcome: reranked **Top 5** still had **4/5** molecules with a hinge H-bond.

---

## Architecture lesson learned
For string-based VAEs in this workflow, **hard domain KPIs are best enforced at inference time**:

- Conditioning can be ignored by an autoregressive decoder (especially when it finds an easier way to minimize loss).
- “Healthy” latent usage does not imply KPI satisfaction.
- Constraint-oriented decoding (e.g., logit bias, constrained sampling, or structured decoding) can align generation with project requirements with less retraining and better controllability.

---

## Takeaway
Cycle 3 establishes **Strategy E (logit bias decoding)** as the default recommended approach when:
- you need a must-have structural motif/interaction,
- you want efficiency better than rejection sampling,
- and you want a method that is easy to audit and tune.
