# IPF / ALK5 (TGFBR1) Cycle 2 — Summary (Cond. Generation Fix Attempts)

**Objective:** fix the Cycle 1 failure mode where the generator produced molecules that dock but do **not** reliably express the required **hinge H-bond** interaction pattern.

## Core problem statement
- In **Cycle 1**, the reranked **Top 5** docked candidates contained only **1/5** molecules with an explicit **hinge H-bond** (interaction-based check on docked poses).
- This indicates the pipeline was selecting the “best of a bad distribution”: docking alone could not compensate for a generator that drifted away from kinase-like hinge-binding chemistry.

## What we tried (Strategy A → D)
All strategies share the same intent: enforce/encourage a hinge-binder motif or fragment during generation and then verify outcomes with downstream filters and interaction checks.

### Results table (high-level)

| Strategy | Short description | Expected mechanism | What happened | Conclusion |
|---|---|---|---|---|
| **A** | Enrich training set with hinge-like examples before training | Shift the learned distribution toward hinge binders | No effective enrichment occurred (no new examples added under the chosen similarity/threshold constraints) | Not actionable under current constraints |
| **B** | Conditioning strength / target-rate tuning | Make conditioning “matter” by increasing its weight | At the tested target rate, training distribution did not meaningfully change | Conditioning signal too weak / easy to ignore |
| **C** | Fragment-prefix / scaffold-prefix conditioning | Force the model to start from a desired substructure | Works as a **trivial solution** when a literal prefix is provided, but fails when prefix is removed or scaled up | Not a robust conditioning method |
| **D** | **Rejection sampling** (generate → filter → keep only those that satisfy hinge constraints) | Enforce hard constraints at sampling time | Constraint satisfaction becomes reliable at the cost of sample efficiency | **Adopted as the Cycle 2 solution** |

## Key finding
**A SELFIES GRU VAE can structurally ignore conditioning.**

Even when latent space diagnostics look “healthy” (e.g., not fully collapsed), a string VAE decoder can learn to reconstruct/generate well while **treating conditioning as optional**, especially when the conditioning is not enforced as a hard constraint.

Practically:
- Training-time conditioning can be silently bypassed.
- Post-hoc selection (interaction checks) reveals the failure.

## Final approach (Cycle 2)
**Strategy D: rejection sampling**

Pipeline:
1. Sample many candidates from the trained generator.
2. Apply strict rules-first filtering (developability + project safety constraints).
3. Dock survivors.
4. Run interaction-based checks on top-ranked docked poses.
5. Keep only candidates that satisfy the hinge H-bond requirement.

## Outcome
- Cycle 2 reranked **Top 5**: **4/5** molecules show an explicit **hinge H-bond**.

## Takeaway
If the KPI is a **must-have structural interaction** (hinge H-bond), then **sampling-time constraints** (even if inefficient) can outperform “soft” training-time conditioning in string VAEs.
