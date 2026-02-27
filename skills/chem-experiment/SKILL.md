---
name: chem-experiment
description: Autonomously plan, execute, and analyze in-silico DMTA (Design-Make-Test-Analyze) cycles. Orchestrates all other chem skills into a coherent experimental workflow — from target selection through candidate delivery. Handles multi-objective optimization, iteration planning, and decision-making under uncertainty.
homepage: https://openclaw.ai
metadata: { "openclaw": { "emoji": "🔄", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "numpy", "pandas"] } } }
---

# Autonomous Experiment Planning

Orchestrate complete DMTA (Design-Make-Test-Analyze) cycles by combining all 11 chem skills into a coherent workflow. This is the "conductor" skill — it decides what to run, in what order, and how to iterate based on results.

## When to Use

- User says "find me drug candidates for target X" (end-to-end request)
- User wants to iterate: "these 5 candidates weren't great, what next?"
- User needs a research plan with timeline and milestones
- User asks "what should I do next?" in a drug design project
- Any task that requires coordinating 3+ chem skills in sequence

## Core Philosophy

1. **Plan before computing.** Write the experimental plan, get user approval, THEN execute. Wasted GPU time is worse than wasted planning time.
2. **Each cycle must produce a decision.** A DMTA cycle that ends with "more data needed" without specifying WHAT data is a failed cycle.
3. **Multi-objective optimization is the reality.** No molecule is "best" — there are trade-offs (activity vs selectivity vs ADMET vs synthesis). Make trade-offs explicit.
4. **Iterate, don't restart.** Each cycle builds on the last. Carry forward the SAR model, the ADMET failures, the synthesis lessons.
5. **The user is the wet-lab chemist.** Our in-silico cycles inform their real experiments. Prioritize actionable recommendations over exhaustive analyses.

## Phase 1: Project Initialization

### 1.1 Target Assessment Checklist

```python
#!/opt/conda/envs/chem/bin/python
"""Project initialization: assess feasibility before starting."""

def assess_target(target_name):
    """Run through the target assessment checklist.

    Returns a structured assessment that guides skill selection.
    """
    assessment = {
        "target": target_name,
        "questions": [
            {
                "q": "Is there a crystal structure with a co-crystallized ligand?",
                "options": ["yes_good_resolution", "yes_poor_resolution",
                            "homology_model_only", "no_structure"],
                "impact": {
                    "yes_good_resolution": "Enable docking (skill 9), pocket-conditioned 3D gen (skill 10)",
                    "yes_poor_resolution": "Docking possible but less reliable, redocking validation critical",
                    "homology_model_only": "Docking with caution, rank-only (no score cutoffs)",
                    "no_structure": "Skip docking, rely on ligand-based methods (skills 2,3,5)",
                },
            },
            {
                "q": "How many known actives with measured affinity?",
                "options": ["0-50", "50-500", "500-5000", "5000+"],
                "impact": {
                    "0-50": "Transfer learning or few-shot; limited QSAR. Literature-heavy (skill 1)",
                    "50-500": "RF baseline feasible (skill 2); GNN likely overfits (skill 5)",
                    "500-5000": "Sweet spot for RF+GNN comparison; scaffold split meaningful",
                    "5000+": "Full modeling suite; can afford scaffold+temporal splits",
                },
            },
            {
                "q": "Is this a CNS target (needs to cross BBB)?",
                "options": ["yes", "no", "unknown"],
                "impact": {
                    "yes": "TPSA < 90, LogP 1-3, MW < 450 are hard constraints in ADMET (skill 6)",
                    "no": "Standard Lipinski/Veber sufficient",
                    "unknown": "Apply both CNS and standard filters, report both",
                },
            },
            {
                "q": "What's the project goal?",
                "options": ["hit_finding", "lead_optimization", "scaffold_hopping",
                            "selectivity_optimization"],
                "impact": {
                    "hit_finding": "Broad generation (skill 3) + diverse screening",
                    "lead_optimization": "Focused generation around known scaffolds + ADMET optimization",
                    "scaffold_hopping": "3D generation (skill 10) + docking validation",
                    "selectivity_optimization": "Multi-target QSAR + differential docking",
                },
            },
        ],
    }
    return assessment

def print_assessment(assessment):
    """Pretty-print target assessment."""
    print(f"\n{'='*60}")
    print(f"Target Assessment: {assessment['target']}")
    print(f"{'='*60}")
    for i, q in enumerate(assessment["questions"]):
        print(f"\n{i+1}. {q['q']}")
        for opt, impact in q["impact"].items():
            print(f"   [{opt}] → {impact}")
```

### 1.2 Skill Selection Matrix

```python
def select_skills(structure_available, n_actives, is_cns, goal):
    """Select which skills to use based on target assessment.

    Returns ordered list of skills for the DMTA cycle.
    """
    skills = []

    # Always start with literature
    skills.append(("chem-literature", "Gather prior knowledge, known scaffolds, SAR"))

    # Modeling
    if n_actives >= 50:
        skills.append(("chem-qsar", "RF + Morgan baseline"))
    if n_actives >= 500:
        skills.append(("chem-gnn", "GCN comparison to RF"))

    # Generation
    if goal in ("hit_finding", "scaffold_hopping"):
        skills.append(("chem-molgen", "SELFIES VAE for diverse candidates"))
    if goal == "scaffold_hopping" and structure_available:
        skills.append(("chem-3dgen", "3D diffusion for shape-based design"))

    # Filtering
    skills.append(("chem-admet", "Rules-first ADMET filter"))

    # Structure-based
    if structure_available:
        skills.append(("chem-docking", "Vina docking against crystal structure"))
        skills.append(("chem-mlff", "MACE refinement of docking poses"))

    # Synthesis
    skills.append(("chem-retrosynthesis", "Route planning + SA score"))
    skills.append(("chem-rxn-conditions", "Conditions for each synthesis step"))

    # Reporting
    skills.append(("chem-llm", "Generate medchem descriptions and analysis"))

    return skills
```

## Phase 2: DMTA Cycle Execution

### 2.1 Cycle Plan Template

```markdown
# DMTA Cycle [N] Plan
**Target**: [name]
**Date**: [date]
**Goal**: [what this cycle aims to achieve]
**Previous cycle learnings**: [if applicable]

## Design Phase
- [ ] Literature review (skill 1): [specific queries]
- [ ] Model update (skills 2/5): [retrain or reuse?]
- [ ] Generate candidates (skill 3/10): [method, N molecules, constraints]

## Make Phase (in silico)
- [ ] ADMET filter (skill 6): [which rules, thresholds]
- [ ] Docking (skill 9): [PDB ID, box, exhaustiveness]
- [ ] Pose refinement (skill 11): [MACE strain check]
- [ ] Retrosynthesis (skill 4): [SA threshold]
- [ ] Conditions (skill 7): [annotate top candidates]

## Test Phase (scoring)
- [ ] Rank by: [scoring function: RF_pred × QED × (1/SA) × dock_score]
- [ ] Gate criteria: [SA ≤ 6, |RF−GNN| ≤ 1.0, QED > 0.3, dock < -7.0]

## Analyze Phase
- [ ] Top N candidates with full profiles
- [ ] SAR analysis: what worked, what didn't
- [ ] Decision: next cycle direction

## Expected Output
- Final report: `research/ai4chem/closed_loop/[target]/cycle_[N]/report.md`
- Candidate pool: `...cycle_[N]/candidates.csv`
- Top 5 with routes: `...cycle_[N]/top_5.md`
```

### 2.2 Multi-Objective Scoring

```python
import numpy as np

def multi_objective_score(candidate, weights=None):
    """Score a candidate across multiple objectives.

    Default weights balance activity, drug-likeness, synthesizability, and binding.
    All components normalized to [0, 1] where 1 = best.
    """
    if weights is None:
        weights = {
            "activity": 0.30,    # predicted pChEMBL
            "qed": 0.20,         # drug-likeness
            "sa_norm": 0.20,     # synthesizability (inverted: low SA = high score)
            "dock_norm": 0.15,   # docking score (normalized)
            "consistency": 0.15, # model agreement (low RF-GNN diff)
        }

    scores = {}

    # Activity: normalize pChEMBL to [0,1] assuming range [4, 10]
    pred = candidate.get("rf_pred", 6.0)
    scores["activity"] = np.clip((pred - 4.0) / 6.0, 0, 1)

    # QED: already [0, 1]
    scores["qed"] = candidate.get("qed", 0.5)

    # SA: invert and normalize (SA 1=easy → 1.0, SA 10=hard → 0.0)
    sa = candidate.get("sa_score", 5.0)
    scores["sa_norm"] = np.clip(1.0 - (sa - 1.0) / 9.0, 0, 1)

    # Docking: normalize (e.g., -12 to -4 range)
    dock = candidate.get("dock_score")
    if dock is not None:
        scores["dock_norm"] = np.clip((-dock - 4.0) / 8.0, 0, 1)
    else:
        scores["dock_norm"] = 0.5  # neutral if no docking

    # Consistency: low |RF-GNN| = high score
    diff = candidate.get("rf_gnn_diff", 0.5)
    scores["consistency"] = np.clip(1.0 - diff / 2.0, 0, 1)

    # Weighted sum
    total = sum(weights[k] * scores.get(k, 0.5) for k in weights)

    return {
        "total_score": round(total, 4),
        "components": {k: round(v, 4) for k, v in scores.items()},
        "weights": weights,
    }
```

### 2.3 Gate Criteria (Hard Filters)

```python
def apply_gates(candidates, gates=None):
    """Apply hard gates to filter candidates.

    Gates are non-negotiable. A candidate that fails ANY gate is eliminated.
    """
    if gates is None:
        gates = {
            "sa_max": 6.0,
            "qed_min": 0.3,
            "rf_gnn_diff_max": 1.0,
            "lipinski_violations_max": 1,
            "motif_denylist": ["C=C=N", "[SH]"],  # from DRD2 lesson
        }

    passed = []
    rejected = {"sa": 0, "qed": 0, "consistency": 0, "lipinski": 0, "motif": 0}

    for c in candidates:
        # SA gate
        if c.get("sa_score", 10) > gates["sa_max"]:
            rejected["sa"] += 1
            continue
        # QED gate
        if c.get("qed", 0) < gates["qed_min"]:
            rejected["qed"] += 1
            continue
        # Consistency gate
        if c.get("rf_gnn_diff", 999) > gates["rf_gnn_diff_max"]:
            rejected["consistency"] += 1
            continue
        # Lipinski gate
        if c.get("lipinski_violations", 0) > gates["lipinski_violations_max"]:
            rejected["lipinski"] += 1
            continue
        # Motif gate
        from rdkit import Chem
        mol = Chem.MolFromSmiles(c.get("smiles", ""))
        if mol:
            motif_hit = False
            for pattern in gates.get("motif_denylist", []):
                if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                    motif_hit = True
                    break
            if motif_hit:
                rejected["motif"] += 1
                continue

        passed.append(c)

    print(f"\nGate Results: {len(passed)} passed / {len(candidates)} total")
    for gate, count in rejected.items():
        if count > 0:
            print(f"  Rejected by {gate}: {count}")

    return passed, rejected
```

## Phase 3: Iteration Planning

### 3.1 Analyze Cycle Results & Plan Next

```python
def analyze_and_plan_next(cycle_results, cycle_number):
    """Analyze current cycle results and generate plan for next cycle.

    This is the "Analyze" phase of DMTA — the most important phase.
    """
    analysis = {
        "cycle": cycle_number,
        "summary": {},
        "sar_insights": [],
        "failures": [],
        "next_cycle_plan": {},
    }

    candidates = cycle_results.get("top_candidates", [])
    pool = cycle_results.get("candidate_pool", [])

    # Summary statistics
    analysis["summary"] = {
        "pool_size": len(pool),
        "passed_gates": len(candidates),
        "mean_score": np.mean([c.get("total_score", 0) for c in candidates]) if candidates else 0,
        "score_range": (
            min(c.get("total_score", 0) for c in candidates) if candidates else 0,
            max(c.get("total_score", 0) for c in candidates) if candidates else 0,
        ),
    }

    # SAR insights (what structural features correlate with high scores?)
    # This would be extracted from the data in a real implementation
    analysis["sar_insights"] = [
        "Analyze: which scaffolds appear in top 10% vs bottom 10%?",
        "Analyze: which substituents improve activity while maintaining ADMET?",
        "Analyze: are there cluster of similar high-scorers (exploitation opportunity)?",
    ]

    # Common failure modes
    rejected = cycle_results.get("rejected", {})
    if rejected.get("sa", 0) > len(pool) * 0.3:
        analysis["failures"].append("Too many candidates fail SA gate → constrain generation")
    if rejected.get("qed", 0) > len(pool) * 0.3:
        analysis["failures"].append("Too many candidates fail QED → add QED to generation objective")
    if rejected.get("consistency", 0) > len(pool) * 0.3:
        analysis["failures"].append("High model disagreement → retrain with more data or use ensemble")

    # Next cycle plan
    if not candidates:
        analysis["next_cycle_plan"] = {
            "action": "RESTART with relaxed gates or different scaffold",
            "priority": "Investigate why no candidates passed all gates",
        }
    elif analysis["summary"]["mean_score"] < 0.5:
        analysis["next_cycle_plan"] = {
            "action": "OPTIMIZE: focused generation around best scaffolds from this cycle",
            "priority": "Improve activity predictions with expanded training data",
        }
    else:
        analysis["next_cycle_plan"] = {
            "action": "REFINE: narrow around top hits, improve ADMET/synthesis",
            "priority": "Move to lead optimization phase",
        }

    return analysis
```

### 3.2 Iteration History

```python
class DMTAProject:
    """Track multi-cycle DMTA project state."""

    def __init__(self, target_name, project_dir):
        self.target = target_name
        self.project_dir = project_dir
        self.cycles = []
        self.all_candidates = []  # cumulative
        self.sar_model = None     # evolves across cycles

    def start_cycle(self, cycle_plan):
        """Begin a new DMTA cycle."""
        n = len(self.cycles) + 1
        cycle = {
            "number": n,
            "plan": cycle_plan,
            "status": "in_progress",
            "results": None,
            "analysis": None,
        }
        self.cycles.append(cycle)
        print(f"\n{'='*60}")
        print(f"Starting DMTA Cycle {n} for {self.target}")
        print(f"Goal: {cycle_plan.get('goal', 'unspecified')}")
        print(f"{'='*60}")
        return cycle

    def complete_cycle(self, results):
        """Complete current cycle with results and analysis."""
        cycle = self.cycles[-1]
        cycle["results"] = results
        cycle["status"] = "complete"

        # Analyze
        cycle["analysis"] = analyze_and_plan_next(results, cycle["number"])

        # Accumulate candidates
        self.all_candidates.extend(results.get("top_candidates", []))

        return cycle["analysis"]

    def get_project_summary(self):
        """Generate project-level summary across all cycles."""
        summary = {
            "target": self.target,
            "n_cycles": len(self.cycles),
            "total_candidates_evaluated": sum(
                len(c.get("results", {}).get("candidate_pool", []))
                for c in self.cycles
            ),
            "total_top_candidates": len(self.all_candidates),
            "cycles": [],
        }
        for c in self.cycles:
            summary["cycles"].append({
                "number": c["number"],
                "goal": c["plan"].get("goal"),
                "status": c["status"],
                "n_passed": len(c.get("results", {}).get("top_candidates", [])),
            })
        return summary
```

## Phase 4: Report Generation

### 4.1 Final Deliverable Template

```markdown
# DMTA Cycle [N] Report: [Target]
**Date**: [date]
**Status**: Complete

## Executive Summary
[2-3 sentences: what was done, key finding, recommendation]

## Candidates (Top 5)
| Rank | SMILES | Activity | QED | SA | Dock | Score | Route |
|------|--------|----------|-----|----|------|-------|-------|
| 1 | ... | ... | ... | ... | ... | ... | ... |

## SAR Insights
[What structural features drive activity? What fails?]

## ADMET Profiles
[Traffic light table for top 5]

## Synthesis Routes
[Conditions table for top 5]

## Decision for Next Cycle
[What to do next and why]

## Reproducibility
```bash
# Full pipeline
cd /home/node/.openclaw/workspace-chemicalexpert
python research/ai4chem/closed_loop/[target]/cycle_[N]/run_cycle.py
```
```

### 4.2 Skill Usage Log

Always document which skills were used and how:

```markdown
## Skill Usage
| Skill | Purpose | Key Result |
|-------|---------|------------|
| chem-literature | DRD2 prior knowledge | 5 known scaffold families identified |
| chem-qsar | RF baseline | R²=0.72 on scaffold split |
| chem-gnn | GCN comparison | Lost to RF → used as consistency check |
| chem-molgen | SELFIES VAE | 4984 candidates generated |
| chem-admet | Rules filter | 3304 passed Lipinski+Veber+QED |
| chem-retrosynthesis | Route planning | 2035 have single-step routes |
| chem-rxn-conditions | Conditions | EDC/HOBt + DIPEA for amide formation |
| chem-docking | 6CM4 docking | Top 1 scores -8.48 kcal/mol |
| chem-mlff | Pose refinement | Strain < 3 kcal/mol for top 3 |
| chem-llm | Medchem description | Benzothiazole-piperazine SAR context |
```

## Checklist Before Reporting

- [ ] **Plan approved**: User confirmed cycle plan before execution?
- [ ] **All gates documented**: Hard filters and thresholds stated upfront?
- [ ] **Multi-objective scores**: Not just one metric — balanced scorecard?
- [ ] **Iteration history**: This cycle connected to previous cycles?
- [ ] **SAR extracted**: Not just candidates — learnings for next cycle?
- [ ] **Actionable recommendation**: Clear "do this next" statement?
- [ ] **Reproducible**: Script path + git commit in report?
- [ ] **Skill attribution**: Which skills used, in what order, why?
