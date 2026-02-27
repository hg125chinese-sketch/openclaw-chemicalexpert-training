---
name: chem-kinase-sar
description: Kinase-specific SAR prior knowledge, hinge binder fragment libraries, and scaffold coverage KPIs. Provides SMARTS-based classification of kinase-relevant motifs, coverage statistics for candidate pools, and hard gates for kinase-targeted campaigns.
homepage: https://www.ebi.ac.uk/chembl/
metadata: { "openclaw": { "emoji": "🔑", "requires": { "bins": ["python3"], "python": ["rdkit", "numpy", "pandas"] } } }
---

# Kinase SAR: Hinge Binders, Scaffold Coverage & Chemical Space KPIs

Kinase inhibitors are the largest class of targeted cancer and anti-inflammatory drugs. Their design is constrained by well-understood SAR — especially the hinge region interaction. This skill turns kinase medicinal chemistry knowledge into computable, auditable gates.

## When to Use

- Generating molecules targeting **any kinase** (ALK5, EGFR, BRAF, JAK, CDK, etc.)
- Evaluating whether a candidate pool "looks like kinase inhibitors"
- Setting coverage KPIs for generative models in kinase campaigns
- Filtering candidates by kinase-relevant structural features
- Analyzing training set composition for kinase QSAR models

## Core Philosophy

1. **Kinase SAR is not optional knowledge — it's a hard filter.** Decades of kinase drug discovery have established which motifs work and which don't. Ignoring this prior is like ignoring Lipinski rules for oral drugs.
2. **Hinge binding is the minimum bar.** Most kinase inhibitors form 1-3 hydrogen bonds with the hinge region. A molecule without a plausible hinge-binding motif is overwhelmingly unlikely to be a potent kinase inhibitor.
3. **Coverage is a generative model KPI.** If your generator produces <5% hinge-binder-containing molecules when targeting a kinase, the generator is broken — no amount of docking will fix it.
4. **Classification before scoring.** First ask "is this a plausible kinase inhibitor?" (binary), then ask "how good?" (scoring). Don't score molecules that fail the plausibility check.
5. **Target-class knowledge generalizes.** The hinge motif library applies across the kinome. Target-specific SAR (DFG-in/out, gatekeeper interactions) is an overlay, not a replacement.

## Phase 1: Hinge Binder Fragment Library

### 1.1 Core Hinge Binder SMARTS

These fragments represent the most common kinase hinge-binding motifs found across approved kinase inhibitors and clinical candidates. Each provides 1-3 H-bonds to the hinge backbone (typically to the NH and C=O of the +1 and +3 residues).

```python
#!/opt/conda/envs/chem/bin/python
"""Kinase hinge binder SMARTS library.

Each entry: (name, SMARTS, description, example_drugs)
SMARTS are designed to match the core hinge-binding fragment,
not the full molecule.
"""
from rdkit import Chem

HINGE_BINDERS = [
    # --- Pyrimidine/Aminopyrimidine family (most common) ---
    (
        "2-aminopyrimidine",
        "[NH1,NH2]-c1ncccn1",
        "Classic hinge binder. 2 H-bonds (NH donor + ring N acceptor). "
        "Found in imatinib-like scaffolds.",
        ["imatinib", "dasatinib"]
    ),
    (
        "4-aminopyrimidine",
        "[NH1,NH2]-c1ccncn1",
        "Aminopyrimidine with NH at 4-position. Common in CDK/kinase inhibitors.",
        ["palbociclib", "ribociclib"]
    ),
    (
        "2,4-diaminopyrimidine",
        "[NH1,NH2]-c1nc([NH1,NH2])ccn1",
        "Dual NH donors. Strong hinge engagement.",
        ["trimethoprim-like"]
    ),
    (
        "pyrimidine_ring",
        "c1ccncn1",
        "Bare pyrimidine ring (broader match, less specific).",
        []
    ),

    # --- Pyridine family ---
    (
        "2-aminopyridine",
        "[NH1,NH2]-c1ccccn1",
        "Single hinge H-bond donor + ring N acceptor.",
        ["crizotinib"]
    ),
    (
        "3-aminopyridine",
        "[NH1,NH2]-c1cccnc1",
        "NH at 3-position.",
        []
    ),

    # --- Pyrazole family ---
    (
        "pyrazole_nh",
        "[nH]1ccnc1",
        "NH pyrazole — strong hinge H-bond donor. Very common in kinase inhibitors.",
        ["crizotinib"]
    ),
    (
        "pyrazolopyrimidine",
        "c1cnc2[nH]ncc2n1",
        "Fused pyrazole-pyrimidine. Privileged kinase scaffold.",
        ["ibrutinib"]
    ),
    (
        "pyrazole_ring",
        "c1cc[nH]n1",
        "Any NH-pyrazole (broader match).",
        []
    ),

    # --- Quinazoline family ---
    (
        "4-aminoquinazoline",
        "[NH1,NH2]-c1ncnc2ccccc12",
        "Anilinoquinazoline — the EGFR inhibitor scaffold.",
        ["gefitinib", "erlotinib", "lapatinib"]
    ),
    (
        "quinazoline_ring",
        "c1ncnc2ccccc12",
        "Bare quinazoline ring.",
        []
    ),

    # --- Imidazole / Benzimidazole ---
    (
        "imidazole_nh",
        "[nH]1ccnc1",
        "NH-imidazole as hinge binder.",
        []
    ),
    (
        "benzimidazole",
        "c1ccc2[nH]cnc2c1",
        "Benzimidazole — used in some kinase inhibitors.",
        ["dovitinib"]
    ),

    # --- Indazole ---
    (
        "indazole_nh",
        "c1ccc2[nH]ncc2c1",
        "Indazole NH — strong hinge donor.",
        ["axitinib", "linifanib"]
    ),

    # --- Pyrrolopyrimidine ---
    (
        "7-deazapurine",
        "c1cnc2[nH]ccc2n1",
        "7-deazapurine (pyrrolopyrimidine). Privileged kinase scaffold.",
        ["tofacitinib"]
    ),

    # --- Oxindole ---
    (
        "oxindole",
        "O=C1Cc2ccccc2N1",
        "2-oxindole — hinge binder via lactam NH + C=O.",
        ["sunitinib", "nintedanib"]
    ),

    # --- Urea/Thiourea (Type II) ---
    (
        "diaryl_urea",
        "c-[NH]-C(=O)-[NH]-c",
        "Diaryl urea — extends from hinge to DFG (Type II).",
        ["sorafenib", "regorafenib"]
    ),

    # --- Quinoline ---
    (
        "aminoquinoline",
        "[NH1,NH2]-c1ccnc2ccccc12",
        "Aminoquinoline hinge binder.",
        ["bosutinib"]
    ),
]


def compile_hinge_library():
    """Compile SMARTS patterns. Returns list of (name, mol_pattern, desc, examples)."""
    compiled = []
    for name, smarts, desc, examples in HINGE_BINDERS:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            print(f"WARNING: Failed to compile SMARTS for {name}: {smarts}")
            continue
        compiled.append((name, pat, desc, examples))
    return compiled


HINGE_LIBRARY = compile_hinge_library()
```

### 1.2 Classify a Single Molecule

```python
def classify_hinge_binders(smi, library=None):
    """Check which hinge binder motifs a molecule contains.

    Returns:
        dict with keys:
        - smiles: input SMILES
        - valid: bool
        - hits: list of (name, n_matches)
        - has_hinge_binder: bool (at least one hit)
        - n_distinct_motifs: int
    """
    if library is None:
        library = HINGE_LIBRARY

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "hits": [],
                "has_hinge_binder": False, "n_distinct_motifs": 0}

    hits = []
    for name, pat, desc, examples in library:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            hits.append((name, len(matches)))

    return {
        "smiles": smi,
        "valid": True,
        "hits": hits,
        "has_hinge_binder": len(hits) > 0,
        "n_distinct_motifs": len(hits),
    }
```

## Phase 2: Pool Coverage Analysis

### 2.1 Coverage Statistics

```python
import pandas as pd
import numpy as np


def pool_coverage_analysis(smiles_list, library=None, label="pool"):
    """Analyze hinge binder coverage across a molecule pool.

    Args:
        smiles_list: list of SMILES strings
        library: compiled hinge library (default: HINGE_LIBRARY)
        label: name for this pool (e.g., "candidate_pool", "training_set")

    Returns:
        dict with:
        - summary: overall stats
        - per_motif: DataFrame with per-motif coverage
        - molecules: per-molecule classification
    """
    if library is None:
        library = HINGE_LIBRARY

    results = [classify_hinge_binders(smi, library) for smi in smiles_list]

    n_total = len(results)
    n_valid = sum(1 for r in results if r["valid"])
    n_with_hinge = sum(1 for r in results if r["has_hinge_binder"])

    # Per-motif stats
    motif_counts = {}
    for name, _, _, _ in library:
        motif_counts[name] = 0
    for r in results:
        for name, count in r["hits"]:
            motif_counts[name] += 1

    per_motif = pd.DataFrame([
        {"motif": name, "count": motif_counts[name],
         "pct": 100.0 * motif_counts[name] / n_valid if n_valid > 0 else 0}
        for name in motif_counts
    ]).sort_values("count", ascending=False)

    summary = {
        "label": label,
        "n_total": n_total,
        "n_valid": n_valid,
        "n_with_any_hinge_binder": n_with_hinge,
        "pct_with_hinge_binder": 100.0 * n_with_hinge / n_valid if n_valid > 0 else 0,
        "top3_motifs": per_motif.head(3)[["motif", "pct"]].to_dict("records"),
    }

    return {
        "summary": summary,
        "per_motif": per_motif,
        "molecules": results,
    }


def compare_pools(pool_a_smiles, pool_b_smiles, label_a="training", label_b="generated"):
    """Compare hinge binder coverage between two pools.

    Typical use: compare training set vs generated candidates.
    Red flag: if generated pool has much lower coverage than training set,
    the generator is not exploring kinase-relevant chemical space.

    Returns:
        dict with comparison table and diagnosis
    """
    a = pool_coverage_analysis(pool_a_smiles, label=label_a)
    b = pool_coverage_analysis(pool_b_smiles, label=label_b)

    # Merge per-motif tables
    merged = a["per_motif"][["motif", "pct"]].rename(columns={"pct": f"pct_{label_a}"})
    merged = merged.merge(
        b["per_motif"][["motif", "pct"]].rename(columns={"pct": f"pct_{label_b}"}),
        on="motif", how="outer"
    ).fillna(0)
    merged["ratio"] = merged[f"pct_{label_b}"] / merged[f"pct_{label_a}"].replace(0, np.nan)

    # Diagnosis
    overall_ratio = (b["summary"]["pct_with_hinge_binder"] /
                     a["summary"]["pct_with_hinge_binder"]
                     if a["summary"]["pct_with_hinge_binder"] > 0 else 0)

    if overall_ratio >= 0.5:
        diagnosis = "ACCEPTABLE: Generated pool has reasonable kinase motif coverage"
    elif overall_ratio >= 0.1:
        diagnosis = "WARNING: Generated pool has significantly lower kinase motif coverage. Consider scaffold-conditioned generation."
    else:
        diagnosis = "CRITICAL: Generated pool has almost no kinase motif coverage. Generator is exploring wrong chemical space. Fix generation strategy before filtering."

    return {
        "summary_a": a["summary"],
        "summary_b": b["summary"],
        "comparison": merged,
        "overall_ratio": overall_ratio,
        "diagnosis": diagnosis,
    }
```

### 2.2 Coverage KPI Gate

```python
def kinase_coverage_gate(pool_smiles, min_coverage_pct=5.0, label="candidate_pool"):
    """Hard gate: does the pool have enough kinase-relevant molecules?

    Args:
        pool_smiles: list of SMILES
        min_coverage_pct: minimum % of molecules with any hinge binder motif
            Default 5% is deliberately low — even this catches total generator failure.
            For a focused kinase campaign, expect 20-50%.

    Returns:
        dict with pass/fail and stats
    """
    result = pool_coverage_analysis(pool_smiles, label=label)
    pct = result["summary"]["pct_with_hinge_binder"]

    return {
        "gate": "kinase_coverage",
        "pass": pct >= min_coverage_pct,
        "pct_with_hinge_binder": pct,
        "threshold": min_coverage_pct,
        "n_with_hinge": result["summary"]["n_with_any_hinge_binder"],
        "n_total": result["summary"]["n_valid"],
        "action_if_fail": "Do NOT proceed to docking. Fix generator first: "
                          "scaffold-conditioned VAE, fragment-constrained generation, "
                          "or augment training data with known kinase actives.",
    }
```

## Phase 3: Kinase Inhibitor Type Classification

### 3.1 Type I / II / III / IV / V Classification

```python
# DFG-motif and allosteric pocket indicators
KINASE_TYPE_INDICATORS = {
    "type_I_hinge": [
        # Occupies ATP pocket, interacts with hinge in DFG-in conformation
        ("aminopyrimidine_hinge", "[NH]-c1ncccn1"),
        ("aminoquinazoline_hinge", "[NH]-c1ncnc2ccccc12"),
    ],
    "type_II_dfg_ext": [
        # Extends into DFG-out pocket (allosteric adjacent to ATP site)
        ("diaryl_urea", "c-[NH]-C(=O)-[NH]-c"),
        ("diaryl_amide", "c-[NH]-C(=O)-c"),
    ],
    "type_III_allosteric": [
        # Binds allosteric pocket adjacent to ATP site, no hinge contact
        # Hard to detect by SMARTS alone — flag as "possible" if no hinge hit
    ],
    "covalent_warhead": [
        ("acrylamide", "C=CC(=O)[NH]"),
        ("propynamide", "C#CC(=O)[NH]"),
        ("chloroacetamide", "ClCC(=O)[NH]"),
        ("vinylsulfonamide", "C=CS(=O)(=O)[NH]"),
    ],
}


def classify_kinase_type(smi):
    """Rough classification of kinase inhibitor type based on structural features.

    NOTE: This is heuristic. True classification requires binding mode (X-ray/cryo-EM).
    Use for triage, not as ground truth.

    Returns:
        dict with type indicators and suggested classification
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False}

    hinge_result = classify_hinge_binders(smi)
    has_hinge = hinge_result["has_hinge_binder"]

    # Check covalent warheads
    covalent_hits = []
    for name, smarts in KINASE_TYPE_INDICATORS["covalent_warhead"]:
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            covalent_hits.append(name)

    # Check type II extensions
    type2_hits = []
    for name, smarts in KINASE_TYPE_INDICATORS["type_II_dfg_ext"]:
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            type2_hits.append(name)

    # Classify
    if covalent_hits and has_hinge:
        suggested = "Type VI (covalent, hinge-binding)"
    elif covalent_hits:
        suggested = "Covalent (warhead detected, check binding site)"
    elif has_hinge and type2_hits:
        suggested = "Type II (hinge + DFG-out extension)"
    elif has_hinge:
        suggested = "Type I (hinge-binding, DFG-in)"
    elif type2_hits:
        suggested = "Type II-like (DFG extension, no clear hinge motif)"
    else:
        suggested = "Type III/IV/Unknown (no hinge or DFG motif detected)"

    return {
        "smiles": smi,
        "valid": True,
        "has_hinge": has_hinge,
        "hinge_motifs": [h[0] for h in hinge_result["hits"]],
        "covalent_warheads": covalent_hits,
        "type2_extensions": type2_hits,
        "suggested_type": suggested,
    }
```

## Phase 4: Integration with DMTA Pipeline

### 4.1 Pre-Docking Kinase Filter

Use this BEFORE docking (Phase 9) to avoid wasting compute on non-kinase molecules:

```python
def kinase_pre_filter(smiles_list, require_hinge=True, allow_covalent=True):
    """Filter candidates for kinase plausibility before docking.

    Args:
        smiles_list: list of SMILES
        require_hinge: if True, reject molecules without any hinge binder motif
        allow_covalent: if True, accept covalent warhead molecules

    Returns:
        passed: list of SMILES that pass
        failed: list of (SMILES, reason) tuples
        stats: summary dict
    """
    passed = []
    failed = []

    for smi in smiles_list:
        result = classify_kinase_type(smi)
        if not result["valid"]:
            failed.append((smi, "invalid_smiles"))
            continue

        if require_hinge and not result["has_hinge"]:
            if not (allow_covalent and result["covalent_warheads"]):
                failed.append((smi, f"no_hinge_motif (suggested: {result['suggested_type']})"))
                continue

        passed.append(smi)

    stats = {
        "n_input": len(smiles_list),
        "n_passed": len(passed),
        "n_failed": len(failed),
        "pass_rate": len(passed) / len(smiles_list) * 100 if smiles_list else 0,
    }

    return passed, failed, stats
```

### 4.2 Decision Tree: When to Apply This Skill

```
Is the target a kinase?
├── YES
│   ├── Before generation (skill 3):
│   │   → Run coverage analysis on training set
│   │   → Set coverage KPI: generated pool should have ≥50% of training set coverage
│   │
│   ├── After generation, before ADMET (skill 6):
│   │   → Run kinase_coverage_gate(pool, min_coverage_pct=5.0)
│   │   → If FAIL: stop, fix generator. Do NOT proceed to docking.
│   │   → If PASS: run kinase_pre_filter to remove non-kinase molecules
│   │
│   ├── After docking (skill 9), before final ranking:
│   │   → Classify top candidates by kinase type (I/II/covalent)
│   │   → Report type distribution in final results
│   │
│   └── In cycle summary (skill 12):
│       → Report coverage ratio (generated vs training) as a KPI
│       → If coverage ratio < 0.5: recommend generator improvement in next cycle
│
└── NO
    → Skip this skill entirely. It adds no value for non-kinase targets.
```

## Phase 5: Target-Specific SAR Overlays

### 5.1 ALK5/TGFBR1-Specific Knowledge

```python
ALK5_SAR = {
    "target": "ALK5/TGFBR1",
    "target_class": "Serine/threonine kinase",
    "chembl_id": "CHEMBL260",
    "key_residues": {
        "hinge": ["His283", "Asp351 (DFG)"],
        "gatekeeper": "Thr204",
        "p_loop": "Gly-rich loop (Gly185-Thr189)",
    },
    "known_scaffolds": [
        ("pyrazolopyridine", "c1c[nH]c2ncccc12",
         "GW788388, LY2157299 (galunisertib) class"),
        ("imidazopyridine", "c1cnc2[nH]ccn2c1",
         "SB-431542, SB-525334 class"),
        ("dihydropyrrolopyrazole", None,
         "Vactosertib (EW-7197) class — complex ring system"),
        ("quinazoline_4amine", "[NH]-c1ncnc2ccccc12",
         "Some ALK5 series use quinazoline scaffold"),
    ],
    "selectivity_concerns": [
        "ALK4 (ACVR1B) — very close active site",
        "ALK7 (ACVR1C) — related TGF-β family",
        "p38 MAPK — similar ATP pocket shape",
    ],
    "chronic_disease_flags": [
        "IPF requires long-term dosing → strict safety",
        "Avoid reactive metabolites (N-N bonds, Michael acceptors)",
        "Prefer Type I inhibitors (reversible, ATP-competitive)",
        "hERG liability must be checked (chem-admet skill 6)",
    ],
}
```

### 5.2 Adding New Target SAR

To add SAR for a new kinase target, create a dict with the same structure:

```python
# Template for new kinase target
NEW_TARGET_SAR = {
    "target": "<name>",
    "target_class": "<kinase family>",
    "chembl_id": "<CHEMBLXXXXX>",
    "key_residues": {
        "hinge": ["<residue1>", "<residue2>"],
        "gatekeeper": "<residue>",
    },
    "known_scaffolds": [
        ("<name>", "<SMARTS or None>", "<reference>"),
    ],
    "selectivity_concerns": ["<related kinases>"],
    "chronic_disease_flags": ["<if applicable>"],
}
```

## Failure Modes

1. **SMARTS too broad → false positives.** A bare pyrimidine ring will match molecules that are not kinase inhibitors. Use the specific aminopyrimidine/aminoquinazoline patterns for hard gates; use broad patterns only for soft signals.

2. **SMARTS too narrow → false negatives.** Novel scaffolds won't match any pattern. If coverage is low, check whether the molecules have *other* plausible hinge-binding motifs not in the library. Expand the library when validated.

3. **Confusing coverage with quality.** A molecule can have a hinge binder motif and still be a terrible kinase inhibitor (wrong vector, too bulky, bad pharmacokinetics). This skill answers "is it a plausible kinase inhibitor?", not "is it a good one?" Quality requires docking + ADMET.

4. **Applying to non-kinase targets.** This skill is kinase-specific. Using hinge binder coverage as a KPI for a GPCR or ion channel campaign is meaningless.

5. **Over-relying on Type classification.** The I/II/III/IV classification from SMARTS alone is heuristic. True classification requires structural data. Use it for triage, never as definitive.

## One-Sentence Rule

**If the candidate pool doesn't contain kinase hinge binder motifs at ≥50% of training set frequency, the generator is broken — fix generation before filtering.**
