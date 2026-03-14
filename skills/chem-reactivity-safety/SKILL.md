---
name: chem-reactivity-safety
description: Systematic structural alert screening for drug candidates. Covers PAINS, reactive metabolites, genotoxicity alerts, metabolic liabilities, and chronic disease safety flags. Outputs auditable hit reports with SMARTS evidence. Use as a hard gate before any candidate enters final ranking.
homepage: https://www.rdkit.org/docs/source/rdkit.Chem.FilterCatalog.html
metadata: { "openclaw": { "emoji": "🚨", "requires": { "bins": ["python3"], "python": ["rdkit", "numpy", "pandas"] } } }
---

# Reactivity & Safety: Structural Alert Screening

Drug candidates can look great on paper — good activity, decent QED, passing Lipinski — and still fail because they contain structural features that cause toxicity, metabolic instability, or assay interference. This skill turns medicinal chemistry safety knowledge into computable, auditable gates.

## When to Use

- After molecular generation (skill 3) and ADMET filtering (skill 6), before final ranking
- When evaluating candidates for **chronic disease** targets (IPF, diabetes, autoimmune) — safety bar is higher
- When a candidate contains unusual functional groups or unfamiliar motifs
- When building or updating a project-specific denylist
- As part of the chem-experiment (skill 12) DMTA cycle gate system

## Core Philosophy

1. **Structural alerts are not predictions — they are red flags.** A hit doesn't mean the molecule is toxic. It means the molecule needs extra scrutiny before advancing. But in early-stage virtual screening, extra scrutiny means "remove it" unless there's a compelling reason to keep it.
2. **Layer your filters: universal → context-specific.** PAINS and reactive metabolite alerts apply to everything. Chronic disease flags (metabolic stability, covalent liability) are context-dependent. Don't apply chronic disease gates to an oncology program where covalent inhibitors are intentional.
3. **Every flag must be auditable.** The output is not just "failed" — it's "failed rule X because SMARTS Y matched atoms Z." The chemist (or agent) must be able to inspect and override with justification.
4. **This skill complements ADMET (skill 6), not replaces it.** ADMET covers drug-likeness properties (LogP, TPSA, QED). This skill covers structural liabilities that ADMET doesn't catch — reactive groups, toxicophores, metabolic soft spots.
5. **When in doubt, flag it.** False positives are cheap (you just review them). False negatives are expensive (a toxic compound advances).

## Phase 1: Universal Structural Alerts (Always Apply)

### 1.1 PAINS (Pan-Assay Interference Compounds)

PAINS are molecules that give false positives in many assay types. They're not necessarily toxic, but they waste experimental resources.

```python
#!/opt/conda/envs/chem/bin/python
"""PAINS filter using RDKit's built-in FilterCatalog."""
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def get_pains_catalog():
    """Load RDKit's built-in PAINS filter catalog."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


PAINS_CATALOG = get_pains_catalog()


def check_pains(smi):
    """Check a molecule against PAINS filters.

    Returns:
        dict with keys:
        - smiles: input
        - valid: bool
        - is_pains: bool
        - hits: list of (filter_name, description)
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "is_pains": False, "hits": []}

    entries = PAINS_CATALOG.GetMatches(mol)
    hits = [(e.GetDescription(), str(e.GetFilterMatches(mol)))
            for e in entries]

    return {
        "smiles": smi,
        "valid": True,
        "is_pains": len(hits) > 0,
        "hits": hits,
    }
```

### 1.2 Reactive Metabolite Alerts

These motifs can form reactive intermediates (epoxides, quinones, nitroso compounds) that covalently modify proteins, causing idiosyncratic toxicity.

```python
REACTIVE_METABOLITE_ALERTS = [
    # --- Nitrogen-based ---
    ("aniline_unsubst", "[NH2]-c1ccccc1",
     "Unsubstituted aniline → hydroxylamine → nitroso metabolite. "
     "Risk: idiosyncratic hepatotoxicity."),
    ("hydrazine", "[NX3H2]-[NX3H]",
     "Hydrazine → diazene metabolites. "
     "Risk: hepatotoxicity, genotoxicity."),
    ("hydrazone", "[CX3]=[NX2]-[NX3]",
     "Hydrazone → hydrolysis + reactive aldehyde/ketone. "
     "Risk: metabolic instability, reactive intermediates."),
    ("n_n_single_bond", "[NX3,NX2;!$([N-])]-[NX3,NX2;!$([N-])]",
     "N-N single bond → oxidative cleavage risk. "
     "Risk: metabolic instability (especially for chronic dosing)."),
    ("nitro_aromatic", "[N+](=O)[O-]-c",
     "Aromatic nitro → nitroso → hydroxylamine cascade. "
     "Risk: genotoxicity, methemoglobinemia."),
    ("nitroso", "[NX2]=O",
     "Nitroso group — directly reactive with proteins/DNA."),
    ("hydroxamic_acid", "[OH]-[NX3]-C(=O)",
     "Hydroxamic acid → reactive nitroso metabolite. "
     "Note: some approved drugs (e.g., HDAC inhibitors) use this deliberately."),

    # --- Sulfur-based ---
    ("thiol_free", "[SH]",
     "Free thiol → oxidation, disulfide formation, metal chelation. "
     "Risk: assay interference, metabolic instability."),
    ("thiourea", "[NX3]-C(=S)-[NX3]",
     "Thiourea → oxidative desulfurization → reactive sulfoxide. "
     "Risk: thyroid toxicity, hepatotoxicity."),
    ("thioamide", "C(=S)-[NX3]",
     "Thioamide → similar to thiourea liability."),

    # --- Carbon electrophiles ---
    ("michael_acceptor", "[CX3]=[CX3]-[CX3](=O)",
     "α,β-unsaturated carbonyl (Michael acceptor). "
     "Risk: covalent protein modification. "
     "Note: intentional in some covalent drugs (ibrutinib acrylamide)."),
    ("acyl_halide", "C(=O)-[F,Cl,Br,I]",
     "Acyl halide — highly reactive electrophile."),
    ("aldehyde", "[CX3H1](=O)[#6]",
     "Aldehyde → Schiff base with lysine, metabolic oxidation. "
     "Risk: mutagenicity, protein adducts."),
    ("epoxide", "C1OC1",
     "Epoxide — alkylating agent. "
     "Risk: genotoxicity (unless rapidly detoxified)."),
    ("aziridine", "C1NC1",
     "Aziridine — alkylating agent like epoxide."),

    # --- Quinone/redox ---
    ("quinone", "O=C1C=CC(=O)C=C1",
     "Quinone — redox cycling, protein adducts. "
     "Risk: oxidative stress, hepatotoxicity."),
    ("hydroquinone", "Oc1ccc(O)cc1",
     "Hydroquinone → quinone via oxidation. "
     "Risk: if unprotected, will form quinone in vivo."),
    ("catechol", "Oc1ccccc1O",
     "Catechol → ortho-quinone via COMT/oxidation. "
     "Risk: reactive metabolite, COMT inhibition."),

    # --- Halide-based ---
    ("alkyl_halide", "[CX4]-[Cl,Br,I]",
     "Alkyl halide — SN2 electrophile. "
     "Risk: DNA alkylation, mutagenicity. "
     "Note: -F is generally safe (C-F bond too strong)."),
    ("benzyl_halide", "c-[CH2]-[Cl,Br,I]",
     "Benzylic halide — activated alkyl halide. Higher reactivity."),

    # --- Miscellaneous ---
    ("acyl_hydrazide", "C(=O)-[NH]-[NH2]",
     "Acyl hydrazide → hydrolysis + hydrazine release."),
    ("peroxide", "OO",
     "Peroxide — homolytic cleavage risk."),
    ("azide_organic", "[N-]=[N+]=N-[#6]",
     "Organic azide — potential explosive, nitrene precursor."),
]


def compile_alerts(alert_list):
    """Compile SMARTS patterns. Returns list of (name, pattern, description)."""
    compiled = []
    for name, smarts, desc in alert_list:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            print(f"WARNING: Failed to compile SMARTS for {name}: {smarts}")
            continue
        compiled.append((name, pat, desc))
    return compiled


COMPILED_REACTIVE = compile_alerts(REACTIVE_METABOLITE_ALERTS)
```

### 1.3 Check a Single Molecule

```python
def check_reactive_alerts(smi, alerts=None):
    """Screen a molecule against reactive metabolite alerts.

    Returns:
        dict with:
        - smiles: input
        - valid: bool
        - n_alerts: int
        - alerts: list of (name, description, matched_atoms)
        - severity: "clean" | "warning" | "danger"
    """
    if alerts is None:
        alerts = COMPILED_REACTIVE

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "n_alerts": 0,
                "alerts": [], "severity": "invalid"}

    hits = []
    for name, pat, desc in alerts:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            hits.append({
                "name": name,
                "description": desc,
                "n_matches": len(matches),
                "atom_indices": [list(m) for m in matches],
            })

    n = len(hits)
    if n == 0:
        severity = "clean"
    elif n <= 2:
        severity = "warning"
    else:
        severity = "danger"

    return {
        "smiles": smi,
        "valid": True,
        "n_alerts": n,
        "alerts": hits,
        "severity": severity,
    }
```

## Phase 2: Genotoxicity / Mutagenicity Alerts

### 2.1 Beilstein/Lhasa Alerts (Simplified)

These patterns are associated with Ames-positive (mutagenic) outcomes:

```python
GENOTOXICITY_ALERTS = [
    ("aromatic_nitro", "[N+](=O)[O-]-c",
     "Aromatic nitro — strong Ames predictor. "
     "Metabolic activation to hydroxylamine → nitroso → nitrenium."),
    ("aromatic_azo", "c-N=N-c",
     "Azo compound — reductive cleavage to aromatic amines."),
    ("aromatic_amine_primary", "[NH2]-c1ccccc1",
     "Primary aromatic amine — N-oxidation → nitrenium ion. "
     "The #1 structural alert for Ames mutagenicity."),
    ("n_mustard", "[NX3](-[CH2]-[Cl,Br])(-[CH2]-[Cl,Br])",
     "Nitrogen mustard — bis-alkylating agent. "
     "Intentional in some chemotherapeutics, never in chronic disease drugs."),
    ("epoxide_genotox", "C1OC1",
     "Epoxide — direct DNA alkylation."),
    ("alkyl_nitrosamine", "[NX3]-N=O",
     "N-nitrosamine — potent mutagen. "
     "ppb-level concern (see FDA NDSRI guidance)."),
    ("polycyclic_flat", "c1ccc2c(c1)ccc1ccccc12",
     "Tricyclic aromatic — potential intercalator / DNA adduct. "
     "Flag for review, not automatic rejection."),
]

COMPILED_GENOTOX = compile_alerts(GENOTOXICITY_ALERTS)


def check_genotoxicity(smi):
    """Screen for genotoxicity/mutagenicity structural alerts.

    Returns same format as check_reactive_alerts.
    Severity thresholds are stricter: any hit = "danger" for chronic disease.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "n_alerts": 0,
                "alerts": [], "severity": "invalid"}

    hits = []
    for name, pat, desc in COMPILED_GENOTOX:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            hits.append({
                "name": name,
                "description": desc,
                "n_matches": len(matches),
                "atom_indices": [list(m) for m in matches],
            })

    n = len(hits)
    # Genotoxicity is more serious — any hit is at least "warning"
    if n == 0:
        severity = "clean"
    elif n == 1:
        severity = "warning"
    else:
        severity = "danger"

    return {
        "smiles": smi,
        "valid": True,
        "n_alerts": n,
        "alerts": hits,
        "severity": severity,
    }
```

## Phase 3: Chronic Disease Safety Layer

### 3.1 Context-Dependent Flags

These alerts are specifically relevant for **chronic disease drug candidates** (IPF, diabetes, autoimmune, cardiovascular) where patients take the drug for months to years.

```python
CHRONIC_DISEASE_FLAGS = [
    # Metabolic instability (long-term exposure = more metabolite formation)
    ("n_n_bond_chronic", "[NX3,NX2;!$([N-])]-[NX3,NX2;!$([N-])]",
     "N-N bond in chronic context → cumulative metabolic instability risk."),
    ("hydrazone_chronic", "[CX3]=[NX2]-[NX3]",
     "Hydrazone in chronic context → hydrolysis over repeated dosing."),
    ("ester_prodrug_risk", "C(=O)-O-[CX4]",
     "Alkyl ester → rapid hydrolysis by esterases. "
     "May be intentional (prodrug) but flag for chronic dosing."),

    # Covalent liability (acceptable in oncology, not in chronic disease)
    ("acrylamide_chronic", "C=CC(=O)[NH]",
     "Acrylamide warhead — covalent modifier. "
     "Acceptable in oncology (ibrutinib), problematic for chronic dosing."),
    ("vinyl_sulfone_chronic", "C=CS(=O)(=O)",
     "Vinyl sulfone — covalent warhead. Not for chronic use."),
    ("chloroacetamide_chronic", "ClCC(=O)[NH]",
     "Chloroacetamide — covalent warhead. Not for chronic use."),

    # Photoactivation risk (long-term outdoor exposure)
    ("quinolone_phototox", "O=c1cc(-c2ccccc2)oc2ccccc12",
     "Quinolone/fluoroquinolone-like → phototoxicity risk. "
     "Relevant for long-term use with sun exposure."),

    # Accumulation risk
    ("high_logp_chronic", None,  # Computed, not SMARTS
     "LogP > 5 with chronic dosing → tissue accumulation risk. "
     "Already partially caught by Lipinski, but flag explicitly here."),
]

# Only compile SMARTS-based alerts (skip computed ones like LogP)
COMPILED_CHRONIC = compile_alerts(
    [(n, s, d) for n, s, d in CHRONIC_DISEASE_FLAGS if s is not None]
)


def check_chronic_safety(smi, is_chronic=True):
    """Screen for chronic disease-specific safety concerns.

    Args:
        smi: SMILES string
        is_chronic: whether this is a chronic disease target
                    If False, this function returns clean for everything.

    Returns:
        same format as check_reactive_alerts
    """
    if not is_chronic:
        return {"smiles": smi, "valid": True, "n_alerts": 0,
                "alerts": [], "severity": "clean",
                "note": "Chronic safety layer skipped (not a chronic target)"}

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "n_alerts": 0,
                "alerts": [], "severity": "invalid"}

    hits = []
    for name, pat, desc in COMPILED_CHRONIC:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            hits.append({
                "name": name,
                "description": desc,
                "n_matches": len(matches),
                "atom_indices": [list(m) for m in matches],
            })

    # Also check LogP for chronic accumulation risk
    from rdkit.Chem import Descriptors
    logp = Descriptors.MolLogP(mol)
    if logp > 5.0:
        hits.append({
            "name": "high_logp_chronic",
            "description": f"LogP = {logp:.2f} > 5.0 → tissue accumulation risk "
                           "with chronic dosing.",
            "n_matches": 1,
            "atom_indices": [],
        })

    n = len(hits)
    if n == 0:
        severity = "clean"
    elif n <= 1:
        severity = "warning"
    else:
        severity = "danger"

    return {
        "smiles": smi,
        "valid": True,
        "n_alerts": n,
        "alerts": hits,
        "severity": severity,
    }
```

## Phase 4: Unified Safety Screening Pipeline

### 4.1 Full Safety Report

```python
def full_safety_screen(smi, is_chronic=True):
    """Run all safety screens and produce a unified report.

    Args:
        smi: SMILES string
        is_chronic: whether chronic disease flags should be applied

    Returns:
        dict with:
        - smiles: input
        - valid: bool
        - overall_severity: "clean" | "warning" | "danger"
        - pass_hard_gate: bool (True if overall_severity != "danger")
        - layers: dict of individual screen results
        - all_alerts: flat list of all alerts with layer attribution
        - summary: human-readable one-line summary
    """
    pains = check_pains(smi)
    reactive = check_reactive_alerts(smi)
    genotox = check_genotoxicity(smi)
    chronic = check_chronic_safety(smi, is_chronic=is_chronic)

    if not all(r.get("valid", False) for r in [pains, reactive, genotox, chronic]):
        return {"smiles": smi, "valid": False, "overall_severity": "invalid",
                "pass_hard_gate": False, "layers": {}, "all_alerts": [],
                "summary": "Invalid SMILES"}

    # Collect all alerts with layer attribution
    all_alerts = []
    for layer_name, result in [("pains", pains), ("reactive", reactive),
                                ("genotox", genotox), ("chronic", chronic)]:
        if layer_name == "pains":
            for hit_name, hit_detail in result.get("hits", []):
                all_alerts.append({"layer": layer_name, "name": hit_name,
                                   "detail": hit_detail})
        else:
            for alert in result.get("alerts", []):
                all_alerts.append({"layer": layer_name, **alert})

    # Overall severity: worst of all layers
    severities = [reactive["severity"], genotox["severity"], chronic["severity"]]
    if pains["is_pains"]:
        severities.append("warning")

    if "danger" in severities:
        overall = "danger"
    elif "warning" in severities:
        overall = "warning"
    else:
        overall = "clean"

    n_total = len(all_alerts)
    summary = (f"{n_total} alert(s): {overall}" if n_total > 0
               else "Clean — no structural alerts detected")

    return {
        "smiles": smi,
        "valid": True,
        "overall_severity": overall,
        "pass_hard_gate": overall != "danger",
        "n_total_alerts": n_total,
        "layers": {
            "pains": pains,
            "reactive": reactive,
            "genotox": genotox,
            "chronic": chronic,
        },
        "all_alerts": all_alerts,
        "summary": summary,
    }
```

### 4.2 Batch Screening with Statistics

```python
import pandas as pd


def batch_safety_screen(smiles_list, is_chronic=True, label="pool"):
    """Screen a pool of molecules and produce summary statistics.

    Returns:
        dict with:
        - summary: pool-level stats
        - per_molecule: list of full_safety_screen results
        - alert_frequency: DataFrame of which alerts fire most often
    """
    results = [full_safety_screen(smi, is_chronic) for smi in smiles_list]

    n_total = len(results)
    n_valid = sum(1 for r in results if r["valid"])
    n_clean = sum(1 for r in results if r["overall_severity"] == "clean")
    n_warning = sum(1 for r in results if r["overall_severity"] == "warning")
    n_danger = sum(1 for r in results if r["overall_severity"] == "danger")
    n_pass = sum(1 for r in results if r["pass_hard_gate"])

    # Alert frequency
    alert_counts = {}
    for r in results:
        for a in r.get("all_alerts", []):
            key = f"{a['layer']}:{a['name']}"
            alert_counts[key] = alert_counts.get(key, 0) + 1

    alert_freq = pd.DataFrame([
        {"alert": k, "count": v, "pct": 100.0 * v / n_valid if n_valid > 0 else 0}
        for k, v in alert_counts.items()
    ]).sort_values("count", ascending=False) if alert_counts else pd.DataFrame()

    summary = {
        "label": label,
        "n_total": n_total,
        "n_valid": n_valid,
        "n_clean": n_clean,
        "n_warning": n_warning,
        "n_danger": n_danger,
        "n_pass_hard_gate": n_pass,
        "pass_rate": 100.0 * n_pass / n_valid if n_valid > 0 else 0,
        "top3_alerts": alert_freq.head(3).to_dict("records") if len(alert_freq) > 0 else [],
    }

    return {
        "summary": summary,
        "per_molecule": results,
        "alert_frequency": alert_freq,
    }
```

## Phase 5: Integration with DMTA Pipeline

### 5.1 Decision Tree

```
After ADMET gates (skill 6), before docking (skill 9):
│
├── Run full_safety_screen(smi, is_chronic=<target_context>)
│
├── overall_severity == "danger"?
│   ├── YES → REJECT. Do not dock. Log reason.
│   └── NO → Continue to docking.
│
├── overall_severity == "warning"?
│   ├── Add to "flagged" list with alert details
│   ├── Continue to docking, but tag in final report
│   └── In Top 5 output: explicitly note warnings for chemist review
│
└── In cycle summary (skill 12):
    ├── Report: N molecules rejected by safety screen
    ├── Report: top 3 most frequent alerts in pool
    └── If >20% rejected → flag generation quality issue
```

### 5.2 Project-Specific Denylist Extension

Every project should maintain its own denylist on top of the universal alerts:

```python
def make_project_denylist(additional_alerts):
    """Create a project-specific denylist by extending universal alerts.

    Args:
        additional_alerts: list of (name, SMARTS, description) tuples

    Returns:
        compiled alert list that includes both universal and project-specific alerts

    Example:
        ipf_extras = [
            ("ipf_n_oxide_risk", "[NX3](=O)", "N-oxide risk for IPF chronic dosing"),
        ]
        ipf_alerts = make_project_denylist(ipf_extras)
    """
    combined = REACTIVE_METABOLITE_ALERTS + additional_alerts
    return compile_alerts(combined)
```

## Failure Modes

1. **SMARTS too broad → excessive false positives.** For example, a generic `[NH2]` pattern matches amino acids and simple amines that are perfectly safe. Always use the most specific SMARTS that captures the actual liability (e.g., `[NH2]-c` for aromatic amines, not `[NH2]`).

2. **Missing context → wrong severity.** A Michael acceptor acrylamide is "danger" for an IPF chronic drug but "by design" for a covalent kinase inhibitor. Always pass `is_chronic` correctly. When unsure, flag and let the chemist decide.

3. **Confusing structural alerts with toxicity.** An alert means "investigate further," not "this molecule is toxic." Many approved drugs contain structural alerts (metformin has an N-N bond, vorinostat has a hydroxamic acid). The difference is that those went through extensive safety testing.

4. **Alert fatigue from noisy screens.** If >50% of your pool gets flagged, the alerts are too broad or your generator is producing low-quality molecules. Check which specific alerts fire most and tighten or remove the noisy ones.

5. **Not updating project denylist.** If you discover a new problematic motif during a project (like the N-N bond in IPF Cycle 1), add it to the project denylist immediately. Don't rely on memory.

## Relationship to Other Skills

| Skill | Relationship |
|-------|-------------|
| chem-admet (6) | ADMET covers drug-likeness properties; this skill covers structural liabilities. Run ADMET first, then safety screen. |
| chem-kinase-sar (13) | Kinase SAR checks "does it look like a kinase inhibitor?"; this skill checks "is it safe?" Both are pre-docking gates. |
| chem-molgen (3) | If safety screen rejects >20% of generated pool, the generator is producing too many problematic motifs → generation quality issue. |
| chem-experiment (12) | Safety screen is a hard gate in the DMTA cycle. "danger" = reject, "warning" = flag. |

## One-Sentence Rule

**Every candidate must pass structural alert screening before entering final ranking — "danger" is a hard reject, "warning" is a flag for chemist review.**

## Phase 2: ADMET-AI ML Enhancement

Keep the existing SMARTS workflow exactly as-is.

Use ADMET-AI as a **second layer**, not a replacement:
- **Layer 1 (hard filter):** SMARTS / PAINS / reactive / chronic alerts
- **Layer 2 (ML enhancement):** ADMET-AI predictions via ToolUniverse

### Policy

1. **SMARTS remains the hard gate.**
   - If the SMARTS layer says hard reject, reject.
2. **ADMET-AI is a review layer.**
   - If SMARTS says PASS but ADMET-AI says high risk, do **not** auto-reject.
   - Mark as **`flag_for_review`**.
3. **Backward compatibility is required.**
   - Existing SMARTS logic and outputs remain valid.
   - This phase only adds `admet_scores` and `combined_risk`.

### What to predict

Use ToolUniverse ADMET-AI tools when available for signals such as:
- BBB penetration
- CYP inhibition / metabolic liability
- oral bioavailability
- hERG risk
- (optionally) hepatotoxicity / clearance / plasma protein binding

### Minimal ToolUniverse pattern

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

result = tu.run({
    "name": "<ADMET_AI_tool_name>",
    "arguments": {"smiles": "CCO"}
})
```

### Example code (ML enhancement wrapper)

```python
#!/opt/conda/envs/chem/bin/python
from __future__ import annotations

from typing import Any
from tooluniverse import ToolUniverse


def _safe_run(tu: ToolUniverse, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    try:
        return tu.run({"name": name, "arguments": arguments})
    except Exception as e:
        return {"status": "error", "error": str(e), "tool": name}


def admet_ai_layer(smiles: str) -> dict[str, Any]:
    """Second-layer ML enhancement using ToolUniverse ADMET-AI tools.

    Adjust tool names to match your ToolUniverse install.
    Returns a normalized score bundle.
    """
    tu = ToolUniverse()
    tu.load_tools()

    # Replace these with the exact ADMET-AI tool names exposed in your install.
    bbb = _safe_run(tu, "ADMETAI_predict_bbb", {"smiles": smiles})
    cyp = _safe_run(tu, "ADMETAI_predict_cyp_inhibition", {"smiles": smiles})
    oral = _safe_run(tu, "ADMETAI_predict_oral_bioavailability", {"smiles": smiles})
    herg = _safe_run(tu, "ADMETAI_predict_herg", {"smiles": smiles})

    # Normalize into a stable dictionary even if raw tool schemas differ.
    return {
        "bbb": bbb,
        "cyp": cyp,
        "oral_bioavailability": oral,
        "herg": herg,
    }


def combine_structural_and_ml(*, structural_result: dict[str, Any], admet_scores: dict[str, Any]) -> dict[str, Any]:
    """Combine rule-based safety with ADMET-AI predictions.

    combined_risk semantics:
    - high   : SMARTS layer already says danger / hard reject
    - medium : SMARTS passes but ML flags high risk -> review
    - low    : no strong rule hits, no strong ML concern detected
    """
    severity = structural_result.get("overall_severity", "clean")
    hard_reject = structural_result.get("hard_reject", False)

    # Placeholder ML interpretation logic.
    # Replace with schema-aware parsing for your ToolUniverse ADMET-AI outputs.
    ml_high_risk = False
    review_reasons = []

    text_blob = str(admet_scores).lower()
    if "herg" in text_blob and ("high" in text_blob or "risk" in text_blob):
        ml_high_risk = True
        review_reasons.append("ml:herg")
    if "bbb" in text_blob and "high" in text_blob:
        review_reasons.append("ml:bbb")
    if "cyp" in text_blob and "inhib" in text_blob:
        review_reasons.append("ml:cyp")

    flag_for_review = False
    if hard_reject or severity == "danger":
        combined_risk = "high"
    elif ml_high_risk:
        combined_risk = "medium"
        flag_for_review = True
    else:
        combined_risk = "low"

    return {
        "alerts": structural_result.get("alerts", []),
        "admet_scores": admet_scores,
        "combined_risk": combined_risk,
        "flag_for_review": flag_for_review,
        "review_reasons": review_reasons,
    }
```

### Upgraded output format

For each molecule, emit:

```python
{
  "smiles": "CCOc1ccc(...)",
  "alerts": [...],
  "admet_scores": {
    "bbb": {...},
    "cyp": {...},
    "oral_bioavailability": {...},
    "herg": {...}
  },
  "combined_risk": "low"
}
```

### Conflict handling

The critical case is:
- **SMARTS = PASS**
- **ADMET-AI = high risk**

Action:
- do **not** auto reject
- return `combined_risk = "medium"`
- return `flag_for_review = True`
- keep the molecule visible in downstream reports

Why:
- SMARTS are explicit mechanistic red flags
- ADMET-AI is probabilistic and should guide triage, not silently override hard chemistry logic
