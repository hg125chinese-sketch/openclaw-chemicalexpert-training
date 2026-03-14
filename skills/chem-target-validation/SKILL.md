---
name: chem-target-validation
description: Validate a therapeutic target for drug discovery using a practical 10-phase framework powered by ToolUniverse APIs (Open Targets, ChEMBL, UniProt, Pharos, DrugBank/clinical proxies, PDB/AlphaFold, PubMed). Input: target (gene symbol or Ensembl), optional disease and modality. Output: validation markdown report + JSON scorecard with GO/NO-GO tiering. Hard gate: target ID disambiguation must succeed.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🎯", "requires": { "bins": ["python3"], "python": ["tooluniverse", "pandas", "requests"] } } }
---

# chem-target-validation — should we work on this target at all?

This skill is for **early target triage**.

Use it when you need to answer questions like:
- “Is this target worth a chemistry campaign?”
- “How strong is the disease link?”
- “Are there already tool compounds / clinical precedents?”
- “Is this a GO / conditional GO / NO-GO target?”

It is built for **practical CE workflows**, not a maximal literature review.

Core rule:
- **Phase 0 must succeed.** If target identity is ambiguous, stop. Do not score garbage.

---

## Workspace variables

```bash
OPENCLAW_WORKSPACE="/home/node/.openclaw/workspace-chemicalexpert"
CHEM_PY="/opt/conda/envs/chem/bin/python"
QMD="/home/node/.openclaw/.npm-global/bin/qmd"
```

---

## Inputs

Required:
- `target`: gene symbol (`TGFBR1`) or Ensembl ID (`ENSG...`)

Optional:
- `disease`: disease name / indication (`idiopathic pulmonary fibrosis`)
- `modality`: small molecule / PROTAC / antibody / RNA / unspecified

---

## Outputs

Write two artifacts:
1. `[TARGET]_[DISEASE]_validation_report.md`
2. `validation_scorecard.json`

Suggested workspace location:

```bash
research/ai4chem/target_validation/<target_slug>/
```

---

## Ten-phase framework (CE practical version)

## Phase 0 — Target disambiguation (**hard gate**)

Goal: confirm that all major databases point to the **same biological target**.

Use ToolUniverse to cross-check:
- Open Targets target record
- ChEMBL target search
- UniProt entry

Minimum success criteria:
- stable primary gene symbol
- stable target name / protein name
- unambiguous Ensembl / UniProt / ChEMBL target mapping

If Phase 0 fails:
- stop
- return `GO/NO-GO = NO-GO`
- report `reason = target_identity_ambiguous`

## Phase 1 — Disease association

Goal: estimate how strongly the target is linked to the disease.

Use:
- Open Targets evidence / association score
- quick literature check (PubMed)

Questions:
- is the association strong and replicated?
- is the directionality coherent? (up/down, gain/loss, activation/inhibition)
- does the literature conflict?

## Phase 2 — Druggability / tractability

Goal: can this target plausibly be modulated by the chosen modality?

Use:
- Open Targets tractability
- Pharos TDL / family context

Interpretation:
- for small molecules, kinase / enzyme / GPCR / ligandable pocket helps
- for other modalities, adapt expectations accordingly

## Phase 3 — Chemical starting points

Goal: determine whether chemistry can start immediately.

Use:
- ChEMBL bioactivity search
- known tool compounds / probes / benchmark ligands

Questions:
- are there sub-µM / low-nM compounds?
- are there known selective tools?
- is this chemistry tractable or still “from zero”?

## Phase 4 — Clinical precedent

Goal: is there translational precedent for the target / pathway?

Use:
- DrugBank-style target-drug relationships available through ToolUniverse
- ClinicalTrials-style target / pathway / intervention search available through ToolUniverse

Questions:
- are there approved / clinical-stage molecules on the same target or pathway?
- were they stopped for efficacy or safety?

## Phase 5 — Safety prefilter

Goal: catch obvious biology-risk before chemistry effort scales.

Use:
- Open Targets safety signals
- tissue expression
- knockout / phenotype evidence

Questions:
- is broad expression likely to create liabilities?
- do KO phenotypes imply severe on-target risk?
- are there obvious contraindication tissues (heart, CNS, immune, liver)?

## Phase 6 — Structural insight

Goal: determine whether structure-based design is realistic.

Use:
- PDB structures
- AlphaFold coverage when no experimental structure exists

Questions:
- do we have a relevant domain / pocket structure?
- is there ligand-bound evidence?
- does this support docking / 3D gen / pocket-conditioned design?

## Phase 7 — Literature deep dive

Goal: move beyond headline counts into contradictions and caveats.

Use:
- PubMed search via ToolUniverse

Look for:
- target role disagreements across disease contexts
- species / model dependence
- pathway compensation / redundancy
- papers that argue *against* the target

## Phase 8 — Integrated score

Score on 0–100:
- **Disease Association**: 0–30
- **Druggability**: 0–25
- **Safety**: 0–20
- **Clinical Precedent**: 0–15
- **Validation Evidence**: 0–10

Tiering:
- **T1**: >75 → **Strong GO**
- **T2**: 50–75 → **Conditional GO**
- **T3**: 25–50 → **Explore**
- **T4**: <25 → **NO-GO**

## Phase 9 — Validation roadmap

Produce a practical next-step plan:
- key biology experiment(s)
- chemical starting point(s)
- structural work
- biomarkers / readouts
- main kill-risk experiment

---

## ToolUniverse integration pattern

All external data access should go through ToolUniverse:

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

result = tu.run({
    "name": "tool_name_here",
    "arguments": {"key": "value"}
})
```

Practical note:
- tool names may vary slightly by ToolUniverse release
- list available tools first, then map them to this framework

Example discovery pattern:

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Adjust to the SDK surface available in your install.
# The key point: discover tool names before hardcoding them.
print(len(tu.tools))
for t in list(tu.tools)[:20]:
    print(t)
```

---

## Minimal working Python example

This is a **template**: adapt exact ToolUniverse tool names to your install.

```python
#!/opt/conda/envs/chem/bin/python
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from tooluniverse import ToolUniverse


def slug(s: str) -> str:
    s = (s or "unknown").strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_") or "unknown"


def safe_run(tu: ToolUniverse, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    return tu.run({"name": name, "arguments": arguments})


def build_scorecard() -> dict[str, Any]:
    return {
        "phase0_ok": False,
        "target": {},
        "disease": {},
        "scores": {
            "disease_association": 0,
            "druggability": 0,
            "safety": 0,
            "clinical_precedent": 0,
            "validation_evidence": 0,
            "total": 0,
        },
        "tier": "T4",
        "decision": "NO-GO",
        "notes": [],
        "roadmap": [],
        "raw": {},
    }


def tier_from_total(total: float) -> tuple[str, str]:
    if total > 75:
        return "T1", "Strong GO"
    if total >= 50:
        return "T2", "Conditional GO"
    if total >= 25:
        return "T3", "Explore"
    return "T4", "NO-GO"


def main(target: str, disease: str | None = None, modality: str | None = None, out_dir: str = "research/ai4chem/target_validation"):
    tu = ToolUniverse()
    tu.load_tools()

    scorecard = build_scorecard()
    scorecard["target"]["query"] = target
    scorecard["disease"]["query"] = disease
    scorecard["target"]["modality"] = modality or "unspecified"

    # ------------------------------
    # Phase 0: disambiguation
    # ------------------------------
    # Replace these tool names with the names exposed by your ToolUniverse install.
    # Typical intent:
    # - Open Targets target lookup
    # - ChEMBL target search
    # - UniProt lookup
    ot_target = safe_run(tu, "open_targets_target_search", {"query": target})
    chembl_target = safe_run(tu, "chembl_target_search", {"query": target})
    uniprot_target = safe_run(tu, "uniprot_search", {"query": target})

    scorecard["raw"]["phase0"] = {
        "open_targets": ot_target,
        "chembl": chembl_target,
        "uniprot": uniprot_target,
    }

    # Minimal pseudo-resolution logic — adapt to returned schemas.
    # You want a stable tuple like: symbol / ensembl / uniprot / chembl_target_id.
    resolved = {
        "symbol": target,
        "ensembl_id": None,
        "uniprot_id": None,
        "chembl_target_id": None,
    }

    # TODO: parse actual tool outputs here.
    # Example placeholder:
    # resolved["ensembl_id"] = ot_target["results"][0]["id"]
    # resolved["uniprot_id"] = uniprot_target["results"][0]["primaryAccession"]
    # resolved["chembl_target_id"] = chembl_target["targets"][0]["target_chembl_id"]

    # Hard gate:
    if not target:
        scorecard["phase0_ok"] = False
        scorecard["notes"].append("Phase 0 failed: target identity ambiguous")
    else:
        scorecard["phase0_ok"] = True
        scorecard["target"].update(resolved)

    if not scorecard["phase0_ok"]:
        total = 0
        tier, decision = tier_from_total(total)
        scorecard["scores"]["total"] = total
        scorecard["tier"] = tier
        scorecard["decision"] = decision
        write_outputs(scorecard, out_dir=out_dir, target=target, disease=disease)
        return

    # ------------------------------
    # Phase 1: disease association
    # ------------------------------
    if disease:
        ot_assoc = safe_run(tu, "open_targets_disease_association", {"target": target, "disease": disease})
        lit_assoc = safe_run(tu, "pubmed_search", {"query": f"{target} AND {disease}"})
        scorecard["raw"]["phase1"] = {"association": ot_assoc, "literature": lit_assoc}
        # TODO: map to 0-30
        scorecard["scores"]["disease_association"] = 20
    else:
        scorecard["notes"].append("No disease provided; disease association score is conservative")
        scorecard["scores"]["disease_association"] = 10

    # ------------------------------
    # Phase 2: druggability
    # ------------------------------
    tract = safe_run(tu, "open_targets_tractability", {"target": target})
    pharos = safe_run(tu, "pharos_target_info", {"query": target})
    scorecard["raw"]["phase2"] = {"tractability": tract, "pharos": pharos}
    # TODO: map to 0-25
    scorecard["scores"]["druggability"] = 15

    # ------------------------------
    # Phase 3: chemical starting points
    # ------------------------------
    chembl_act = safe_run(tu, "chembl_activity_search", {"target": target})
    scorecard["raw"]["phase3"] = {"chembl_bioactivity": chembl_act}
    scorecard["scores"]["validation_evidence"] = 6

    # ------------------------------
    # Phase 4: clinical precedent
    # ------------------------------
    drugbank = safe_run(tu, "drugbank_target_lookup", {"target": target})
    trials = safe_run(tu, "clinicaltrials_search", {"query": f"{target} {disease or ''}".strip()})
    scorecard["raw"]["phase4"] = {"drugbank": drugbank, "clinical_trials": trials}
    scorecard["scores"]["clinical_precedent"] = 8

    # ------------------------------
    # Phase 5: safety
    # ------------------------------
    ot_safety = safe_run(tu, "open_targets_safety", {"target": target})
    expression = safe_run(tu, "open_targets_expression", {"target": target})
    ko = safe_run(tu, "open_targets_phenotypes", {"target": target})
    scorecard["raw"]["phase5"] = {"safety": ot_safety, "expression": expression, "ko": ko}
    scorecard["scores"]["safety"] = 12

    # ------------------------------
    # Phase 6: structural insight
    # ------------------------------
    pdb_hits = safe_run(tu, "pdb_search", {"query": target})
    af = safe_run(tu, "alphafold_lookup", {"target": target})
    scorecard["raw"]["phase6"] = {"pdb": pdb_hits, "alphafold": af}

    # ------------------------------
    # Phase 7: literature deep dive
    # ------------------------------
    deep_lit = safe_run(tu, "pubmed_search", {"query": f"{target} mechanism resistance safety biomarker"})
    conflict_lit = safe_run(tu, "pubmed_search", {"query": f"{target} contradictory OR conflicting OR controversy"})
    scorecard["raw"]["phase7"] = {"deep_literature": deep_lit, "conflicts": conflict_lit}

    # ------------------------------
    # Phase 8: integrated score
    # ------------------------------
    total = sum(scorecard["scores"][k] for k in [
        "disease_association",
        "druggability",
        "safety",
        "clinical_precedent",
        "validation_evidence",
    ])
    scorecard["scores"]["total"] = total
    tier, decision = tier_from_total(total)
    scorecard["tier"] = tier
    scorecard["decision"] = decision

    # ------------------------------
    # Phase 9: validation roadmap
    # ------------------------------
    scorecard["roadmap"] = [
        "Confirm disease-directional biology in the most relevant model.",
        "Nominate at least one tractable assay and one kill-risk assay.",
        "List known tool compounds / starting ligands from ChEMBL.",
        "Define a biomarker or target-engagement readout before chemistry scales.",
    ]

    write_outputs(scorecard, out_dir=out_dir, target=target, disease=disease)


def write_outputs(scorecard: dict[str, Any], *, out_dir: str, target: str, disease: str | None):
    tslug = slug(target)
    dslug = slug(disease or "general")
    root = Path(out_dir) / tslug
    root.mkdir(parents=True, exist_ok=True)

    report_md = root / f"{tslug}_{dslug}_validation_report.md"
    score_json = root / "validation_scorecard.json"

    score_json.write_text(json.dumps(scorecard, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"# Target validation — {target}")
    if disease:
        lines.append(f" / {disease}")
    lines.append("\n\n")

    lines.append("## Decision\n\n")
    lines.append(f"- Tier: **{scorecard['tier']}**\n")
    lines.append(f"- Recommendation: **{scorecard['decision']}**\n")
    lines.append(f"- Total score: **{scorecard['scores']['total']}/100**\n\n")

    lines.append("## Score breakdown\n\n")
    lines.append(f"- Disease Association: {scorecard['scores']['disease_association']}/30\n")
    lines.append(f"- Druggability: {scorecard['scores']['druggability']}/25\n")
    lines.append(f"- Safety: {scorecard['scores']['safety']}/20\n")
    lines.append(f"- Clinical Precedent: {scorecard['scores']['clinical_precedent']}/15\n")
    lines.append(f"- Validation Evidence: {scorecard['scores']['validation_evidence']}/10\n\n")

    lines.append("## Notes\n\n")
    if scorecard['notes']:
        for n in scorecard['notes']:
            lines.append(f"- {n}\n")
    else:
        lines.append("- No major blocking notes captured in template run.\n")

    lines.append("\n## Validation roadmap\n\n")
    for i, step in enumerate(scorecard['roadmap'], start=1):
        lines.append(f"{i}. {step}\n")

    report_md.write_text("".join(lines), encoding="utf-8")
    print(report_md)
    print(score_json)
```

---

## Practical scoring guidance

Use the score as a **decision aid**, not fake precision.

### Disease Association (0–30)
- 25–30: strong genetics / multiple orthogonal evidence streams / disease-relevant direction clear
- 15–24: decent association, but some gaps or inconsistency
- 5–14: weak or indirect
- 0–4: speculative

### Druggability (0–25)
- 20–25: highly tractable target class, good modality fit
- 10–19: plausible but not easy
- 0–9: hard / unclear tractability

### Safety (0–20)
- 15–20: manageable known biology risk
- 8–14: moderate concern
- 0–7: severe on-target biology risk

### Clinical Precedent (0–15)
- 12–15: strong same-target or same-pathway precedent
- 6–11: partial precedent
- 0–5: little/no precedent

### Validation Evidence (0–10)
- 8–10: strong chemistry + literature + mechanistic support
- 4–7: partial support
- 0–3: little direct support

---

## Reporting checklist

- [ ] Phase 0 disambiguation explicitly documented
- [ ] disease link explained, not just score copied
- [ ] tractability tied to modality
- [ ] known compounds / tool molecules listed when available
- [ ] safety caveats written plainly
- [ ] structural evidence summarized (PDB / AlphaFold / pocket status)
- [ ] conflict literature noted
- [ ] scorecard JSON written
- [ ] final GO / NO-GO rationale auditable

---

## Style rule

Be a **lab-mate**, not a brochure.

Good:
- “Disease link is decent, but directionality is noisy across tissues.”
- “Chemistry is not starting from zero; ChEMBL suggests tool-like matter exists.”
- “Safety risk is probably the main kill factor, not tractability.”

Bad:
- “This target is highly promising and transformative.”
- vague hype without kill-risk analysis
