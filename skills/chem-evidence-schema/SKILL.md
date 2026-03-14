---
name: chem-evidence-schema
description: Standardize chemical biology and DMTA evidence into a common Evidence / EvidenceCollection schema for cross-skill interoperability, conflict detection, ranking, and auditability. Use when normalizing outputs from docking, safety, PLIF, PoseBusters, Boltz-2, target validation, or panel selection into a shared JSON structure.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧾", "requires": { "bins": ["python3"], "python": ["pydantic", "pandas"] } } }
---

# chem-evidence-schema — make evidence composable

If every skill emits its own private format, the pipeline becomes glue-code soup.

This skill defines a **shared evidence object model** so outputs from docking, safety, PLIF, Boltz‑2, PoseBusters, and target validation can live in one container and be compared, filtered, and audited.

Use it when you need to:
- merge signals from multiple skills
- track provenance (`source_tool`, `source_db`)
- assign evidence grades
- detect conflicts across tools
- hand one normalized bundle into **chem-panel-selection** or downstream reports

Core rule:
- **Do not overwrite original outputs.** Wrap them into Evidence objects.

---

## Output objects

### Evidence

Minimum fields:
- `evidence_id`
- `source_tool`
- `source_db`
- `evidence_type`
- `evidence_grade`
- `target_id`
- `mol_id` (optional)
- `value`
- `unit`
- `confidence`
- `timestamp`
- `metadata`
- `conflicts`

### EvidenceCollection

A list-like container of Evidence with helper methods:
- `filter_by(...)`
- `summarize()`
- `detect_conflicts()`
- `to_json()` / `from_json()`

---

## Evidence JSON schema (practical)

```json
{
  "evidence_id": "cycle6_mol_0064_vina_1vjy",
  "source_tool": "vina_docking",
  "source_db": "PDB",
  "evidence_type": "binding",
  "evidence_grade": "T3",
  "target_id": "ENSG00000106799",
  "mol_id": "mol_0064",
  "value": -10.04,
  "unit": "kcal/mol",
  "confidence": 0.45,
  "timestamp": "2026-03-14T12:00:00Z",
  "metadata": {
    "pdb_id": "1VJY",
    "box_source": "redock_best",
    "seed": 0
  },
  "conflicts": []
}
```

---

## Grade semantics

- **T1** — direct experimental evidence
  - e.g. measured IC50, crystallography, in vitro assay, clinical evidence
- **T2** — computation with strong literature / reference support
  - e.g. hinge H-bond vs co-crystal, PoseBusters valid, curated safety rule hit
- **T3** — computation / model prediction
  - e.g. Vina score, MACE strain, Boltz-2 binder probability, ADMET model output
- **T4** — inferred / indirect / weakly supported
  - e.g. qualitative synthesis note, indirect pathway extrapolation

Practical rule:
- do not inflate grades just because a number looks precise

---

## Evidence types

Use a controlled set unless you have a strong reason not to:
- `binding`
- `safety`
- `disease_association`
- `structural`
- `literature`
- `clinical`
- `expression`

You may add others in `metadata`, but keep `evidence_type` stable for aggregation.

---

## Pydantic reference implementation

```python
#!/opt/conda/envs/chem/bin/python
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from typing import Any, Iterable

from pydantic import BaseModel, Field


VALID_TYPES = {
    "binding",
    "safety",
    "disease_association",
    "structural",
    "literature",
    "clinical",
    "expression",
}
VALID_GRADES = {"T1", "T2", "T3", "T4"}


class Evidence(BaseModel):
    evidence_id: str
    source_tool: str
    source_db: str
    evidence_type: str
    evidence_grade: str
    target_id: str
    mol_id: str | None = None
    value: float | int | str | bool | None = None
    unit: str | None = None
    confidence: float = Field(ge=0.0, le=1.0)
    timestamp: str
    metadata: dict[str, Any] = Field(default_factory=dict)
    conflicts: list[str] = Field(default_factory=list)

    @classmethod
    def now(cls, **kwargs):
        return cls(timestamp=datetime.now(timezone.utc).isoformat(), **kwargs)


class EvidenceCollection(BaseModel):
    items: list[Evidence] = Field(default_factory=list)

    def add(self, ev: Evidence) -> None:
        self.items.append(ev)

    def extend(self, evs: Iterable[Evidence]) -> None:
        self.items.extend(evs)

    def filter_by(
        self,
        *,
        source_tool: str | None = None,
        evidence_type: str | None = None,
        grade_at_most: str | None = None,
        grade_at_least: str | None = None,
    ) -> "EvidenceCollection":
        order = {"T1": 1, "T2": 2, "T3": 3, "T4": 4}
        out = []
        for ev in self.items:
            if source_tool and ev.source_tool != source_tool:
                continue
            if evidence_type and ev.evidence_type != evidence_type:
                continue
            if grade_at_most and order[ev.evidence_grade] > order[grade_at_most]:
                continue
            if grade_at_least and order[ev.evidence_grade] < order[grade_at_least]:
                continue
            out.append(ev)
        return EvidenceCollection(items=out)

    def summarize(self) -> dict[str, Any]:
        by_type = defaultdict(list)
        for ev in self.items:
            by_type[ev.evidence_type].append(ev)
        summary = {}
        for et, xs in by_type.items():
            summary[et] = {
                "count": len(xs),
                "sources": sorted({x.source_tool for x in xs}),
                "molecules": sorted({x.mol_id for x in xs if x.mol_id}),
            }
        return summary

    def detect_conflicts(self) -> list[dict[str, Any]]:
        """Detect simple same-molecule same-type disagreements.

        Rule:
        - same `mol_id` + same `evidence_type`
        - from different `source_tool`
        - contradictory interpretation
        """
        groups = defaultdict(list)
        for ev in self.items:
            if ev.mol_id is None:
                continue
            groups[(ev.mol_id, ev.evidence_type)].append(ev)

        conflicts = []
        for (mol_id, et), xs in groups.items():
            if len(xs) < 2:
                continue
            # Simple domain-specific conflict heuristics.
            if et == "binding":
                strong = [x for x in xs if _binding_positive(x)]
                weak = [x for x in xs if _binding_negative(x)]
                if strong and weak:
                    conflicts.append({
                        "mol_id": mol_id,
                        "evidence_type": et,
                        "evidence_ids": [x.evidence_id for x in xs],
                        "reason": "binding_disagreement",
                    })
            elif et == "safety":
                risky = [x for x in xs if _safety_high(x)]
                clean = [x for x in xs if _safety_low(x)]
                if risky and clean:
                    conflicts.append({
                        "mol_id": mol_id,
                        "evidence_type": et,
                        "evidence_ids": [x.evidence_id for x in xs],
                        "reason": "safety_disagreement",
                    })
        return conflicts

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    @classmethod
    def from_json(cls, text: str) -> "EvidenceCollection":
        return cls.model_validate_json(text)


def _binding_positive(ev: Evidence) -> bool:
    if ev.source_tool == "vina_docking" and isinstance(ev.value, (int, float)):
        return float(ev.value) <= -9.0
    if ev.source_tool == "boltz2_affinity" and isinstance(ev.value, (int, float)):
        # interpret as binder_prob in metadata or direct numeric value if explicitly used
        return float(ev.value) >= 0.5
    if ev.source_tool == "prolif_interactions" and isinstance(ev.value, bool):
        return bool(ev.value)
    return False


def _binding_negative(ev: Evidence) -> bool:
    if ev.source_tool == "vina_docking" and isinstance(ev.value, (int, float)):
        return float(ev.value) > -7.0
    if ev.source_tool == "boltz2_affinity" and isinstance(ev.value, (int, float)):
        return float(ev.value) < 0.2
    if ev.source_tool == "prolif_interactions" and isinstance(ev.value, bool):
        return not bool(ev.value)
    return False


def _safety_high(ev: Evidence) -> bool:
    if isinstance(ev.value, str):
        return ev.value.lower() in {"high", "danger", "fail"}
    if isinstance(ev.value, bool):
        return bool(ev.value)
    return False


def _safety_low(ev: Evidence) -> bool:
    if isinstance(ev.value, str):
        return ev.value.lower() in {"low", "clean", "pass"}
    if isinstance(ev.value, bool):
        return not bool(ev.value)
    return False
```

---

## Mapping guidance from existing skills

### skill 9 — chem-docking

Map Vina score to:

```python
Evidence.now(
    evidence_id="cycle6_mol_0064_vina",
    source_tool="vina_docking",
    source_db="PDB",
    evidence_type="binding",
    evidence_grade="T3",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value=-10.04,
    unit="kcal/mol",
    confidence=0.45,
    metadata={"pdb_id": "1VJY"},
)
```

### skill 10 — chem-mlff / MACE prescreen

```python
Evidence.now(
    evidence_id="cycle6_mol_0064_mace_strain",
    source_tool="mace_prescreen",
    source_db="local_model",
    evidence_type="structural",
    evidence_grade="T3",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value=18.2,
    unit="kcal/mol",
    confidence=0.5,
)
```

### skill 14 — chem-reactivity-safety

- SMARTS alert → `type=safety`, `grade=T2`
- ADMET model output → `type=safety`, `grade=T3`

```python
Evidence.now(
    evidence_id="mol_0064_smarts_safety",
    source_tool="smarts_safety",
    source_db="SMARTS_rules",
    evidence_type="safety",
    evidence_grade="T2",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value="pass",
    unit=None,
    confidence=0.8,
)
```

### skill 16 — chem-docking-interactions

```python
Evidence.now(
    evidence_id="mol_0064_hinge_hbond",
    source_tool="prolif_interactions",
    source_db="PDB",
    evidence_type="binding",
    evidence_grade="T2",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value=True,
    unit=None,
    confidence=0.7,
    metadata={"residues": [280, 283]},
)
```

### skill 19 — chem-affinity-prediction (Boltz-2)

For ALK5 specifically, confidence should be **down-weighted**.

```python
Evidence.now(
    evidence_id="mol_0064_boltz_binder_prob",
    source_tool="boltz2_affinity",
    source_db="Boltz-2",
    evidence_type="binding",
    evidence_grade="T3",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value=0.12,
    unit="binder_prob",
    confidence=0.3,  # ALK5-calibrated weak confidence
)
```

### skill 20 — chem-panel-selection

Panel selection should **consume EvidenceCollection**, not invent isolated ad hoc columns.

Recommended pattern:
- load collection
- filter binding/safety/structural evidence for shortlist
- detect conflicts
- summarize per molecule before ranking

### skill 21 — chem-structure-qc-lite (PoseBusters)

```python
Evidence.now(
    evidence_id="mol_0064_pb_valid_postdock",
    source_tool="posebusters_qc",
    source_db="PoseBusters",
    evidence_type="structural",
    evidence_grade="T2",
    target_id="ENSG00000106799",
    mol_id="mol_0064",
    value=True,
    unit=None,
    confidence=0.75,
    metadata={"stage": "post_dock"},
)
```

### skill 22 — chem-target-validation

Each phase can emit evidence with different types / grades.

Examples:
- Open Targets disease association → `disease_association`, often `T2`
- curated clinical records → `clinical`, often `T1/T2`
- structural coverage / AlphaFold → `structural`, usually `T2/T3`
- literature conflicts → `literature`, often `T4`

---

## Conflict rules

Default conflict rule:
- same `mol_id`
- same `evidence_type`
- different `source_tool`
- contradictory interpretation

Example:
- Vina says strong binding (`-10.0`)
- Boltz-2 says weak binding (`binder_prob=0.12`)
- → emit conflict

Represent it as:

```json
{
  "mol_id": "mol_0064",
  "evidence_type": "binding",
  "evidence_ids": [
    "cycle6_mol_0064_vina",
    "mol_0064_boltz_binder_prob"
  ],
  "reason": "binding_disagreement"
}
```

Important:
- conflict does **not** mean one source is wrong
- conflict means the molecule needs panel-style adjudication

---

## Practical conventions

### Confidence

Use confidence to encode method reliability **for this target/context**, not generic faith.

Examples:
- experimental IC50: `0.9–1.0`
- curated SMARTS alert: `0.7–0.9`
- ProLIF hinge vs co-crystal: `0.6–0.8`
- Vina score: `0.3–0.5`
- Boltz-2 on ALK5: `0.3`

### Units

Keep units explicit:
- Vina: `kcal/mol`
- IC50: `nM`, `uM`
- binder probability: `binder_prob`
- boolean pass/fail: `None`

### Metadata

Use `metadata` for everything awkward but useful:
- pdb id
- box definition
- seed
- model version
- disease id
- residue list
- exact raw response fragment if needed

---

## Reporting checklist

- [ ] each emitted evidence has a stable `evidence_id`
- [ ] `source_tool` and `source_db` are explicit
- [ ] evidence grade is not inflated
- [ ] confidence reflects real method reliability
- [ ] conflicts are detected, not hidden
- [ ] original raw outputs remain separately stored

---

## One-sentence rule

**If two skills disagree, do not flatten the disagreement away — represent both as Evidence, attach confidence, and let the pipeline reason over the conflict explicitly.**
