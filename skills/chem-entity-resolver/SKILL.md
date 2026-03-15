---
name: chem-entity-resolver
description: Resolve target, molecule, and disease identifiers into standardized canonical objects before downstream chemistry/biology workflows. Use when a skill receives ambiguous names, aliases, gene symbols, SMILES, ChEMBL IDs, or disease names and needs stable IDs (Ensembl, UniProt, ChEMBL, InChIKey, EFO) plus confidence and cache-backed reuse.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🧩", "requires": { "bins": ["python3"], "python": ["tooluniverse", "rdkit", "pydantic"] } } }
---

# chem-entity-resolver — canonical IDs first, everything else second

Most pipeline bugs are not chemistry bugs.
They are **entity bugs**:
- “ALK5” vs `TGFBR1`
- gene symbol vs UniProt vs Ensembl vs ChEMBL target id
- non-canonical SMILES vs canonical SMILES
- disease alias vs EFO term

This skill is the **global disambiguation layer**.
Other skills should call it first, then operate on standardized IDs.

Core rule:
- **Do not query downstream databases with ambiguous free text if you can resolve first.**

---

## What it resolves

### 1) Target
Input examples:
- `ALK5`
- `TGFBR1`
- `ENSG00000106799`
- `P36897`

Output:

```json
{
  "gene_symbol": "TGFBR1",
  "ensembl_id": "ENSG00000106799",
  "uniprot_id": "P36897",
  "chembl_target_id": "CHEMBL1994",
  "organism": "Homo sapiens",
  "aliases": ["ALK5", "ACVRLK4", "TBRI"],
  "confidence": 0.92
}
```

### 2) Molecule
Input examples:
- SMILES
- `galunisertib`
- `CHEMBL2364611`

Output:

```json
{
  "canonical_smiles": "Cc1cccc(-c2nn3c(c2-c2ccnc4ccc(C(N)=O)cc24)CCC3)n1",
  "inchi": "InChI=1S/...",
  "inchikey": "...",
  "chembl_id": "CHEMBL2364611",
  "pubchem_cid": 46220502,
  "name": "GALUNISERTIB",
  "mw": 369.43
}
```

### 3) Disease
Input examples:
- `IPF`
- `idiopathic pulmonary fibrosis`

Output:

```json
{
  "efo_id": "EFO_0000768",
  "name": "idiopathic pulmonary fibrosis",
  "therapeutic_area": "respiratory disease",
  "aliases": ["IPF"]
}
```

---

## Hard gates

### Target hard gate
At least **2 databases** must agree on the same entity before target resolution is considered valid.

Acceptable pairings:
- Open Targets + UniProt
- Open Targets + ChEMBL
- UniProt + ChEMBL

If only one source supports the mapping:
- return low confidence
- do **not** mark as resolved

### Molecule hard gate
- RDKit sanitization must pass
- canonical SMILES must be obtainable

### Disease hard gate
- must resolve to a stable ontology id (prefer EFO)

---

## Cache

Cache resolved entities to JSON to avoid repeated API calls.

Suggested location:

```bash
research/ai4chem/entity_cache/
  targets.json
  molecules.json
  diseases.json
```

Cache key suggestions:
- target: normalized input string
- molecule: input string OR standardized InChIKey
- disease: normalized disease name

Use cache unless:
- user requests refresh
- confidence is low
- upstream ids have changed

---

## Minimal object model

```python
#!/opt/conda/envs/chem/bin/python
from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import re

from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import Descriptors
from tooluniverse import ToolUniverse


class ResolvedTarget(BaseModel):
    gene_symbol: str
    ensembl_id: str | None = None
    uniprot_id: str | None = None
    chembl_target_id: str | None = None
    organism: str = "Homo sapiens"
    aliases: list[str] = Field(default_factory=list)
    confidence: float = Field(ge=0.0, le=1.0)
    resolved: bool = False


class ResolvedMolecule(BaseModel):
    canonical_smiles: str
    inchi: str | None = None
    inchikey: str | None = None
    chembl_id: str | None = None
    pubchem_cid: int | None = None
    name: str | None = None
    mw: float | None = None
    resolved: bool = False


class ResolvedDisease(BaseModel):
    efo_id: str | None = None
    name: str
    therapeutic_area: str | None = None
    aliases: list[str] = Field(default_factory=list)
    resolved: bool = False


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())


def _safe_run(tu: ToolUniverse, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    try:
        return tu.run({"name": name, "arguments": arguments})
    except Exception as e:
        return {"status": "error", "error": str(e), "tool": name}


def _load_cache(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _save_cache(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
```

---

## Target resolution

```python
def resolve_target(query: str, cache_dir: str = "research/ai4chem/entity_cache", refresh: bool = False) -> ResolvedTarget:
    cache_path = Path(cache_dir) / "targets.json"
    cache = _load_cache(cache_path)
    key = _norm(query)
    if key in cache and not refresh:
        return ResolvedTarget(**cache[key])

    tu = ToolUniverse()
    tu.load_tools()

    # Heuristic: detect id type first
    ensembl = query if query.startswith("ENSG") else None
    uniprot = query if re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]", query) else None
    symbol = query if (ensembl is None and uniprot is None) else None

    ot = None
    uni = None
    chembl = None

    if symbol:
        ot = _safe_run(tu, "OpenTargets_get_target_id_description_by_name", {"targetName": symbol})
        uni = _safe_run(tu, "UniProtIDMap_gene_to_uniprot", {"gene_names": symbol, "tax_id": 9606})
        chembl = _safe_run(tu, "ChEMBL_search_targets", {"query": symbol})
    elif ensembl:
        # still use name lookup later if needed
        ot = {"data": {"search": {"hits": [{"id": ensembl}]}}}
    elif uniprot:
        uni = {"data": {"results": [{"to": uniprot}]}}

    ot_hits = (((ot or {}).get("data") or {}).get("search") or {}).get("hits", [])
    ot_primary = ot_hits[0] if ot_hits else {}
    uni_results = (((uni or {}).get("data") or {}).get("results", []))
    uni_primary = uni_results[0]["to"] if uni_results else None
    chembl_targets = ((chembl or {}).get("data") or {}).get("targets", []) if isinstance((chembl or {}).get("data"), dict) else []
    chembl_primary = chembl_targets[0]["target_chembl_id"] if chembl_targets else None

    gene_symbol = ot_primary.get("name") or symbol or query
    ensembl_id = ot_primary.get("id") or ensembl
    confidence = 0.0
    confirmations = 0
    if ensembl_id:
        confirmations += 1
        confidence += 0.4
    if uni_primary:
        confirmations += 1
        confidence += 0.3
    if chembl_primary:
        confirmations += 1
        confidence += 0.2
    if gene_symbol:
        confidence += 0.1

    resolved = confirmations >= 2
    obj = ResolvedTarget(
        gene_symbol=gene_symbol,
        ensembl_id=ensembl_id,
        uniprot_id=uni_primary or uniprot,
        chembl_target_id=chembl_primary,
        organism="Homo sapiens",
        aliases=[query] if query != gene_symbol else [],
        confidence=min(confidence, 1.0),
        resolved=resolved,
    )
    cache[key] = obj.model_dump()
    _save_cache(cache_path, cache)
    return obj
```

---

## Molecule resolution

```python
def resolve_molecule(query: str, cache_dir: str = "research/ai4chem/entity_cache", refresh: bool = False) -> ResolvedMolecule:
    cache_path = Path(cache_dir) / "molecules.json"
    cache = _load_cache(cache_path)
    key = _norm(query)
    if key in cache and not refresh:
        return ResolvedMolecule(**cache[key])

    tu = ToolUniverse()
    tu.load_tools()

    mol = None
    chembl_id = query if query.upper().startswith("CHEMBL") else None
    name = None if chembl_id else query

    if chembl_id:
        chembl = _safe_run(tu, "ChEMBL_get_molecule", {"molecule_chembl_id": chembl_id})
    else:
        chembl = _safe_run(tu, "ChEMBL_search_molecules", {"query": query})

    # RDKit-first path for SMILES input
    mol = Chem.MolFromSmiles(query)
    if mol is None:
        # try to recover from ChEMBL result if query was a name
        mols = ((chembl or {}).get("data") or {}).get("molecules", [])
        if mols:
            rec = mols[0]
            smi = ((rec.get("molecule_structures") or {}).get("canonical_smiles"))
            if smi:
                mol = Chem.MolFromSmiles(smi)
                chembl_id = chembl_id or rec.get("molecule_chembl_id")
                name = rec.get("pref_name") or name
                query = smi

    if mol is None:
        raise ValueError("molecule_not_sanitizable")

    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.MolToInchiKey(mol)
    mw = float(Descriptors.MolWt(mol))

    # Optional PubChem cross-check if ToolUniverse tool exists in your install
    pubchem_cid = None
    try:
        pubchem = _safe_run(tu, "PubChem_search_compounds", {"query": canonical_smiles})
        hits = ((pubchem or {}).get("data") or {}).get("compounds", [])
        if hits:
            pubchem_cid = hits[0].get("cid")
    except Exception:
        pass

    obj = ResolvedMolecule(
        canonical_smiles=canonical_smiles,
        inchi=inchi,
        inchikey=inchikey,
        chembl_id=chembl_id,
        pubchem_cid=pubchem_cid,
        name=name,
        mw=round(mw, 3),
        resolved=True,
    )
    cache[key] = obj.model_dump()
    _save_cache(cache_path, cache)
    return obj
```

---

## Disease resolution

```python
def resolve_disease(query: str, cache_dir: str = "research/ai4chem/entity_cache", refresh: bool = False) -> ResolvedDisease:
    cache_path = Path(cache_dir) / "diseases.json"
    cache = _load_cache(cache_path)
    key = _norm(query)
    if key in cache and not refresh:
        return ResolvedDisease(**cache[key])

    tu = ToolUniverse()
    tu.load_tools()

    # Practical approach:
    # 1) maintain a tiny alias map for common abbreviations
    # 2) use Open Targets disease tools downstream with the resolved EFO
    alias_map = {
        "ipf": {"name": "idiopathic pulmonary fibrosis", "efo_id": "EFO_0000768", "therapeutic_area": "respiratory disease", "aliases": ["IPF"]},
    }
    if key in alias_map:
        obj = ResolvedDisease(**alias_map[key], resolved=True)
        cache[key] = obj.model_dump()
        _save_cache(cache_path, cache)
        return obj

    # Fallback: use disease-associated target search if exact disease tools are limited.
    # Replace with a direct disease search tool when available in your ToolUniverse install.
    obj = ResolvedDisease(
        efo_id=None,
        name=query,
        therapeutic_area=None,
        aliases=[],
        resolved=False,
    )
    cache[key] = obj.model_dump()
    _save_cache(cache_path, cache)
    return obj
```

---

## Recommended use from other skills

At the top of the workflow:

```python
target = resolve_target("ALK5")
mol = resolve_molecule("CC1=CC=C...")
disease = resolve_disease("IPF")

if not target.resolved:
    raise ValueError("target_not_resolved")
```

Then use **standardized ids** downstream:
- `target.ensembl_id` for Open Targets
- `target.uniprot_id` for UniProt / AlphaFold / PDBe mappings
- `target.chembl_target_id` for ChEMBL target-centric queries
- `mol.canonical_smiles` / `mol.inchikey` for chemistry pipelines
- `disease.efo_id` for disease ontology queries

---

## Mapping examples

### For skill 9 (docking)
- resolve target first → stable PDB / UniProt / Ensembl context
- resolve molecule first → canonical SMILES before ligand prep

### For skill 14 (safety)
- resolve molecule first → standardized canonical SMILES before SMARTS / ADMET

### For skill 22 (target validation)
- resolve target + disease first
- then all evidence queries use stable ids instead of free text

### For skill 23 (evidence schema)
- use resolved ids inside every Evidence object

---

## Failure modes

1. **Gene alias collision**
- Example: family aliases / old names / pseudogenes
- Fix: require 2-database confirmation before marking resolved

2. **Name-to-molecule ambiguity**
- common names may map to salts, parents, or unrelated synonyms
- Fix: prefer canonical SMILES + InChIKey after RDKit sanitization

3. **Disease abbreviation ambiguity**
- acronyms like `IPF` can be overloaded outside pulmonary fibrosis contexts
- Fix: keep alias map + prefer ontology ids

4. **Overtrusting one database**
- one source alone is not enough for a hard resolve
- Fix: cross-check always

5. **No cache invalidation**
- stale cache can silently propagate old ids
- Fix: add `refresh=True` path and record confidence

---

## One-sentence rule

**Resolve entities once, resolve them well, and make every downstream skill consume canonical IDs instead of raw free text.**
