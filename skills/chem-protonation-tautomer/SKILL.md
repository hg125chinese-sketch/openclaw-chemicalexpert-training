---
name: chem-protonation-tautomer
description: Enumerate physiologically relevant protonation states and tautomers for drug candidates before docking. Ensures docking inputs reflect the actual molecular species present at biological pH. Covers pKa estimation, protomer/tautomer generation, and multi-state docking input preparation.
homepage: https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.html
metadata: { "openclaw": { "emoji": "⚗️", "requires": { "bins": ["python3"], "python": ["rdkit", "numpy", "pandas"] } } }
---

# Protonation & Tautomer Standardization

A molecule's docking score can swing by 2-3 kcal/mol based solely on its protonation state. If you dock the wrong protomer, your score is wrong — not "approximately wrong," but "comparison-invalidatingly wrong." This skill ensures that docking inputs reflect the molecular species actually present at biological pH.

## When to Use

- **Before any docking** (skill 9): standardize all ligands
- Before docking against targets with **pH-sensitive binding** (kinase hinge H-bonds, protease catalytic residues)
- When a candidate has **ionizable groups** (amines, carboxylic acids, phenols, heterocycles)
- When comparing docking scores across molecules with **different ionization profiles**
- When a docking result seems suspicious (unexpectedly good/bad score for a reasonable molecule)

## Core Philosophy

1. **Wrong protonation = wrong score.** A protonated amine makes a salt bridge; a neutral one doesn't. A deprotonated acid forms an H-bond acceptor; a neutral one donates. These aren't subtle effects — they dominate binding.
2. **Enumerate, don't guess.** At pH 7.4, many groups are partially ionized. Don't pick one state — generate the 2-3 most relevant states and dock them all. Report the best score with the state that produced it.
3. **Tautomers are as important as protomers.** A pyridone vs hydroxypyridine tautomer has completely different H-bond donor/acceptor patterns. Docking one and not the other misses half the story.
4. **Standardize first, then protonate.** Remove salts, neutralize charges from drawing conventions, canonicalize, THEN enumerate physiological states. Order matters.
5. **Document the state.** Every docking result must record which protomer/tautomer produced that score. "Vina = -8.5" without specifying the ionization state is incomplete data.

## Phase 1: Molecule Standardization

### 1.1 Clean and Normalize Input

Always run this before protonation/tautomer enumeration:

```python
#!/opt/conda/envs/chem/bin/python
"""Molecule standardization pipeline."""
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def standardize_mol(smi):
    """Standardize a SMILES string: remove salts, neutralize, canonicalize.

    Steps:
    1. Parse SMILES
    2. Remove salt fragments (keep largest)
    3. Uncharge (neutralize drawing-convention charges)
    4. Canonicalize tautomer (to a consistent starting point)
    5. Return standardized SMILES

    Returns:
        dict with original, standardized SMILES, and changes made
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"original": smi, "standardized": None,
                "valid": False, "changes": ["invalid SMILES"]}

    changes = []

    # Step 1: Remove fragments (salts, counterions)
    remover = rdMolStandardize.FragmentParentNormalizer()
    largest = rdMolStandardize.LargestFragmentChooser()
    mol_clean = largest.choose(mol)
    smi_clean = Chem.MolToSmiles(mol_clean)
    if smi_clean != Chem.MolToSmiles(mol):
        changes.append("removed_salt_fragment")

    # Step 2: Uncharge
    uncharger = rdMolStandardize.Uncharger()
    mol_neutral = uncharger.uncharge(mol_clean)
    smi_neutral = Chem.MolToSmiles(mol_neutral)
    if smi_neutral != smi_clean:
        changes.append("neutralized")

    # Step 3: Canonical tautomer
    enumerator = rdMolStandardize.TautomerEnumerator()
    mol_canon = enumerator.Canonicalize(mol_neutral)
    smi_canon = Chem.MolToSmiles(mol_canon)
    if smi_canon != smi_neutral:
        changes.append("canonical_tautomer")

    return {
        "original": smi,
        "standardized": smi_canon,
        "valid": True,
        "changes": changes if changes else ["no_changes"],
    }
```

## Phase 2: Protonation State Enumeration

### 2.1 Common Ionizable Groups and Their pKa Ranges

Understanding which groups ionize at physiological pH (7.4 ± 1):

```python
IONIZABLE_GROUPS = {
    # Group: (SMARTS, typical pKa range, behavior at pH 7.4)
    "aliphatic_amine": {
        "smarts": "[NX3;H2,H1;!$(NC=O);!$(NS=O)]",
        "pka_range": (9.0, 11.0),
        "at_ph74": "mostly protonated (cationic)",
        "note": "Primary/secondary aliphatic amines — almost always +1 at pH 7.4",
    },
    "aromatic_amine": {
        "smarts": "[NH2]-c",
        "pka_range": (2.0, 5.0),
        "at_ph74": "mostly neutral",
        "note": "Anilines — usually neutral at pH 7.4 (pKa too low)",
    },
    "pyridine_n": {
        "smarts": "n1ccccc1",
        "pka_range": (4.5, 6.5),
        "at_ph74": "mostly neutral, but check",
        "note": "Substitution can raise pKa to ~7. If pKa > 6, enumerate both states.",
    },
    "imidazole_n": {
        "smarts": "[nH]1ccnc1",
        "pka_range": (6.0, 8.0),
        "at_ph74": "PARTIALLY ionized — must enumerate both states",
        "note": "This is the most important group for docking! ~50% protonated at pH 7.4.",
    },
    "carboxylic_acid": {
        "smarts": "C(=O)[OH]",
        "pka_range": (3.5, 5.0),
        "at_ph74": "mostly deprotonated (anionic)",
        "note": "Almost always COO⁻ at physiological pH",
    },
    "phenol": {
        "smarts": "c-[OH]",
        "pka_range": (8.0, 11.0),
        "at_ph74": "mostly neutral",
        "note": "Usually neutral. Electron-withdrawing groups lower pKa — check.",
    },
    "sulfonamide_nh": {
        "smarts": "[NH]-S(=O)(=O)",
        "pka_range": (8.0, 10.0),
        "at_ph74": "mostly neutral",
        "note": "Usually neutral at pH 7.4, but some are partially ionized.",
    },
    "tetrazole": {
        "smarts": "c1nn[nH]n1",
        "pka_range": (4.5, 5.5),
        "at_ph74": "mostly deprotonated (anionic)",
        "note": "Common carboxylic acid bioisostere — treat like COO⁻.",
    },
    "guanidine": {
        "smarts": "[NX3]-C(=[NX3])-[NX3]",
        "pka_range": (12.0, 14.0),
        "at_ph74": "always protonated (cationic)",
        "note": "Strongest base in drug molecules — always +1.",
    },
}
```

### 2.2 Identify Ionizable Sites

```python
def identify_ionizable_groups(smi):
    """Find all ionizable groups in a molecule and predict their state at pH 7.4.

    Returns:
        dict with:
        - smiles: input
        - valid: bool
        - groups: list of identified ionizable groups
        - needs_enumeration: bool (True if any group is partially ionized at pH 7.4)
        - recommendation: str
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"smiles": smi, "valid": False, "groups": [],
                "needs_enumeration": False}

    groups = []
    needs_enum = False

    for name, info in IONIZABLE_GROUPS.items():
        pat = Chem.MolFromSmarts(info["smarts"])
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        if matches:
            pka_lo, pka_hi = info["pka_range"]
            # If pKa range overlaps pH 6.4-8.4, both states are relevant
            partially_ionized = (pka_lo <= 8.4) and (pka_hi >= 6.4)
            if partially_ionized:
                needs_enum = True

            groups.append({
                "name": name,
                "n_matches": len(matches),
                "pka_range": info["pka_range"],
                "at_ph74": info["at_ph74"],
                "partially_ionized": partially_ionized,
                "note": info["note"],
            })

    if needs_enum:
        recommendation = "ENUMERATE: molecule has groups with pKa near 7.4. Dock multiple protomers."
    elif groups:
        recommendation = "ASSIGN: all ionizable groups have clear state at pH 7.4. Use dominant protomer."
    else:
        recommendation = "NEUTRAL: no ionizable groups detected. Dock as-is."

    return {
        "smiles": smi,
        "valid": True,
        "groups": groups,
        "needs_enumeration": needs_enum,
        "recommendation": recommendation,
    }
```

### 2.3 Generate Protomers

```python
def enumerate_protomers(smi, ph=7.4, ph_tolerance=1.0):
    """Generate physiologically relevant protonation states.

    Strategy:
    1. Start from standardized (neutral) molecule
    2. For groups clearly ionized at pH 7.4 (pKa far from 7.4): assign dominant state
    3. For groups partially ionized (pKa within ±1 of pH): generate both states

    Uses RDKit's built-in protomer enumeration when available,
    falls back to rule-based approach.

    Returns:
        dict with:
        - original: input SMILES
        - protomers: list of (SMILES, description) tuples
        - n_protomers: int
    """
    # Try RDKit's Enumerator first
    try:
        from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return {"original": smi, "protomers": [], "n_protomers": 0}

        # Use MolStandardize's protomer enumerator
        enumerator = rdMolStandardize.TautomerEnumerator()
        # Note: RDKit's TautomerEnumerator handles tautomers, not protomers
        # For true protomer enumeration, we use the rule-based approach below
    except ImportError:
        pass

    # Rule-based protomer generation
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"original": smi, "protomers": [], "n_protomers": 0}

    protomers = [(smi, "neutral/input form")]

    # Apply dominant ionization states
    ionizable = identify_ionizable_groups(smi)

    # For each group clearly dominant at pH 7.4, generate that state
    # For each group partially ionized, generate both states
    rw_mol = Chem.RWMol(mol)

    # Aliphatic amines: protonate (add H+)
    pat_aliph_amine = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O);!$(NS=O)]")
    if pat_aliph_amine and mol.HasSubstructMatch(pat_aliph_amine):
        # Generate protonated form
        try:
            prot_smi = _protonate_amines(smi)
            if prot_smi and prot_smi != smi:
                protomers.append((prot_smi, "aliphatic amine protonated (+1)"))
        except Exception:
            pass

    # Carboxylic acids: deprotonate
    pat_cooh = Chem.MolFromSmarts("C(=O)[OH]")
    if pat_cooh and mol.HasSubstructMatch(pat_cooh):
        try:
            deprot_smi = _deprotonate_acid(smi)
            if deprot_smi and deprot_smi != smi:
                protomers.append((deprot_smi, "carboxylate deprotonated (-1)"))
        except Exception:
            pass

    # Imidazole: enumerate both states (critical for kinases)
    pat_imidazole = Chem.MolFromSmarts("[nH]1ccnc1")
    if pat_imidazole and mol.HasSubstructMatch(pat_imidazole):
        try:
            prot_smi = _protonate_imidazole(smi)
            if prot_smi and prot_smi != smi:
                protomers.append((prot_smi, "imidazole protonated (+1) — critical for hinge binding"))
        except Exception:
            pass

    # Deduplicate
    seen = set()
    unique = []
    for s, desc in protomers:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else s
        if canon not in seen:
            seen.add(canon)
            unique.append((canon, desc))

    return {
        "original": smi,
        "protomers": unique,
        "n_protomers": len(unique),
    }


def _protonate_amines(smi):
    """Protonate primary/secondary aliphatic amines."""
    mol = Chem.MolFromSmiles(smi)
    pat = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O);!$(NS=O)]")
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return smi
    rw = Chem.RWMol(mol)
    for match in matches:
        n_idx = match[0]
        atom = rw.GetAtomWithIdx(n_idx)
        atom.SetFormalCharge(1)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        return None


def _deprotonate_acid(smi):
    """Deprotonate carboxylic acids."""
    mol = Chem.MolFromSmiles(smi)
    pat = Chem.MolFromSmarts("C(=O)[OH]")
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return smi
    rw = Chem.RWMol(mol)
    for match in matches:
        oh_idx = match[2]  # The [OH] atom
        atom = rw.GetAtomWithIdx(oh_idx)
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        return None


def _protonate_imidazole(smi):
    """Protonate the unprotonated N in imidazole ring."""
    mol = Chem.MolFromSmiles(smi)
    pat = Chem.MolFromSmarts("[n;!H]1cc[nH]c1")
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return smi
    rw = Chem.RWMol(mol)
    for match in matches:
        n_idx = match[0]  # The non-H nitrogen
        atom = rw.GetAtomWithIdx(n_idx)
        atom.SetFormalCharge(1)
        atom.SetNumExplicitHs(1)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        return None
```

## Phase 3: Tautomer Enumeration

### 3.1 Generate Relevant Tautomers

```python
def enumerate_tautomers(smi, max_tautomers=10):
    """Enumerate tautomers using RDKit's built-in enumerator.

    Key tautomeric systems for drug molecules:
    - Keto/enol (amide/imidic acid)
    - Pyridone/hydroxypyridine
    - Thione/thiol
    - Nitroso/oxime
    - Imine/enamine

    Returns:
        dict with:
        - original: input SMILES
        - canonical: RDKit's canonical tautomer
        - tautomers: list of SMILES strings
        - n_tautomers: int
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"original": smi, "canonical": None,
                "tautomers": [], "n_tautomers": 0}

    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(mol)
    canon_smi = Chem.MolToSmiles(canon)

    tautomers_mols = enumerator.Enumerate(mol)
    tautomers = list(set(
        Chem.MolToSmiles(t) for t in tautomers_mols
    ))[:max_tautomers]

    return {
        "original": smi,
        "canonical": canon_smi,
        "tautomers": tautomers,
        "n_tautomers": len(tautomers),
    }
```

## Phase 4: Multi-State Docking Input Preparation

### 4.1 Prepare Docking Input Set

This is the main entry point — call this before docking:

```python
def prepare_docking_inputs(smi, ph=7.4, max_tautomers=5,
                           include_neutral=True):
    """Prepare a set of protomers and tautomers for docking.

    Pipeline:
    1. Standardize input
    2. Enumerate protomers at pH 7.4
    3. For each protomer, enumerate tautomers
    4. Deduplicate
    5. Return all unique states for docking

    Args:
        smi: input SMILES
        ph: target pH (default 7.4)
        max_tautomers: max tautomers per protomer
        include_neutral: always include the neutral form

    Returns:
        dict with:
        - original: input SMILES
        - standardized: cleaned SMILES
        - docking_inputs: list of (SMILES, state_description) tuples
        - n_states: total number of states to dock
        - ionizable_groups: identified ionizable groups
    """
    # Step 1: Standardize
    std = standardize_mol(smi)
    if not std["valid"]:
        return {"original": smi, "standardized": None,
                "docking_inputs": [], "n_states": 0}

    std_smi = std["standardized"]

    # Step 2: Identify ionizable groups
    ionizable = identify_ionizable_groups(std_smi)

    # Step 3: Enumerate protomers
    prot = enumerate_protomers(std_smi, ph=ph)

    # Step 4: For each protomer, enumerate tautomers
    all_states = []
    seen = set()

    for p_smi, p_desc in prot["protomers"]:
        tauts = enumerate_tautomers(p_smi, max_tautomers=max_tautomers)
        for t_smi in tauts["tautomers"]:
            if t_smi not in seen:
                seen.add(t_smi)
                desc = f"{p_desc} | tautomer"
                all_states.append((t_smi, desc))

        # Also include the protomer itself (canonical tautomer)
        if p_smi not in seen:
            seen.add(p_smi)
            all_states.append((p_smi, p_desc))

    # Step 5: Always include neutral standardized form
    if include_neutral and std_smi not in seen:
        seen.add(std_smi)
        all_states.append((std_smi, "standardized neutral"))

    return {
        "original": smi,
        "standardized": std_smi,
        "docking_inputs": all_states,
        "n_states": len(all_states),
        "ionizable_groups": ionizable,
    }
```

### 4.2 Batch Preparation

```python
import pandas as pd


def batch_prepare_docking(smiles_list, ph=7.4, max_tautomers=3):
    """Prepare docking inputs for a pool of molecules.

    Returns:
        dict with:
        - summary: pool-level statistics
        - per_molecule: list of prepare_docking_inputs results
        - expanded_list: flat list of all (original_smi, dock_smi, state_desc)
                         ready for docking pipeline
    """
    results = [prepare_docking_inputs(smi, ph=ph, max_tautomers=max_tautomers)
               for smi in smiles_list]

    n_input = len(smiles_list)
    n_valid = sum(1 for r in results if r["n_states"] > 0)
    total_states = sum(r["n_states"] for r in results)
    n_need_enum = sum(1 for r in results
                      if r.get("ionizable_groups", {}).get("needs_enumeration", False))

    # Expansion ratio: how many docking runs per input molecule
    expansion_ratio = total_states / n_valid if n_valid > 0 else 0

    # Flat list for docking pipeline
    expanded = []
    for r in results:
        for dock_smi, desc in r.get("docking_inputs", []):
            expanded.append({
                "original_smi": r["original"],
                "dock_smi": dock_smi,
                "state_description": desc,
            })

    summary = {
        "n_input": n_input,
        "n_valid": n_valid,
        "total_docking_states": total_states,
        "expansion_ratio": round(expansion_ratio, 2),
        "n_need_enumeration": n_need_enum,
        "pct_need_enumeration": round(100.0 * n_need_enum / n_valid, 1) if n_valid > 0 else 0,
    }

    return {
        "summary": summary,
        "per_molecule": results,
        "expanded_list": expanded,
    }
```

## Phase 5: Integration with Docking Pipeline

### 5.1 Decision Tree

```
Before docking (skill 9):
│
├── For each candidate SMILES:
│   ├── standardize_mol() → clean input
│   ├── identify_ionizable_groups() → check if enumeration needed
│   ├── prepare_docking_inputs() → get all states
│   └── Dock ALL states, report BEST score + which state produced it
│
├── In docking results:
│   ├── Record: best_score, best_protomer, best_tautomer
│   ├── If best state ≠ neutral form: FLAG for review
│   │   (score may depend critically on ionization)
│   └── Report expansion ratio in docking summary
│
└── Compute budget check:
    ├── expansion_ratio < 3: normal, proceed
    ├── expansion_ratio 3-5: acceptable, note in report
    └── expansion_ratio > 5: too many states, cap max_tautomers
```

### 5.2 Post-Docking State Attribution

```python
def attribute_best_state(docking_results):
    """Given docking results for multiple states of one molecule,
    identify the best-scoring state and attribute it.

    Args:
        docking_results: list of dicts with keys:
            - dock_smi: SMILES of the state docked
            - state_description: from prepare_docking_inputs
            - vina_score: docking score (more negative = better)

    Returns:
        dict with best state info
    """
    if not docking_results:
        return None

    best = min(docking_results, key=lambda x: x["vina_score"])

    return {
        "best_smi": best["dock_smi"],
        "best_score": best["vina_score"],
        "best_state": best["state_description"],
        "score_range": (min(r["vina_score"] for r in docking_results),
                        max(r["vina_score"] for r in docking_results)),
        "score_spread": (max(r["vina_score"] for r in docking_results) -
                         min(r["vina_score"] for r in docking_results)),
        "n_states_docked": len(docking_results),
        "state_sensitive": abs(max(r["vina_score"] for r in docking_results) -
                               min(r["vina_score"] for r in docking_results)) > 1.0,
    }
```

## Failure Modes

1. **Skipping standardization → inconsistent comparisons.** If molecule A is drawn as a salt and molecule B as free base, their docking scores aren't comparable. Always standardize first.

2. **Docking only the neutral form.** For molecules with aliphatic amines (pKa ~10), the dominant species at pH 7.4 is protonated (+1). Docking only the neutral form misses salt bridges that may contribute 2+ kcal/mol.

3. **Combinatorial explosion.** A molecule with 3 ionizable groups and 4 tautomers each = 64 states. Cap enumeration (max_tautomers=3-5) and focus on groups near pH 7.4.

4. **Trusting pKa predictions too much.** The pKa ranges in this skill are heuristic. Real pKa depends on the full molecular context (inductive effects, solvent, etc.). For critical candidates, use dedicated pKa prediction tools (Epik, MoKa) if available.

5. **Not recording which state was docked.** If you can't trace "score = -8.5 came from the protonated imidazole tautomer," the result is not reproducible. Always log the state.

6. **Over-enumerating for virtual screening.** For a 3000-molecule pool, 5× expansion = 15000 docking runs. For initial screening, enumerate only the top candidates (e.g., Top 50 after first-pass docking with canonical forms).

## Relationship to Other Skills

| Skill | Relationship |
|-------|-------------|
| chem-docking (9) | This skill prepares inputs for docking. Run this BEFORE docking. |
| chem-kinase-sar (13) | Kinase hinge binding is highly sensitive to protonation (imidazole, aminopyrimidine). Always enumerate for kinase targets. |
| chem-admet (6) | LogP and TPSA change with ionization state. Use standardized neutral form for ADMET calculations. |
| chem-experiment (12) | Report expansion ratio and state sensitivity in cycle summaries. |

## One-Sentence Rule

**Never dock a single protonation state — enumerate protomers and tautomers at pH 7.4, dock them all, and report the best score with the state that produced it.**
