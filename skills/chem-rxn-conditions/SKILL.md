---
name: chem-rxn-conditions
description: Predict and recommend reaction conditions (solvent, catalyst, temperature, reagents, yield) for organic reactions. Complements chem-retrosynthesis (skill 4) which answers "what to disconnect" — this skill answers "how to actually run the reaction". Covers condition recommendation, yield prediction, and solvent/catalyst selection.
homepage: https://www.rdkit.org
metadata: { "openclaw": { "emoji": "🌡️", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "numpy", "pandas", "sklearn"] } } }
---

# Reaction Condition Prediction

Given a reaction (reactants → product), predict the optimal conditions: solvent, catalyst, temperature, reagents, and expected yield. This skill closes the gap between retrosynthetic planning (skill 4: "what bonds to break") and practical synthesis ("how to actually run it").

## When to Use

- User has a retrosynthetic route (from skill 4) and needs conditions for each step
- User asks "what solvent/catalyst/temperature should I use for reaction X?"
- User wants to predict yield for a proposed reaction
- User wants to compare conditions across similar reactions
- User asks about reaction feasibility beyond SA Score

## Core Philosophy

1. **Conditions are not an afterthought.** A reaction that "works on paper" (retro template matches) can fail completely with wrong solvent, temperature, or catalyst. Conditions determine whether a reaction gives 90% or 0% yield.
2. **Analogy is the chemist's primary tool.** The best predictor of conditions for a new reaction is conditions that worked for similar reactions. Similarity = same reaction type + similar substrates.
3. **Yield prediction is probabilistic, not deterministic.** Report ranges, not point estimates. A model that says "yield = 73%" is lying; "yield likely 60-85%" is honest.
4. **Solvent is half the battle.** Solvent affects rate, selectivity, solubility, workup, and safety. Getting solvent right often matters more than getting temperature right.
5. **Data quality trumps model complexity.** Reaction databases are noisy (incomplete conditions, unreported failures). Acknowledge limitations; don't pretend precision you don't have.

## Phase 1: Reaction Type Classification

### 1.1 Common Reaction Types & Typical Conditions

```python
#!/opt/conda/envs/chem/bin/python
"""Reaction type → default condition lookup table."""

# Curated from standard organic chemistry practice.
# Each entry: reaction_type → {solvents, catalysts, temp_range, typical_yield_range, key_notes}

REACTION_CONDITIONS_DB = {
    "amide_formation": {
        "subtypes": {
            "acid_chloride_amine": {
                "solvents": ["DCM", "THF", "DMF"],
                "base": ["Et3N", "DIPEA", "pyridine"],
                "catalyst": None,
                "temp_C": (0, 25),
                "yield_range": (70, 95),
                "time_h": (1, 4),
                "notes": "Add acid chloride slowly at 0°C, warm to RT. Base scavenges HCl."
            },
            "coupling_reagent": {
                "solvents": ["DMF", "DCM", "NMP"],
                "base": ["DIPEA", "NMM"],
                "catalyst": ["HATU", "HBTU", "EDC/HOBt", "T3P"],
                "temp_C": (0, 25),
                "yield_range": (60, 95),
                "time_h": (2, 12),
                "notes": "EDC/HOBt cheapest; HATU best for sterically hindered. Pre-activate acid 5-10 min."
            },
        },
    },
    "suzuki_coupling": {
        "solvents": ["DMF", "dioxane", "THF/H2O", "toluene/H2O"],
        "base": ["K2CO3", "Cs2CO3", "K3PO4", "Na2CO3"],
        "catalyst": ["Pd(PPh3)4", "Pd(dppf)Cl2", "Pd(OAc)2/SPhos"],
        "temp_C": (80, 110),
        "yield_range": (50, 95),
        "time_h": (2, 18),
        "atmosphere": "N2 or Ar",
        "notes": "Degas solvents. Electron-poor ArBr couple faster. Pd(dppf)Cl2 is general-purpose."
    },
    "buchwald_hartwig": {
        "solvents": ["toluene", "dioxane", "THF"],
        "base": ["NaOtBu", "Cs2CO3", "K3PO4"],
        "catalyst": ["Pd2(dba)3/BINAP", "Pd(OAc)2/XPhos", "Pd-PEPPSI"],
        "temp_C": (80, 110),
        "yield_range": (40, 90),
        "time_h": (4, 24),
        "atmosphere": "N2 or Ar",
        "notes": "Ligand choice critical for selectivity. XPhos/BrettPhos for difficult substrates."
    },
    "reductive_amination": {
        "solvents": ["DCM", "MeOH", "DCE", "THF"],
        "reducing_agent": ["NaBH(OAc)3", "NaBH3CN", "NaBH4"],
        "catalyst": None,
        "temp_C": (0, 25),
        "yield_range": (50, 90),
        "time_h": (2, 18),
        "pH": "slightly acidic (AcOH catalyst) or neutral",
        "notes": "NaBH(OAc)3 preferred (selective for iminium over aldehyde). Pre-form imine 30 min."
    },
    "fischer_esterification": {
        "solvents": ["toluene", "benzene", "neat"],
        "catalyst": ["H2SO4", "p-TsOH", "Amberlyst-15"],
        "temp_C": (60, 120),
        "yield_range": (50, 90),
        "time_h": (4, 24),
        "notes": "Dean-Stark trap to remove water drives equilibrium. Not for sterically hindered acids."
    },
    "acylation_phenol": {
        "solvents": ["DCM", "pyridine", "neat"],
        "reagent": ["Ac2O", "AcCl"],
        "base": ["Et3N", "pyridine", "DMAP (cat.)"],
        "temp_C": (0, 25),
        "yield_range": (80, 98),
        "time_h": (0.5, 4),
        "notes": "Ac2O + DMAP is mild and high-yielding. AcCl faster but more vigorous."
    },
    "SNAr": {
        "solvents": ["DMF", "DMSO", "NMP"],
        "base": ["K2CO3", "Cs2CO3", "Et3N"],
        "catalyst": None,
        "temp_C": (60, 120),
        "yield_range": (40, 90),
        "time_h": (4, 24),
        "notes": "Requires electron-withdrawing group ortho/para to leaving group. F > Cl as leaving group."
    },
    "wittig_reaction": {
        "solvents": ["THF", "DCM", "toluene"],
        "base": ["n-BuLi", "NaHMDS", "KOtBu"],
        "catalyst": None,
        "temp_C": (-78, 25),
        "yield_range": (40, 85),
        "time_h": (1, 12),
        "notes": "Stabilized ylides → E-alkene; non-stabilized → Z-alkene. Generate ylide at low T."
    },
    "aldol_reaction": {
        "solvents": ["THF", "DCM", "Et2O"],
        "base": ["LDA", "NaHMDS", "Et3N (Mukaiyama)"],
        "catalyst": ["proline (organocatalytic)", "Ti(OiPr)4 (Mukaiyama)"],
        "temp_C": (-78, 25),
        "yield_range": (40, 85),
        "time_h": (1, 12),
        "notes": "Control of syn/anti diastereoselectivity depends on enolization conditions."
    },
    "hydrogenation": {
        "solvents": ["MeOH", "EtOAc", "EtOH", "THF"],
        "catalyst": ["Pd/C", "PtO2", "Rh/Al2O3", "Raney Ni"],
        "temp_C": (20, 50),
        "pressure": "1-4 atm H2 (balloon to Parr shaker)",
        "yield_range": (80, 99),
        "time_h": (1, 12),
        "notes": "Pd/C is general-purpose. Filter through Celite. Watch for over-reduction."
    },
}
```

### 1.2 Classify a Reaction

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def classify_reaction(rxn_smiles):
    """Classify a reaction by type using substructure matching.

    Args:
        rxn_smiles: "reactants>>products" format
    Returns:
        list of matching reaction types
    """
    parts = rxn_smiles.split(">>")
    if len(parts) != 2:
        return ["unknown"]

    reactants_smi, products_smi = parts
    products = [Chem.MolFromSmiles(s) for s in products_smi.split(".") if Chem.MolFromSmiles(s)]
    reactants = [Chem.MolFromSmiles(s) for s in reactants_smi.split(".") if Chem.MolFromSmiles(s)]

    if not products or not reactants:
        return ["unknown"]

    matches = []

    # Check product for characteristic functional groups
    for prod in products:
        # Amide bond in product
        if prod.HasSubstructMatch(Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")):
            matches.append("amide_formation")
        # Biaryl in product (Suzuki indicator)
        if prod.HasSubstructMatch(Chem.MolFromSmarts("[c]-[c]")):
            # Check if reactants had boronic acid
            for r in reactants:
                if r.HasSubstructMatch(Chem.MolFromSmarts("[#6]B(O)O")):
                    matches.append("suzuki_coupling")
                    break
        # Ester in product
        if prod.HasSubstructMatch(Chem.MolFromSmarts("[#6][OX2][CX3](=[OX1])")):
            matches.append("fischer_esterification")
        # Amine (reductive amination indicator)
        if prod.HasSubstructMatch(Chem.MolFromSmarts("[CX4][NX3]")):
            matches.append("reductive_amination")

    return matches if matches else ["unknown"]
```

## Phase 2: Condition Recommendation

### 2.1 Rule-Based Recommendation

```python
def recommend_conditions(reaction_type, substrate_info=None):
    """Recommend conditions for a given reaction type.

    Args:
        reaction_type: key from REACTION_CONDITIONS_DB
        substrate_info: optional dict with substrate properties
    Returns:
        dict with recommended conditions
    """
    if reaction_type not in REACTION_CONDITIONS_DB:
        return {"error": f"Unknown reaction type: {reaction_type}",
                "suggestion": "Check REACTION_CONDITIONS_DB keys or classify reaction first"}

    conditions = REACTION_CONDITIONS_DB[reaction_type]

    # Handle subtypes
    if "subtypes" in conditions:
        # Default to most general subtype
        subtype = list(conditions["subtypes"].keys())[0]
        if substrate_info and "subtype" in substrate_info:
            subtype = substrate_info["subtype"]
        conditions = conditions["subtypes"][subtype]

    recommendation = {
        "reaction_type": reaction_type,
        "recommended_solvent": conditions.get("solvents", ["N/A"])[0],
        "alternative_solvents": conditions.get("solvents", [])[1:],
        "catalyst": conditions.get("catalyst"),
        "base": conditions.get("base"),
        "reagent": conditions.get("reagent"),
        "reducing_agent": conditions.get("reducing_agent"),
        "temp_range_C": conditions.get("temp_C"),
        "yield_range_%": conditions.get("yield_range"),
        "time_range_h": conditions.get("time_h"),
        "atmosphere": conditions.get("atmosphere"),
        "notes": conditions.get("notes"),
    }

    return recommendation

def print_recommendation(rec):
    """Pretty-print a condition recommendation."""
    print(f"\n=== Conditions for: {rec['reaction_type']} ===")
    print(f"  Solvent:    {rec['recommended_solvent']} (alt: {', '.join(rec.get('alternative_solvents', []))})")
    if rec.get("catalyst"):
        print(f"  Catalyst:   {rec['catalyst']}")
    if rec.get("base"):
        print(f"  Base:       {rec['base']}")
    if rec.get("reagent"):
        print(f"  Reagent:    {rec['reagent']}")
    if rec.get("reducing_agent"):
        print(f"  Red. agent: {rec['reducing_agent']}")
    print(f"  Temp:       {rec['temp_range_C']}°C")
    print(f"  Time:       {rec['time_range_h']} h")
    print(f"  Yield:      {rec['yield_range_%']}%")
    if rec.get("atmosphere"):
        print(f"  Atmosphere: {rec['atmosphere']}")
    print(f"  Notes:      {rec['notes']}")
```

### 2.2 Substrate-Aware Adjustments

```python
def adjust_for_substrate(recommendation, substrate_smiles):
    """Adjust conditions based on substrate properties."""
    mol = Chem.MolFromSmiles(substrate_smiles)
    if mol is None:
        return recommendation

    from rdkit.Chem import Descriptors

    adjustments = []

    # Solubility concerns
    logp = Descriptors.MolLogP(mol)
    if logp > 4:
        adjustments.append("High LogP → consider less polar solvent (toluene, DCM)")
    if logp < -1:
        adjustments.append("Low LogP → consider protic solvent (MeOH, H2O co-solvent)")

    # Steric hindrance
    mw = Descriptors.MolWt(mol)
    if mw > 400:
        adjustments.append("Large substrate → may need longer reaction time or higher temp")

    # Sensitive functional groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
        adjustments.append("Free OH present → may need protection or milder conditions")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
        adjustments.append("Free NH2 present → may compete in nucleophilic reactions")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]")):
        adjustments.append("Nitro group → compatible with most conditions, may reduce under hydrogenation")

    recommendation["substrate_adjustments"] = adjustments
    return recommendation
```

## Phase 3: Yield Prediction

### 3.1 Heuristic Yield Estimation

```python
def estimate_yield_heuristic(reaction_type, substrate_smiles=None, conditions=None):
    """Estimate yield range using heuristics.

    Returns (low, high, confidence) where confidence is 'low'/'medium'/'high'.
    """
    base_conditions = REACTION_CONDITIONS_DB.get(reaction_type)
    if not base_conditions:
        return (20, 80, "low")

    # Get base yield range
    if "subtypes" in base_conditions:
        subtype = list(base_conditions["subtypes"].keys())[0]
        base_yield = base_conditions["subtypes"][subtype].get("yield_range", (30, 80))
    else:
        base_yield = base_conditions.get("yield_range", (30, 80))

    low, high = base_yield
    confidence = "medium"

    # Adjust for substrate complexity
    if substrate_smiles:
        mol = Chem.MolFromSmiles(substrate_smiles)
        if mol:
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            rings = Descriptors.RingCount(mol)
            stereo = len(Chem.FindMolChiralCenters(mol))

            # Complex substrates → lower yield
            if mw > 500: low -= 10; high -= 5
            if rings > 4: low -= 5
            if stereo > 2: low -= 10; confidence = "low"

    # Clamp
    low = max(5, low)
    high = min(99, high)

    return (low, high, confidence)
```

### 3.2 ML Yield Prediction (Reaction Fingerprints)

```python
import numpy as np
from rdkit.Chem import AllChem

def reaction_fingerprint(rxn_smiles, nbits=2048):
    """Compute a reaction difference fingerprint.

    product_fp - reactant_fp (Morgan difference fingerprint).
    """
    parts = rxn_smiles.split(">>")
    if len(parts) != 2:
        return None

    reactant_mols = [Chem.MolFromSmiles(s) for s in parts[0].split(".")]
    product_mols = [Chem.MolFromSmiles(s) for s in parts[1].split(".")]

    reactant_mols = [m for m in reactant_mols if m is not None]
    product_mols = [m for m in product_mols if m is not None]

    if not reactant_mols or not product_mols:
        return None

    # Sum fingerprints for multi-component
    def sum_fps(mols):
        fp_sum = np.zeros(nbits)
        for mol in mols:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)
            fp_sum += np.array(fp)
        return fp_sum

    rxn_fp = sum_fps(product_mols) - sum_fps(reactant_mols)
    return rxn_fp

# Train a yield predictor (sketch — needs reaction dataset with yields)
# from sklearn.ensemble import GradientBoostingRegressor
# X = np.array([reaction_fingerprint(rxn) for rxn in rxn_smiles_list])
# y = np.array(yields)
# model = GradientBoostingRegressor(n_estimators=500).fit(X_train, y_train)
```

## Phase 4: Solvent Selection Guide

### 4.1 Solvent Properties Database

```python
SOLVENT_DB = {
    # name: {bp, polarity(ET30), miscible_water, green_score(1-5, 5=greenest), common_use}
    "DCM":      {"bp": 40,  "ET30": 40.7, "water_misc": False, "green": 1, "use": "general organic"},
    "THF":      {"bp": 66,  "ET30": 37.4, "water_misc": True,  "green": 2, "use": "organometallics, Grignard"},
    "DMF":      {"bp": 153, "ET30": 43.2, "water_misc": True,  "green": 2, "use": "polar aprotic, SNAr, coupling"},
    "DMSO":     {"bp": 189, "ET30": 45.1, "water_misc": True,  "green": 3, "use": "polar aprotic, high-T"},
    "MeOH":     {"bp": 65,  "ET30": 55.4, "water_misc": True,  "green": 4, "use": "protic, reductions, crystallization"},
    "EtOH":     {"bp": 78,  "ET30": 51.9, "water_misc": True,  "green": 5, "use": "protic, green solvent"},
    "toluene":  {"bp": 111, "ET30": 33.9, "water_misc": False, "green": 3, "use": "Pd catalysis, high-T"},
    "dioxane":  {"bp": 101, "ET30": 36.0, "water_misc": True,  "green": 2, "use": "Suzuki, Buchwald"},
    "EtOAc":    {"bp": 77,  "ET30": 38.1, "water_misc": False, "green": 4, "use": "extraction, hydrogenation"},
    "ACN":      {"bp": 82,  "ET30": 45.6, "water_misc": True,  "green": 3, "use": "HPLC, polar reactions"},
    "H2O":      {"bp": 100, "ET30": 63.1, "water_misc": True,  "green": 5, "use": "green, enzyme, some Pd"},
    "NMP":      {"bp": 202, "ET30": 42.2, "water_misc": True,  "green": 2, "use": "high-T amidation, SNAr"},
    "DCE":      {"bp": 83,  "ET30": 41.3, "water_misc": False, "green": 1, "use": "Lewis acid catalysis"},
}

def suggest_solvent(reaction_type, prefer_green=False):
    """Suggest solvents for a reaction type, optionally preferring green solvents."""
    conditions = REACTION_CONDITIONS_DB.get(reaction_type, {})
    candidates = conditions.get("solvents", [])

    if not candidates:
        return ["DMF", "THF"]  # safe defaults

    if prefer_green:
        # Sort by green score descending
        scored = [(s, SOLVENT_DB.get(s, {}).get("green", 3)) for s in candidates]
        scored.sort(key=lambda x: -x[1])
        return [s for s, _ in scored]

    return candidates
```

## Phase 5: Route Annotation (connects to skill 4)

### 5.1 Annotate a Retrosynthetic Route with Conditions

```python
def annotate_route(route_steps):
    """Take a retrosynthetic route (from skill 4) and annotate each step with conditions.

    Args:
        route_steps: list of dicts from skill 4's retro_tree_search
            [{"target": smi, "reaction": type, "precursors": [smi, ...]}, ...]
    Returns:
        annotated route with conditions for each step
    """
    annotated = []
    for step in route_steps:
        rxn_type = step["reaction"]
        conditions = recommend_conditions(rxn_type)

        # Adjust for substrate
        if step.get("precursors"):
            conditions = adjust_for_substrate(conditions, step["precursors"][0])

        # Yield estimate
        yield_est = estimate_yield_heuristic(rxn_type, step.get("target"))

        annotated.append({
            **step,
            "conditions": conditions,
            "estimated_yield": yield_est,
        })

    return annotated
```

### 5.2 Route Conditions Table

Always produce this when annotating a route:

```markdown
| Step | Reaction | Solvent | Catalyst/Reagent | Temp (°C) | Time (h) | Est. Yield | Notes |
|------|----------|---------|-----------------|-----------|----------|------------|-------|
| 1 | suzuki_coupling | dioxane/H2O | Pd(dppf)Cl2, K2CO3 | 90 | 12 | 60-90% | Degas solvents |
| 2 | reductive_amination | DCM | NaBH(OAc)3, AcOH | 0→RT | 4 | 55-85% | Pre-form imine |
| 3 | amide_formation | DMF | HATU, DIPEA | 0→RT | 2 | 70-90% | Pre-activate acid |
```

### 5.3 Overall Route Yield Estimate

```python
def estimate_overall_yield(annotated_route):
    """Estimate overall yield from step yields (multiplicative)."""
    low_total = 1.0
    high_total = 1.0
    for step in annotated_route:
        low, high, _ = step["estimated_yield"]
        low_total *= low / 100
        high_total *= high / 100

    n_steps = len(annotated_route)
    low_pct = low_total * 100
    high_pct = high_total * 100

    print(f"\nOverall route yield estimate ({n_steps} steps):")
    print(f"  Optimistic: {high_pct:.1f}%")
    print(f"  Pessimistic: {low_pct:.1f}%")
    print(f"  Geometric mean: {(low_pct * high_pct)**0.5:.1f}%")

    return {"n_steps": n_steps, "yield_low": low_pct, "yield_high": high_pct}
```

## Checklist Before Reporting

- [ ] **Reaction type classified**: Named reaction identified?
- [ ] **Conditions recommended**: Solvent + catalyst/reagent + base + temp + time?
- [ ] **Yield estimated as range**: Not point estimate? Confidence stated?
- [ ] **Substrate adjustments noted**: FG compatibility? Steric issues?
- [ ] **Green alternatives mentioned**: If user cares about sustainability?
- [ ] **Integrated with route**: Each retro step annotated (connects to skill 4)?
- [ ] **Overall yield calculated**: Multiplicative across steps?

## Integration with Knowledge Base

- **Save condition recommendations** to `research/ai4chem/conditions/<reaction-type>.md`
- **Annotate routes** from `research/ai4chem/retrosynthesis/` with conditions
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "conditions: <reaction> annotation"`
