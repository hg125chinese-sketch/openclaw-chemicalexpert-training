---
name: chem-retrosynthesis
description: Perform retrosynthetic analysis on target molecules — disconnect strategies, reaction template extraction, single-step retro prediction, and multi-step route planning. Covers both rule-based (RDKit reaction SMARTS) and ML-based approaches.
homepage: https://www.rdkit.org/docs/Cookbook.html
metadata: { "openclaw": { "emoji": "⚗️", "requires": { "bins": ["python3"], "python": ["rdkit", "torch", "numpy", "pandas"] } } }
---

# Retrosynthetic Analysis

Break down target molecules into purchasable starting materials through systematic disconnection. This skill covers manual analysis, template-based methods, and template-free ML approaches.

## When to Use

- User asks "how would I synthesize molecule X?"
- User provides a target SMILES and wants a synthetic route
- User wants to extract or apply reaction templates
- User asks to evaluate synthesizability of generated molecules
- User wants to build or use a retrosynthesis model

## Core Philosophy

1. **Think backwards.** Retrosynthesis is target → precursors, not forward synthesis. Every step should simplify the molecule or disconnect a strategic bond.
2. **Templates are knowledge, not crutches.** Named reactions encode decades of chemistry. A model that ignores templates is discarding free information; a model that only uses templates can't discover new routes.
3. **Synthesizability is the ultimate filter.** A beautiful molecule from a generative model is worthless if it can't be made. Retrosynthesis closes the loop between design and reality.
4. **Route quality > route existence.** Finding *any* route is easy (enough disconnections reach methane). Finding a route that's short, uses cheap starting materials, and avoids hazardous reagents — that's the real problem.
5. **Validate with forward check.** Every retro prediction should be sanity-checked: do the proposed precursors actually combine to form the target? RDKit reaction SMARTS can verify this.

## Phase 1: Manual Retrosynthetic Analysis

### 1.1 Key Disconnection Strategies

```
Target molecule
├── Identify functional groups (FG analysis)
│   ├── Amide → disconnect to amine + carboxylic acid/acyl chloride
│   ├── Ester → disconnect to alcohol + acid/acyl chloride
│   ├── C-C bond adjacent to carbonyl → aldol / Claisen / Wittig
│   ├── Biaryl → Suzuki / Heck / Buchwald coupling
│   └── C-N bond (aryl) → Buchwald-Hartwig
│
├── Identify strategic bonds
│   ├── Bonds that disconnect rings → simplification
│   ├── Bonds between two complex fragments → convergent
│   └── Bonds adjacent to heteroatoms → polar disconnection
│
└── Apply transforms (retro-reactions)
    └── Target ⟹ Precursor(s) + Reagent/Conditions
```

### 1.2 FG Analysis with RDKit

```python
#!/opt/conda/envs/chem/bin/python
"""Analyze functional groups in a target molecule."""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def analyze_target(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return None

    print(f"Target: {smiles}")
    print(f"  Canonical: {Chem.MolToSmiles(mol)}")
    print(f"  Formula: {rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"  MW: {Descriptors.MolWt(mol):.1f}")
    print(f"  Num rings: {Descriptors.RingCount(mol)}")
    print(f"  Num aromatic rings: {rdMolDescriptors.CalcNumAromaticRings(mol)}")
    print(f"  Num rotatable bonds: {Descriptors.NumRotatableBonds(mol)}")
    print(f"  Num H-bond donors: {Descriptors.NumHDonors(mol)}")
    print(f"  Num H-bond acceptors: {Descriptors.NumHAcceptors(mol)}")
    print(f"  Num chiral centers: {len(Chem.FindMolChiralCenters(mol))}")

    # Functional group SMARTS patterns
    fg_patterns = {
        "amide":          "[NX3][CX3](=[OX1])",
        "ester":          "[OX2][CX3](=[OX1])",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "alcohol":        "[OX2H]",
        "amine_primary":  "[NX3H2]",
        "amine_secondary":"[NX3H1]([#6])[#6]",
        "ketone":         "[#6][CX3](=O)[#6]",
        "aldehyde":       "[CX3H1](=O)[#6]",
        "nitrile":        "[NX1]#[CX2]",
        "nitro":          "[$([NX3](=O)=O)]",
        "sulfonamide":    "[SX4](=[OX1])(=[OX1])([NX3])",
        "biaryl":         "c-c",
        "aryl_halide":    "[cX3][F,Cl,Br,I]",
        "vinyl":          "[CX3]=[CX3]",
    }

    print("\n  Functional groups found:")
    found_any = False
    for name, smarts in fg_patterns.items():
        pat = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pat)
        if matches:
            print(f"    {name}: {len(matches)} occurrence(s)")
            found_any = True
    if not found_any:
        print("    (none from standard set)")

    return mol

# Example
analyze_target("CC(=O)Oc1ccc(CC(=O)O)cc1")  # a simple ester-acid
```

### 1.3 Identify Disconnection Sites

```python
def suggest_disconnections(mol):
    """Suggest strategic bonds to disconnect based on FG analysis."""
    suggestions = []

    # Pattern: bond type → retro-reaction name → what it produces
    disconnections = [
        ("[NX3:1][CX3:2](=[OX1])", "amide_hydrolysis",
         "Disconnect amide → amine + carboxylic acid"),
        ("[OX2:1][CX3:2](=[OX1])", "ester_hydrolysis",
         "Disconnect ester → alcohol + carboxylic acid"),
        ("[cX3:1]-[cX3:2]", "suzuki_coupling",
         "Disconnect biaryl → aryl boronic acid + aryl halide"),
        ("[CX3:1](=O)-[CX4:2]", "aldol_retro",
         "Retro-aldol → two carbonyl fragments"),
        ("[NX3:1]-[cX3:2]", "buchwald_hartwig",
         "Disconnect C(ar)-N → amine + aryl halide"),
        ("[CX3:1]=[CX3:2]", "wittig",
         "Retro-Wittig → aldehyde/ketone + phosphonium ylide"),
    ]

    for smarts, name, description in disconnections:
        pat = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pat)
        if matches:
            for match in matches:
                suggestions.append({
                    "reaction": name,
                    "atoms": match,
                    "description": description
                })

    print(f"\nSuggested disconnections ({len(suggestions)} found):")
    for i, s in enumerate(suggestions):
        print(f"  {i+1}. {s['reaction']} at atoms {s['atoms']}: {s['description']}")

    return suggestions
```

## Phase 2: Reaction Templates (SMARTS)

### 2.1 Defining Reaction Templates

```python
from rdkit.Chem import AllChem

# Retro-reaction templates: product SMARTS >> reactant SMARTS
RETRO_TEMPLATES = {
    "amide_formation": "[C:1](=[O:2])-[N:3]>>[C:1](=[O:2])O.[N:3]",
    "ester_formation": "[C:1](=[O:2])-[O:3]>>[C:1](=[O:2])O.[O:3]",
    "suzuki_coupling": "[c:1]-[c:2]>>[c:1]B(O)O.[c:2]Br",
    "reductive_amination": "[C:1]-[N:2]>>[C:1]=O.[N:2]",
    "aldol": "[C:1]-[C:2](-[OH:3])>>[C:1]=O.[C:2]=O",
    "wittig": "[C:1]=[C:2]>>[C:1]=O.[C:2]=[P]",
    "fischer_esterification": "[C:1](=[O:2])-[O:3]-[C:4]>>[C:1](=[O:2])[OH].[C:4][OH]",
}

def apply_retro_template(product_smiles, template_name):
    """Apply a retro-reaction template to a product molecule."""
    mol = Chem.MolFromSmiles(product_smiles)
    if mol is None:
        return []

    template = RETRO_TEMPLATES.get(template_name)
    if template is None:
        print(f"Unknown template: {template_name}")
        return []

    rxn = AllChem.ReactionFromSmarts(template)
    results = rxn.RunReactants((mol,))

    precursor_sets = []
    for product_set in results:
        precursors = []
        for p in product_set:
            try:
                Chem.SanitizeMol(p)
                precursors.append(Chem.MolToSmiles(p))
            except:
                precursors.append("[sanitize_failed]")
        precursor_sets.append(precursors)

    return precursor_sets
```

### 2.2 Template Extraction from Reaction Data

```python
from rdkit.Chem import rdChemReactions

def extract_reaction_template(rxn_smiles, radius=1):
    """Extract a reaction template (reaction SMARTS) from a reaction SMILES.

    Args:
        rxn_smiles: "reactants>>products" format
        radius: atoms around reaction center to include (0=minimal, 1=default, 2=extended)
    """
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return None
        rdChemReactions.RemoveMappingNumbersFromReactions(rxn)
        return AllChem.ReactionToSmarts(rxn)
    except Exception as e:
        print(f"Template extraction failed: {e}")
        return None

# Example: aspirin synthesis
rxn = "CC(=O)O.Oc1ccccc1C(=O)O>>CC(=O)Oc1ccccc1C(=O)O"
template = extract_reaction_template(rxn)
print(f"Extracted template: {template}")
```

### 2.3 Forward Validation

```python
def validate_retro_step(precursors, target_smiles, forward_template=None):
    """Verify that precursors can produce the target via forward reaction."""
    target_mol = Chem.MolFromSmiles(target_smiles)
    target_canon = Chem.MolToSmiles(target_mol)

    if forward_template:
        rxn = AllChem.ReactionFromSmarts(forward_template)
        reactant_mols = tuple(Chem.MolFromSmiles(s) for s in precursors)
        if any(m is None for m in reactant_mols):
            return False, "Invalid precursor SMILES"

        products = rxn.RunReactants(reactant_mols)
        for product_set in products:
            for prod in product_set:
                try:
                    Chem.SanitizeMol(prod)
                    if Chem.MolToSmiles(prod) == target_canon:
                        return True, "Forward validation passed"
                except:
                    continue

    return False, "Forward validation failed (may need different template)"
```

## Phase 3: Synthesizability Scoring

### 3.1 SA Score (Ertl & Schuffenhauer)

```python
# SA Score implementation (simplified — for full version use sascorer.py from RDKit contrib)
from rdkit.Chem import Descriptors, rdMolDescriptors

def quick_sa_features(smi):
    """Quick synthesizability heuristics (not the full SA score)."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    features = {
        "num_rings": Descriptors.RingCount(mol),
        "num_stereo": len(Chem.FindMolChiralCenters(mol)),
        "num_spiro": rdMolDescriptors.CalcNumSpiroAtoms(mol),
        "num_bridgehead": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        "ring_complexity": rdMolDescriptors.CalcNumRings(mol),
        "num_macrocycles": sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) > 8),
        "fraction_sp3": rdMolDescriptors.CalcFractionCSP3(mol),
    }

    # Simple heuristic score (lower = easier to synthesize)
    score = (
        features["num_stereo"] * 1.5 +
        features["num_spiro"] * 2.0 +
        features["num_bridgehead"] * 2.0 +
        features["num_macrocycles"] * 3.0 +
        max(0, features["num_rings"] - 3) * 0.5
    )
    features["heuristic_difficulty"] = score

    return features
```

### 3.2 Full SA Score (RDKit Contrib)

```python
import os, sys

def get_sa_score(mol):
    """Compute Ertl's SA Score (1-10, lower = easier to synthesize)."""
    # Try importing from RDKit contrib
    try:
        from rdkit.Chem import RDConfig
        sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_path not in sys.path:
            sys.path.insert(0, sa_path)
        import sascorer
        return sascorer.calculateScore(mol)
    except ImportError:
        # Fallback: download sascorer.py if not available
        print("SA Score module not found. Install from RDKit contrib or use quick_sa_features.")
        return None

def batch_sa_scores(smiles_list):
    """Compute SA scores for a list of SMILES."""
    scores = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            sa = get_sa_score(mol)
            scores.append({"smiles": smi, "sa_score": sa})
        else:
            scores.append({"smiles": smi, "sa_score": None})
    return scores
```

### 3.3 Synthesizability Tiers

```
SA Score interpretation:
├── 1-3:  Easy — simple, drug-like, few steps needed
├── 3-5:  Moderate — achievable but may require specialized chemistry
├── 5-7:  Hard — multiple challenging steps, protecting groups likely
└── 7-10: Very hard — probably impractical for most labs
```

## Phase 4: Multi-Step Route Planning

### 4.1 Recursive Retrosynthesis (Tree Search)

```python
from collections import deque

# Set of "purchasable" building blocks (simplified — use real catalog in production)
BUILDING_BLOCKS = {
    "c1ccccc1", "CC(=O)O", "CCO", "CC=O", "c1ccc(Br)cc1",
    "OB(O)c1ccccc1", "CC(N)=O", "NCC", "OC(=O)CC(=O)O",
    "c1ccc(N)cc1", "c1ccc(O)cc1", "CC(=O)Cl",
}

def is_purchasable(smi):
    """Check if a molecule is in the building block catalog."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    canon = Chem.MolToSmiles(mol)
    return canon in BUILDING_BLOCKS

def retro_tree_search(target_smi, max_depth=5, max_branches=3):
    """BFS-based retrosynthetic tree search."""
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))

    if is_purchasable(target_canon):
        return {"target": target_canon, "route": "purchasable", "steps": 0}

    queue = deque()
    queue.append({"smi": target_canon, "depth": 0, "history": []})
    visited = {target_canon}
    routes_found = []

    while queue and len(routes_found) < 3:
        node = queue.popleft()

        if node["depth"] >= max_depth:
            continue

        # Try each retro template
        for template_name, template_smarts in RETRO_TEMPLATES.items():
            precursor_sets = apply_retro_template(node["smi"], template_name)

            for precursors in precursor_sets[:max_branches]:
                step = {
                    "target": node["smi"],
                    "reaction": template_name,
                    "precursors": precursors
                }
                new_history = node["history"] + [step]

                # Check if all precursors are purchasable
                all_purchasable = all(is_purchasable(p) for p in precursors)

                if all_purchasable:
                    routes_found.append({
                        "steps": new_history,
                        "n_steps": len(new_history),
                        "all_purchasable": True
                    })
                else:
                    # Continue searching for non-purchasable precursors
                    for p in precursors:
                        if not is_purchasable(p) and p not in visited:
                            visited.add(p)
                            queue.append({
                                "smi": p,
                                "depth": node["depth"] + 1,
                                "history": new_history
                            })

    return routes_found

def print_route(route):
    """Pretty-print a synthesis route."""
    print(f"\n=== Route ({route['n_steps']} steps) ===")
    for i, step in enumerate(route["steps"]):
        precursor_str = " + ".join(step["precursors"])
        tag = " ✓" if all(is_purchasable(p) for p in step["precursors"]) else ""
        print(f"  Step {i+1}: {step['target']}")
        print(f"    via {step['reaction']}")
        print(f"    from: {precursor_str}{tag}")
```

### 4.2 Route Scoring

```python
def score_route(route):
    """Score a synthetic route on multiple criteria."""
    n_steps = route["n_steps"]
    all_precursors = set()
    for step in route["steps"]:
        all_precursors.update(step["precursors"])

    # Criteria
    n_purchasable = sum(1 for p in all_precursors if is_purchasable(p))
    convergence = n_purchasable / max(n_steps, 1)  # higher = more convergent

    # SA scores of intermediates
    sa_scores = []
    for step in route["steps"]:
        mol = Chem.MolFromSmiles(step["target"])
        if mol:
            sa = get_sa_score(mol)
            if sa is not None:
                sa_scores.append(sa)
    avg_sa = sum(sa_scores) / len(sa_scores) if sa_scores else 0

    score = {
        "n_steps": n_steps,
        "n_unique_precursors": len(all_precursors),
        "n_purchasable": n_purchasable,
        "convergence": convergence,
        "avg_intermediate_sa": avg_sa,
        # Lower is better
        "composite": n_steps * 0.3 + avg_sa * 0.3 + (1 - convergence) * 0.4
    }

    return score
```

## Phase 5: ML-Based Retrosynthesis

### 5.1 Single-Step Retro as Sequence-to-Sequence

```python
"""
Architecture sketch for ML retrosynthesis (template-free).

Input:  product SMILES (tokenized)
Output: reactant SMILES (tokenized, dot-separated for multi-reactant)

Model: Transformer encoder-decoder
  - Encoder: product tokens → hidden representation
  - Decoder: autoregressive generation of reactant tokens

Training data: USPTO reaction dataset (forward reactions reversed)
  product >> reactant1.reactant2
"""

# Example: preparing USPTO data for retro
def prepare_retro_data(rxn_smiles_list):
    """Convert forward reactions to retro format."""
    retro_pairs = []
    for rxn in rxn_smiles_list:
        parts = rxn.split(">>")
        if len(parts) != 2:
            continue
        reactants, products = parts
        # Retro: product → reactants
        retro_pairs.append({
            "input": products.strip(),
            "target": reactants.strip()
        })
    return retro_pairs
```

### 5.2 Evaluating Retro Predictions

```python
def evaluate_retro_predictions(predictions, ground_truth):
    """Evaluate single-step retrosynthesis predictions.

    Args:
        predictions: list of lists (top-k predictions per target)
        ground_truth: list of ground truth reactant SMILES
    """
    def canonicalize_rxn(smi):
        """Canonicalize a multi-reactant SMILES (dot-separated)."""
        parts = smi.split(".")
        canon = []
        for p in parts:
            mol = Chem.MolFromSmiles(p.strip())
            if mol:
                canon.append(Chem.MolToSmiles(mol))
            else:
                return None
        return ".".join(sorted(canon))

    top_k_acc = {1: 0, 3: 0, 5: 0, 10: 0}
    n = len(ground_truth)

    for preds, gt in zip(predictions, ground_truth):
        gt_canon = canonicalize_rxn(gt)
        if gt_canon is None:
            n -= 1
            continue

        for k in top_k_acc:
            for pred in preds[:k]:
                pred_canon = canonicalize_rxn(pred)
                if pred_canon == gt_canon:
                    top_k_acc[k] += 1
                    break

    print("\n=== Retrosynthesis Evaluation ===")
    for k, count in top_k_acc.items():
        acc = count / n if n > 0 else 0
        print(f"  Top-{k} accuracy: {acc:.4f} ({count}/{n})")

    return {f"top_{k}": count / n for k, count in top_k_acc.items()}
```

## Phase 6: Integration with Generative Models

### 6.1 Synthesizability Filter for Generated Molecules

```python
def filter_by_synthesizability(generated_smiles, sa_threshold=5.0):
    """Filter generated molecules by synthesizability."""
    results = {"easy": [], "moderate": [], "hard": [], "invalid": []}

    for smi in generated_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results["invalid"].append(smi)
            continue

        sa = get_sa_score(mol)
        if sa is None:
            results["moderate"].append(smi)
        elif sa <= 3.0:
            results["easy"].append(smi)
        elif sa <= sa_threshold:
            results["moderate"].append(smi)
        else:
            results["hard"].append(smi)

    total_valid = len(results["easy"]) + len(results["moderate"]) + len(results["hard"])
    print(f"\nSynthesizability filter (threshold={sa_threshold}):")
    print(f"  Easy (SA≤3):       {len(results['easy'])}")
    print(f"  Moderate (SA≤{sa_threshold}): {len(results['moderate'])}")
    print(f"  Hard (SA>{sa_threshold}):     {len(results['hard'])}")
    print(f"  Invalid:           {len(results['invalid'])}")
    print(f"  Pass rate:         {(len(results['easy']) + len(results['moderate'])) / max(total_valid, 1):.2%}")

    return results
```

### 6.2 Route-Guided Molecule Scoring

```python
def score_molecule_for_synthesis(smi):
    """Comprehensive synthesizability assessment for a single molecule."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    # 1. SA Score
    sa = get_sa_score(mol) or 5.0

    # 2. Quick retro check
    routes = retro_tree_search(smi, max_depth=3, max_branches=2)
    has_route = len(routes) > 0
    best_steps = min((r["n_steps"] for r in routes), default=99)

    # 3. Structural complexity
    features = quick_sa_features(smi)

    result = {
        "smiles": smi,
        "sa_score": sa,
        "has_retro_route": has_route,
        "min_retro_steps": best_steps if has_route else None,
        "n_stereocenters": features["num_stereo"],
        "n_rings": features["num_rings"],
        "difficulty": features["heuristic_difficulty"],
    }

    return result
```

## Checklist Before Reporting

- [ ] **Target analyzed**: FG identification + disconnection suggestions done?
- [ ] **Templates applied**: At least the obvious disconnections tried?
- [ ] **Forward validated**: Do the precursors actually produce the target?
- [ ] **Purchasability checked**: Are starting materials in a catalog?
- [ ] **SA Score reported**: Quantitative synthesizability assessment?
- [ ] **Route scored**: n_steps, convergence, intermediate complexity?
- [ ] **If from generative model**: Synthesizability filter applied before reporting "novel molecules"?

## Reporting Template

```markdown
# Retrosynthesis: <target name or SMILES>

## Target Analysis
- SMILES: ...
- MW: ... | Rings: ... | Stereocenters: ... | SA Score: ...
- Key FGs: ...

## Proposed Route

| Step | Target | Reaction | Precursors | Purchasable? |
|------|--------|----------|------------|-------------|
| 1 | ... | ... | ... | ✓/✗ |
| 2 | ... | ... | ... | ✓/✗ |

## Route Assessment
- Total steps: ...
- All starting materials purchasable: Yes/No
- Key challenges: ...
- Alternative routes considered: ...
```

## Integration with Knowledge Base

- **Save analyses** to `research/ai4chem/retrosynthesis/<target-slug>.md`
- **Save template libraries** to `research/ai4chem/retrosynthesis/templates/`
- **Cross-reference** with generated molecules from chem-molgen experiments
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "retro: <target> analysis"`
