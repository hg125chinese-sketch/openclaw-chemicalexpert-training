---
name: chem-llm
description: Use Large Language Models as chemistry tools — reaction prediction, molecule description, experiment planning, literature-grounded hypothesis generation, and SMILES/IUPAC translation. Covers prompt engineering for chemistry tasks, output validation, and when to trust (and not trust) LLM chemistry.
homepage: https://docs.anthropic.com
metadata: { "openclaw": { "emoji": "🤖", "requires": { "bins": ["curl", "python3"], "python": ["rdkit", "numpy"] } } }
---

# Chemical LLM Applications

Use LLMs as chemistry reasoning tools. This skill covers how to query LLMs for chemical tasks, validate their outputs with RDKit, and integrate LLM reasoning into the DMTA cycle (Design-Make-Test-Analyze).

## When to Use

- User asks CE to "describe this molecule", "explain this reaction", or "suggest an experiment"
- User wants natural-language reasoning about chemistry alongside computational results
- User needs SMILES ↔ IUPAC translation or molecular description
- User wants hypothesis generation grounded in literature (connects to skill 1)
- User wants to combine LLM reasoning with quantitative predictions (skills 2-7)

## Core Philosophy

1. **LLMs reason about chemistry, they don't compute it.** An LLM can explain why a reaction works, but can't reliably predict its yield. Use LLMs for reasoning, RDKit/sklearn/PyTorch for computation.
2. **Always validate LLM chemistry output.** If the LLM generates a SMILES, parse it with RDKit. If it names a reaction, verify the template matches. Never trust unvalidated LLM chemistry.
3. **Prompts are protocols.** A well-structured prompt is like a well-designed assay: it controls variables, specifies output format, and includes positive/negative controls.
4. **LLMs excel at integration.** Their unique value is connecting information across domains — linking a molecule's structure to its pharmacology to its synthesis to its safety. No single computational tool does this.
5. **Uncertainty must be explicit.** LLMs sound confident even when wrong. Always ask for confidence level and flag "LLM-generated, pending validation" in reports.

## Phase 1: Chemistry-Specific Prompting

### 1.1 Molecule Description

```python
def prompt_molecule_description(smiles):
    """Generate a prompt for LLM to describe a molecule."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Pre-compute facts the LLM should incorporate
    context = {
        "smiles": smiles,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": f"{Descriptors.MolWt(mol):.1f}",
        "logp": f"{Descriptors.MolLogP(mol):.2f}",
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rings": Descriptors.RingCount(mol),
        "rotbonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": f"{Descriptors.TPSA(mol):.1f}",
    }

    prompt = f"""Describe this molecule for a medicinal chemist. Be concise (3-4 sentences).

SMILES: {context['smiles']}
Molecular formula: {context['formula']}
MW: {context['mw']} | LogP: {context['logp']} | HBD: {context['hbd']} | HBA: {context['hba']}
Rings: {context['rings']} | RotBonds: {context['rotbonds']} | TPSA: {context['tpsa']}

Cover: (1) key functional groups and scaffold, (2) likely drug-likeness assessment,
(3) any notable structural features for SAR. Do NOT speculate about specific targets
unless the scaffold is very well-known."""

    return prompt
```

### 1.2 Reaction Explanation

```python
def prompt_reaction_explanation(rxn_smiles, conditions=None):
    """Generate a prompt for LLM to explain a reaction mechanism."""

    prompt = f"""Explain this reaction at a graduate chemistry level. Be concise.

Reaction: {rxn_smiles}"""

    if conditions:
        prompt += f"""
Conditions: solvent={conditions.get('solvent', 'N/A')}, 
  catalyst={conditions.get('catalyst', 'N/A')}, 
  temp={conditions.get('temp', 'N/A')}°C"""

    prompt += """

Cover: (1) reaction type/name, (2) key mechanistic steps (2-3 sentences),
(3) role of catalyst/reagent, (4) expected regio/stereoselectivity if relevant.
If you're unsure about any aspect, say so explicitly."""

    return prompt
```

### 1.3 Experiment Planning

```python
def prompt_experiment_plan(target_smiles, goal, constraints=None):
    """Generate a prompt for LLM to suggest an experiment plan."""

    prompt = f"""You are helping plan a chemistry experiment.

Target molecule: {target_smiles}
Goal: {goal}"""

    if constraints:
        prompt += f"\nConstraints: {constraints}"

    prompt += """

Provide a structured plan with:
1. **Objective** (1 sentence)
2. **Proposed route** (numbered steps with reaction type)
3. **Key conditions** for each step (solvent, catalyst, temp, time)
4. **Anticipated challenges** (selectivity issues, protecting groups, purification)
5. **Suggested characterization** (TLC, NMR, MS checkpoints)
6. **Estimated timeline** (days)

Be specific and practical. Flag any steps where you're uncertain about feasibility."""

    return prompt
```

### 1.4 Hypothesis Generation

```python
def prompt_hypothesis(observation, context=None):
    """Generate a prompt for LLM to propose hypotheses from experimental data."""

    prompt = f"""Based on this experimental observation, propose 2-3 testable hypotheses.

Observation: {observation}"""

    if context:
        prompt += f"\nContext: {context}"

    prompt += """

For each hypothesis:
1. **Hypothesis** (1 sentence, falsifiable)
2. **Mechanism** (why this could explain the observation, 2-3 sentences)
3. **Test** (specific experiment to confirm/refute, 1-2 sentences)
4. **Confidence** (low/medium/high, with brief justification)

Rank hypotheses by plausibility. Do NOT propose hypotheses that require information you don't have."""

    return prompt
```

## Phase 2: LLM API Integration

### 2.1 Query via curl (simplest)

```bash
# Example: molecule description via Anthropic API
curl -s https://api.anthropic.com/v1/messages \
  -H "Content-Type: application/json" \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -d '{
    "model": "claude-sonnet-4-20250514",
    "max_tokens": 500,
    "messages": [{"role": "user", "content": "PROMPT_HERE"}]
  }' | jq -r '.content[0].text'
```

### 2.2 Query via Python

```python
import json, subprocess

def query_llm(prompt, model="claude-sonnet-4-20250514", max_tokens=500):
    """Query LLM via API. Returns text response."""
    import os
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")

    if not api_key:
        # Fallback: use the agent's own reasoning
        return "[LLM API not available — using built-in reasoning]"

    payload = {
        "model": model,
        "max_tokens": max_tokens,
        "messages": [{"role": "user", "content": prompt}]
    }

    try:
        result = subprocess.run(
            ["curl", "-s", "https://api.anthropic.com/v1/messages",
             "-H", "Content-Type: application/json",
             "-H", f"x-api-key: {api_key}",
             "-H", "anthropic-version: 2023-06-01",
             "-d", json.dumps(payload)],
            capture_output=True, text=True, timeout=30
        )
        data = json.loads(result.stdout)
        return data.get("content", [{}])[0].get("text", "[No response]")
    except Exception as e:
        return f"[LLM query failed: {e}]"
```

### 2.3 Batch Chemistry Queries

```python
def describe_molecules_batch(smiles_list, delay=1.0):
    """Describe multiple molecules via LLM with rate limiting."""
    import time
    results = []
    for smi in smiles_list:
        prompt = prompt_molecule_description(smi)
        if prompt is None:
            results.append({"smiles": smi, "description": "[invalid SMILES]"})
            continue
        response = query_llm(prompt)
        results.append({"smiles": smi, "description": response})
        time.sleep(delay)
    return results
```

## Phase 3: Output Validation

### 3.1 Validate LLM-Generated SMILES

```python
from rdkit import Chem

def validate_llm_smiles(llm_output):
    """Extract and validate SMILES from LLM output.

    LLMs often generate SMILES with errors. Always validate.
    """
    import re

    # Try to find SMILES-like strings in the output
    smiles_pattern = r'[A-Za-z0-9@+\-\[\]\(\)\\/#=.]+'
    candidates = re.findall(smiles_pattern, llm_output)

    valid = []
    invalid = []
    for candidate in candidates:
        if len(candidate) < 3:
            continue
        mol = Chem.MolFromSmiles(candidate)
        if mol is not None:
            canon = Chem.MolToSmiles(mol)
            valid.append({"raw": candidate, "canonical": canon})
        else:
            # Try common LLM mistakes: extra spaces, wrong brackets
            cleaned = candidate.strip().replace(" ", "")
            mol = Chem.MolFromSmiles(cleaned)
            if mol is not None:
                valid.append({"raw": candidate, "canonical": Chem.MolToSmiles(mol), "cleaned": True})
            else:
                invalid.append(candidate)

    return {"valid": valid, "invalid": invalid}
```

### 3.2 Validate LLM Reaction Claims

```python
def validate_reaction_claim(llm_claim, rxn_smiles=None):
    """Cross-check an LLM's claim about a reaction against RDKit.

    Args:
        llm_claim: dict with keys like 'reaction_type', 'product', 'mechanism'
        rxn_smiles: optional reaction SMILES for verification
    """
    checks = []

    # Check if claimed product is valid
    if "product" in llm_claim:
        mol = Chem.MolFromSmiles(llm_claim["product"])
        checks.append({
            "claim": f"Product SMILES: {llm_claim['product']}",
            "valid": mol is not None,
            "note": "RDKit parse check"
        })

    # Check if reaction type matches our classification
    if rxn_smiles and "reaction_type" in llm_claim:
        our_types = classify_reaction(rxn_smiles)  # from skill 7
        match = llm_claim["reaction_type"].lower() in [t.lower() for t in our_types]
        checks.append({
            "claim": f"Reaction type: {llm_claim['reaction_type']}",
            "valid": match,
            "note": f"Our classification: {our_types}"
        })

    return checks
```

### 3.3 Confidence Tagging

```python
def tag_llm_output(text, source="LLM"):
    """Tag LLM-generated content with provenance and confidence."""
    return {
        "text": text,
        "source": source,
        "validated": False,
        "confidence": "pending_review",
        "tag": f"[{source}-generated, pending validation]"
    }
```

## Phase 4: DMTA Cycle Integration

### 4.1 Design Phase (LLM + skills 2-5)

```python
def design_phase_prompt(target_profile, known_actives=None):
    """LLM-assisted molecular design brief."""

    prompt = f"""You are assisting in the Design phase of a DMTA cycle.

Target profile:
{json.dumps(target_profile, indent=2)}"""

    if known_actives:
        prompt += f"\n\nKnown active molecules (SMILES): {known_actives[:5]}"

    prompt += """

Suggest 3 design strategies:
1. **Scaffold hopping** — propose alternative cores that maintain key pharmacophore
2. **Fragment growing** — suggest modifications to improve potency/selectivity
3. **Property optimization** — suggest changes to improve ADMET while maintaining activity

For each strategy: (a) rationale, (b) specific structural change, (c) expected effect on activity,
(d) expected effect on ADMET. Be concrete — give SMILES or structural descriptions, not vague suggestions."""

    return prompt
```

### 4.2 Analyze Phase (LLM interpretation of results)

```python
def analyze_results_prompt(experiment_results):
    """LLM-assisted analysis of experimental results."""

    prompt = f"""Analyze these experimental results from a DMTA cycle.

Results:
{json.dumps(experiment_results, indent=2)}

Provide:
1. **Key findings** (2-3 sentences)
2. **SAR insights** (what structural features correlate with activity?)
3. **Surprise findings** (anything unexpected? Why?)
4. **Recommendations for next cycle** (what to make next, what to drop)
5. **Confidence assessment** (how reliable are these conclusions given the data?)

Ground your analysis in the data. Do NOT speculate beyond what the results support."""

    return prompt
```

## Phase 5: Safety & Limitations

### 5.1 What LLMs Are Good At (use confidently)

- Natural-language molecule/reaction descriptions
- Connecting structural features to known pharmacology
- Generating experiment plans and protocols
- Explaining mechanisms at a conceptual level
- Literature-grounded hypothesis generation
- Summarizing SAR trends from tabular data

### 5.2 What LLMs Are Bad At (always validate)

- **SMILES generation**: ~30-60% error rate for complex molecules. ALWAYS validate with RDKit.
- **Yield prediction**: LLMs cannot reliably predict yields. Use heuristics (skill 7) or ML models.
- **Stereochemistry**: LLMs frequently confuse R/S, cis/trans. Verify with RDKit.
- **Quantitative predictions**: pIC50, LogP, TPSA — always compute with RDKit, never trust LLM numbers.
- **Novel chemistry**: LLMs interpolate training data. They can't reliably predict outcomes for truly novel reaction types.

### 5.3 Validation Protocol

```
For EVERY LLM chemistry output, apply this protocol:

1. If output contains SMILES → validate with RDKit
2. If output names a reaction → verify against REACTION_CONDITIONS_DB (skill 7)
3. If output claims a property value → recompute with RDKit
4. If output proposes a mechanism → check if it's consistent with reaction type
5. If output makes a literature claim → verify with chem-literature (skill 1)
6. Tag all LLM outputs as "[LLM-generated, validated/unvalidated]"
```

## Checklist Before Reporting

- [ ] **All SMILES validated**: RDKit parse check on every LLM-generated SMILES?
- [ ] **Quantitative claims recomputed**: Property values from RDKit, not LLM?
- [ ] **Confidence stated**: Every LLM output tagged with confidence level?
- [ ] **Source marked**: Clear distinction between computational results and LLM reasoning?
- [ ] **Validation protocol followed**: Steps 1-6 from Phase 5.3 applied?
- [ ] **Hallucination check**: Any claims that seem too specific or too confident?

## Integration with Knowledge Base

- **Save LLM-generated descriptions** to `research/ai4chem/descriptions/<molecule>.md`
- **Save experiment plans** to `research/ai4chem/plans/<target>-<date>.md`
- **Tag all files**: Include `[LLM-assisted]` in commit messages for LLM-generated content
- **Git commit**: `cd /home/node/.openclaw/workspace-chemicalexpert && git add -A && git commit -m "llm-assisted: <description>"`
