---
name: chem-docking-interactions
description: Analyze and validate docking poses using protein-ligand interaction fingerprints. Covers PLIP/ProLIF-based interaction profiling, kinase-specific hinge H-bond validation, pose quality control, and interaction-based candidate reranking. Use after docking to convert raw Vina scores into mechanistic evidence.
homepage: https://github.com/pharmai/plip
metadata: { "openclaw": { "emoji": "🔍", "requires": { "bins": ["python3", "plip"], "python": ["rdkit", "prolif", "mdanalysis", "numpy", "pandas", "openbabel"] } } }
---

# Docking Interaction Analysis: From Scores to Mechanistic Evidence

A Vina score of -9.5 means nothing if the ligand isn't making the right contacts. Two molecules can have identical docking scores but completely different interaction profiles — one forming the critical hinge H-bond, the other wedged into the wrong subpocket. This skill turns docking outputs into auditable interaction evidence that a medicinal chemist can evaluate.

## When to Use

- **After every docking run** (skill 9): validate that top-ranked poses make chemically sensible interactions
- When comparing candidates with **similar docking scores** — interaction profiles break ties meaningfully
- When working on **kinase targets**: hinge H-bond is the minimum requirement for Type I inhibitors
- When a docking result seems **suspiciously good or bad** — check if the pose is physically reasonable
- As the final quality gate before presenting candidates to a chemist or writing a cycle report

## Core Philosophy

1. **Scores rank, interactions explain.** A docking score is a noisy single number. The interaction profile tells you WHY the molecule binds (or doesn't). Always report both.
2. **The hinge H-bond is non-negotiable for Type I kinase inhibitors.** If a candidate doesn't form at least one H-bond to the hinge region, it's not a kinase inhibitor — it's an accident. This is the single most important interaction to check.
3. **Compare to the co-crystal, not to perfection.** The reference interaction profile comes from the co-crystal ligand. A good candidate should reproduce most (not necessarily all) of the co-crystal's key contacts.
4. **Interaction fingerprints are more reproducible than scores.** Vina scores fluctuate with minor pose changes. The presence/absence of a specific H-bond or hydrophobic contact is a more stable signal.
5. **Flag, don't just filter.** A missing interaction doesn't always mean a bad molecule. It might mean the wrong protomer was docked (→ skill 15), or the pose needs refinement. Flag it for review.

## Phase 1: Tool Setup

### 1.1 Install PLIP and ProLIF

```bash
# PLIP — rule-based interaction profiler, works on PDB files
pip install plip

# ProLIF — fingerprint-based, works with RDKit molecules and MDAnalysis
pip install prolif

# Both need OpenBabel for some formats
# In the OpenClaw container, openbabel should already be available
# Verify:
python3 -c "from plip.structure.preparation import PDBComplex; print('PLIP OK')"
python3 -c "import prolif; print('ProLIF OK')"
```

### 1.2 When to Use Which Tool

| Tool | Best For | Input | Output |
|------|---------|-------|--------|
| PLIP | Single complex analysis, detailed per-atom report, 8 interaction types | PDB file | XML/TXT report, PyMOL session |
| ProLIF | Batch fingerprinting, docking result comparison, ML features | RDKit Mol + MDAnalysis Universe | Binary fingerprint DataFrame |

**Rule of thumb**: Use PLIP for detailed inspection of top candidates. Use ProLIF for batch comparison across a pool.

## Phase 2: Reference Interaction Profile (Co-Crystal)

### 2.1 Extract Reference from Co-Crystal Structure

Always start by profiling the co-crystal ligand. This is your ground truth.

```python
#!/opt/conda/envs/chem/bin/python
"""Extract reference interaction profile from a co-crystal structure using PLIP."""
from plip.structure.preparation import PDBComplex


def get_reference_interactions(pdb_path, ligand_id=None):
    """Profile all interactions in a co-crystal structure.

    Args:
        pdb_path: path to PDB file (e.g., from RCSB)
        ligand_id: 3-letter ligand code (e.g., "ATP"). If None, profiles all ligands.

    Returns:
        dict with:
        - ligand_id: str
        - interactions: dict keyed by type (hbond, hydrophobic, pi_stack, etc.)
        - residue_contacts: set of residues involved
        - hinge_contacts: list of hinge-region contacts (if kinase)
        - summary: human-readable summary
    """
    mol = PDBComplex()
    mol.load_pdb(pdb_path)
    mol.analyze()

    results = {}
    for bsid, site in sorted(mol.interaction_sets.items()):
        lig_id = bsid.split(":")[0]
        if ligand_id and lig_id != ligand_id:
            continue

        interactions = {
            "hbond_donor": [(h.resnr, h.restype, h.reschain, h.distance_ad,
                             h.type) for h in site.hbonds_pdon],
            "hbond_acceptor": [(h.resnr, h.restype, h.reschain, h.distance_ad,
                                h.type) for h in site.hbonds_ldon],
            "hydrophobic": [(h.resnr, h.restype, h.reschain, h.distance)
                            for h in site.hydrophobic_contacts],
            "pi_stacking": [(s.resnr, s.restype, s.reschain, s.distance,
                             s.type) for s in site.pistacking],
            "pi_cation": [(p.resnr, p.restype, p.reschain, p.distance)
                          for p in site.pication_laro + site.pication_paro],
            "salt_bridge": [(s.resnr, s.restype, s.reschain, s.distance)
                            for s in site.saltbridge_lneg + site.saltbridge_pneg],
            "water_bridge": [(w.resnr, w.restype, w.reschain, w.distance_aw)
                             for w in site.water_bridges],
            "halogen_bond": [(x.resnr, x.restype, x.reschain, x.distance)
                             for x in site.halogen_bonds],
        }

        all_residues = set()
        for itype, ilist in interactions.items():
            for contact in ilist:
                all_residues.add((contact[0], contact[1], contact[2]))

        n_total = sum(len(v) for v in interactions.values())

        results[bsid] = {
            "ligand_id": lig_id,
            "interactions": interactions,
            "residue_contacts": all_residues,
            "n_total_contacts": n_total,
            "n_hbonds": (len(interactions["hbond_donor"]) +
                         len(interactions["hbond_acceptor"])),
            "n_hydrophobic": len(interactions["hydrophobic"]),
        }

    return results
```

### 2.2 Define Kinase-Specific Key Residues

For kinase targets, define the hinge region and other pharmacologically critical residues:

```python
# ALK5/TGFBR1 example (PDB 1VJY)
ALK5_KEY_RESIDUES = {
    "hinge": {
        (280, "HIS"): "hinge backbone NH — primary H-bond donor",
        (281, "ASP"): "hinge+1, H-bond acceptor",
        (283, "ALA"): "hinge backbone CO",
    },
    "gatekeeper": {
        (288, "THR"): "gatekeeper — determines Type I/II selectivity",
    },
    "dfg": {
        (322, "ASP"): "DFG-Asp — catalytic",
        (323, "PHE"): "DFG-Phe — conformational switch",
        (324, "GLY"): "DFG-Gly",
    },
    "catalytic_lysine": {
        (232, "LYS"): "catalytic lysine — salt bridge anchor",
    },
    "p_loop": {
        (210, "GLY"): "P-loop glycine-rich",
    },
}

def make_key_residue_map(key_residues_dict):
    """Flatten the key residue dict into a lookup: (resnr) -> (region, role)."""
    lookup = {}
    for region, residues in key_residues_dict.items():
        for (resnr, restype), role in residues.items():
            lookup[resnr] = {"region": region, "restype": restype, "role": role}
    return lookup
```

## Phase 3: Docking Pose Interaction Analysis

### 3.1 Analyze a Single Docked Pose with PLIP

```python
def analyze_docked_pose(complex_pdb_path, key_residues=None):
    """Analyze a docked ligand-protein complex.

    Args:
        complex_pdb_path: PDB file containing protein + docked ligand
        key_residues: dict from make_key_residue_map() for target-specific checks

    Returns:
        dict with interaction profile + quality assessment
    """
    ref = get_reference_interactions(complex_pdb_path)

    if not ref:
        return {"valid": False, "error": "No ligand interactions detected"}

    # Take first (usually only) ligand
    bsid = list(ref.keys())[0]
    profile = ref[bsid]

    # Key residue coverage
    key_hits = {}
    key_misses = {}
    if key_residues:
        contacted_resnrs = {r[0] for r in profile["residue_contacts"]}
        for resnr, info in key_residues.items():
            if resnr in contacted_resnrs:
                key_hits[resnr] = info
            else:
                key_misses[resnr] = info

    # Hinge H-bond check (critical for kinases)
    has_hinge_hbond = False
    hinge_hbond_details = []
    if key_residues:
        hinge_resnrs = {resnr for resnr, info in key_residues.items()
                        if info["region"] == "hinge"}
        all_hbonds = (profile["interactions"]["hbond_donor"] +
                      profile["interactions"]["hbond_acceptor"])
        for hb in all_hbonds:
            if hb[0] in hinge_resnrs:
                has_hinge_hbond = True
                hinge_hbond_details.append({
                    "resnr": hb[0],
                    "restype": hb[1],
                    "distance": hb[3],
                })

    # Quality assessment
    quality_flags = []
    if key_residues:
        if not has_hinge_hbond:
            quality_flags.append("CRITICAL: No hinge H-bond detected")
        if len(key_hits) < len(key_residues) * 0.3:
            quality_flags.append(
                f"WARNING: Only {len(key_hits)}/{len(key_residues)} "
                "key residues contacted")

    if profile["n_total_contacts"] < 4:
        quality_flags.append("WARNING: Very few total contacts (<4)")
    if profile["n_hbonds"] == 0:
        quality_flags.append("WARNING: No hydrogen bonds detected")

    overall = "good"
    if any("CRITICAL" in f for f in quality_flags):
        overall = "reject"
    elif any("WARNING" in f for f in quality_flags):
        overall = "review"

    return {
        "valid": True,
        "profile": profile,
        "has_hinge_hbond": has_hinge_hbond,
        "hinge_hbond_details": hinge_hbond_details,
        "key_residue_hits": key_hits,
        "key_residue_misses": key_misses,
        "key_residue_coverage": (len(key_hits) / len(key_residues)
                                 if key_residues else None),
        "quality_flags": quality_flags,
        "overall_quality": overall,
    }
```

### 3.2 Batch Analysis with ProLIF Fingerprints

For comparing multiple docked poses across a pool:

```python
import prolif as plf
import MDAnalysis as mda
import pandas as pd


def batch_interaction_fingerprint(protein_pdb, docked_sdf,
                                  ref_ligand_sdf=None):
    """Generate interaction fingerprints for a batch of docked poses.

    Args:
        protein_pdb: path to protein structure (PDB)
        docked_sdf: path to docked ligands (SDF with multiple poses)
        ref_ligand_sdf: optional reference ligand for comparison

    Returns:
        dict with:
        - fingerprint_df: DataFrame of interaction fingerprints
        - similarity_to_ref: Tanimoto similarity to reference (if provided)
        - common_interactions: interactions present in >50% of poses
    """
    # Load protein
    prot = mda.Universe(protein_pdb)
    prot = plf.Molecule.from_mda(prot)

    # Load docked ligands
    lig_suppl = plf.sdf_supplier(docked_sdf)

    # Generate fingerprints
    fp = plf.Fingerprint()
    fp.run_from_iterable(lig_suppl, prot)
    df = fp.to_dataframe()

    result = {
        "fingerprint_df": df,
        "n_poses": len(df),
        "n_interactions": df.shape[1],
    }

    # Reference comparison
    if ref_ligand_sdf:
        ref_suppl = plf.sdf_supplier(ref_ligand_sdf)
        fp_ref = plf.Fingerprint()
        fp_ref.run_from_iterable(ref_suppl, prot)
        ref_bv = fp_ref.to_bitvectors()[0]

        pose_bvs = fp.to_bitvectors()
        similarities = []
        for bv in pose_bvs:
            intersection = (ref_bv & bv).count()
            union = (ref_bv | bv).count()
            sim = intersection / union if union > 0 else 0.0
            similarities.append(sim)

        result["similarity_to_ref"] = similarities
        result["mean_similarity"] = sum(similarities) / len(similarities)

    # Common interactions (present in >50% of poses)
    if len(df) > 1:
        freq = df.mean(axis=0)
        common = freq[freq > 0.5].index.tolist()
        result["common_interactions"] = common
        result["n_common"] = len(common)

    return result
```

## Phase 4: Interaction-Based Reranking

### 4.1 Score Candidates by Interaction Quality

```python
def interaction_rerank(candidates, key_residues, weights=None):
    """Rerank docking candidates by interaction quality.

    Args:
        candidates: list of dicts with keys:
            - smiles: str
            - vina_score: float
            - pose_analysis: output of analyze_docked_pose()
        key_residues: key residue map for the target
        weights: dict of scoring weights (default provided)

    Returns:
        sorted list with interaction_score added
    """
    if weights is None:
        weights = {
            "hinge_hbond": 5.0,       # Critical for kinase
            "key_residue_coverage": 3.0,
            "n_hbonds": 1.0,
            "n_hydrophobic": 0.5,
            "no_critical_flag": 2.0,
        }

    scored = []
    for c in candidates:
        pa = c.get("pose_analysis", {})
        if not pa.get("valid"):
            c["interaction_score"] = -999
            scored.append(c)
            continue

        score = 0.0

        # Hinge H-bond (binary, most important)
        if pa.get("has_hinge_hbond"):
            score += weights["hinge_hbond"]

        # Key residue coverage (0-1 fraction)
        coverage = pa.get("key_residue_coverage", 0) or 0
        score += coverage * weights["key_residue_coverage"]

        # H-bond count (diminishing returns)
        n_hb = pa.get("profile", {}).get("n_hbonds", 0)
        score += min(n_hb, 5) * weights["n_hbonds"]

        # Hydrophobic contacts
        n_hp = pa.get("profile", {}).get("n_hydrophobic", 0)
        score += min(n_hp, 10) * weights["n_hydrophobic"]

        # No critical flags bonus
        if not any("CRITICAL" in f for f in pa.get("quality_flags", [])):
            score += weights["no_critical_flag"]

        c["interaction_score"] = round(score, 2)
        scored.append(c)

    # Sort by interaction score (descending), then by vina score (ascending)
    scored.sort(key=lambda x: (-x["interaction_score"], x.get("vina_score", 0)))
    return scored
```

### 4.2 Composite Score: Vina + Interaction

```python
def composite_score(vina_score, interaction_score,
                    vina_weight=0.5, interaction_weight=0.5):
    """Combine docking score and interaction score into a single ranking metric.

    Both scores are normalized to [0, 1] range before combining.
    More negative vina = better -> normalized to higher = better.

    Args:
        vina_score: Vina docking score (negative, more negative = better)
        interaction_score: from interaction_rerank (positive, higher = better)
        vina_weight: weight for normalized vina score
        interaction_weight: weight for normalized interaction score

    Returns:
        float: composite score (higher = better candidate)
    """
    # Normalize vina: typical range -12 to -4, map to 0-1
    vina_norm = max(0, min(1, (-vina_score - 4.0) / 8.0))

    # Normalize interaction: typical range 0 to 15, map to 0-1
    inter_norm = max(0, min(1, interaction_score / 15.0))

    return round(vina_weight * vina_norm + interaction_weight * inter_norm, 4)
```

## Phase 5: Integration with DMTA Pipeline

### 5.1 Decision Tree

```
After docking (skill 9), before final ranking:
│
├── Step 1: Profile co-crystal -> get_reference_interactions()
│   └── Store as project-level reference (once per target)
│
├── Step 2: For each docked candidate -> analyze_docked_pose()
│   ├── has_hinge_hbond == False?
│   │   ├── For kinase Type I: CRITICAL flag -> reject unless justified
│   │   └── For non-kinase or Type II/III: informational only
│   │
│   ├── key_residue_coverage < 0.3?
│   │   └── WARNING: pose may be in wrong subpocket
│   │
│   └── n_total_contacts < 4?
│       └── WARNING: weak binding, likely bad pose
│
├── Step 3: Rerank using interaction_rerank()
│   └── Composite score = f(vina, interaction_quality)
│
├── Step 4: Top 5 report must include:
│   ├── Vina score
│   ├── Interaction profile (H-bonds, hydrophobic, pi-stack)
│   ├── Hinge H-bond: YES/NO + which residue + distance
│   ├── Key residue coverage fraction
│   └── Quality flags
│
└── In cycle summary (skill 12):
    ├── Report: % of top candidates with hinge H-bond
    ├── Report: mean key residue coverage
    └── Compare interaction profile to co-crystal reference
```

### 5.2 Kinase-Specific Validation Checklist

For any kinase docking project, the Top 5 report MUST answer:

```
[] Does this candidate form a H-bond to the hinge backbone?
   -> If NO: it's not a Type I inhibitor. State this explicitly.
[] Which hinge residue(s) does it contact?
   -> Compare to co-crystal.
[] Does it reach the hydrophobic pocket behind the gatekeeper?
   -> Occupancy of back pocket correlates with potency.
[] Are there any steric clashes (atom overlap < 1.5 A)?
   -> If YES: pose is unreliable. Try redocking or rescoring.
[] Does the interaction profile match the target's known SAR?
   -> Cross-reference with chem-kinase-sar (skill 13).
```

## Failure Modes

1. **Trusting Vina scores without checking interactions.** Vina can give -9.0 to a ligand wedged sideways in the pocket with zero meaningful contacts. Always validate with interaction analysis.

2. **Forgetting to protonate before interaction analysis.** PLIP and ProLIF are sensitive to protonation states. A neutral amine won't show a salt bridge that the protonated form would. Run chem-protonation-tautomer (skill 15) first.

3. **Rigid receptor bias.** Docking into a single crystal conformation may miss interactions that require receptor flexibility (induced fit). If top candidates consistently miss a key contact, consider whether the receptor conformation is appropriate.

4. **Over-weighting interaction count.** More contacts does not equal better binding. A single strong H-bond to the hinge can matter more than 10 weak hydrophobic contacts. Quality over quantity.

5. **Not profiling the co-crystal first.** Without a reference interaction profile, you have no baseline. The co-crystal profile defines "what good looks like" for this target.

6. **PLIP version / hydrogen inconsistency.** PLIP adds hydrogens non-deterministically, which can cause slightly different results between runs. For critical decisions, run analysis multiple times or use ProLIF with explicit hydrogens.

## Relationship to Other Skills

| Skill | Relationship |
|-------|-------------|
| chem-docking (9) | This skill post-processes docking output. Run docking first, then interaction analysis. |
| chem-protonation-tautomer (15) | Interaction detection depends on correct protonation. Always standardize inputs before docking + analysis. |
| chem-kinase-sar (13) | Kinase SAR defines which scaffolds should be present; this skill checks whether those scaffolds make the expected interactions. |
| chem-reactivity-safety (14) | Safety screen filters dangerous molecules; interaction analysis validates that surviving candidates actually bind well. |
| chem-experiment (12) | Interaction profiles are mandatory evidence in cycle reports. "Vina = -9.5" without interaction data is incomplete. |

## One-Sentence Rule

**Never report a docking candidate without its interaction profile — a Vina score without hinge H-bond evidence is not evidence of kinase inhibition.**
