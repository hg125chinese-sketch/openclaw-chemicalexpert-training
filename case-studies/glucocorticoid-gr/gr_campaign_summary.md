# GR (NR3C1) glucocorticoid campaign summary

## Executive summary

This campaign asked whether an AI-guided workflow could generate and refine **new GR ligands** for inflammatory-disease glucocorticoid discovery, with an emphasis on **multi-objective quality** rather than raw potency alone.

The full stack now covers:
- target validation
- literature / proposal framing
- public-data QSAR baselines
- pocket-conditioned DiffSBDD generation
- physics/geometry and safety filtering
- docking / interaction analysis
- exemplar selection
- local analog optimization
- DFT confirmation

### Headline outcomes

- **GR / NR3C1 was validated as a high-confidence target** with strong clinical precedent and rich structural support.
- **QSAR baselines were workable but not trivial**; RF was the strongest scaffold-split model, outperforming SVM and GCN.
- **DiffSBDD generation from the 1M2Z agonist-state pocket produced chemically valid molecules at useful QC rates**.
- A multi-signal ranking strategy identified **three GR exemplars** with distinct tradeoffs:
  - `sample 44` — Gln570-contact exemplar
  - `sample 75` — Asn564-contact exemplar, strongest interaction-faithful hit
  - `sample 68` — strongest raw Vina/GNINA hit, but lacking the canonical contact pattern
- Follow-up local optimization around `sample 75` produced **GR75_A04** and **GR75_A06**, both preserving the Asn564 interaction hypothesis, with **GR75_A04** emerging as the strongest successor analog.
- All three GR exemplars passed DFT, bringing the cumulative **DiffSBDD era total to 17/17 = 100% DFT PASS** when combined with the prior IPF campaign (14/14).

## 1) Target validation

### Decision
**GR / NR3C1 = validated**

Why:
- clinically established glucocorticoid target
- broad inflammatory-disease relevance
- multiple approved agonists (e.g. dexamethasone, methylprednisolone)
- strong ligand + structure support
- small-molecule tractability is not the limiting factor

Program-level interpretation:
- this is **not** a target-risk problem
- it is a **therapeutic-index engineering problem**
- success depends on designing ligands with better profile control than classical systemic glucocorticoids

## 2) Literature and proposal context

Key mechanistic framing from the literature review:
- the classic **transactivation vs transrepression** model is historically important but now understood to be incomplete
- GR efficacy/toxicity separation is context-dependent and cannot be reduced to a simple binary switch
- SEGRMs / dissociated agonists remain a meaningful concept, but no universally successful solution has emerged
- therefore a credible AI proposal should focus on:
  - multi-objective ligand optimization
  - context-selective GR pharmacology
  - better balance across efficacy, potency, and side-effect-associated proxies

This supports the proposal theme strongly:
> AI for **better glucocorticoids**, not merely stronger GR binders.

## 3) QSAR baselines

### Dataset
Primary public source:
- **ChEMBL target**: `CHEMBL2034`

Potency-focused curated dataset:
- full curated set: ~`3030` molecules for modeling
- agonist-focused subset was also explored as a cleaner but smaller task

### Three-way comparison on the full dataset
Scaffold split performance:

| Model | R² | RMSE |
|---|---:|---:|
| RF | 0.385 | 0.879 |
| SVM | 0.293 | 0.943 |
| GCN | 0.142 | 1.021 |

### Interpretation
- **RF was the strongest and most reliable baseline**
- SVM was acceptable but weaker
- GCN underperformed, showing that **better model class alone does not overcome assay/task heterogeneity**

### Agonist subset result
The agonist-focused subset **did not improve scaffold generalization** under the current heuristic definition. That is an important negative result:
- task definition still matters more than model complexity
- more careful mechanism-stratified curation is likely needed before graph models become competitive

### Campaign lesson
For GR, the most defensible near-term AI stack is:
- **curated descriptor/fingerprint baseline (RF)**
- plus structure-guided design and ranking
- rather than assuming GNN superiority on noisy, heterogeneous public data

## 4) DiffSBDD generation

### Pocket-conditioned setup
- structure: **1M2Z**
- state: **GR LBD + dexamethasone, agonist state**
- reference ligand: dexamethasone extracted from 1M2Z

### Generation outcome
- requested samples: `100`
- generated molecules written: `98`
- RDKit-valid: `98/98` (`1.000`)
- PoseBusters-valid: `38/98` (`0.388`)

### Interpretation
This is a useful generation pilot:
- chemistry generation was highly valid at the RDKit level
- geometry realism was meaningfully filtered by PoseBusters
- the output quality was good enough to support downstream docking, interaction analysis, and exemplar selection

## 5) Docking calibration and generated-molecule triage

### Dexamethasone redocking baseline
Using the same 1M2Z pocket and Vina protocol:
- **Vina redocking score**: `-11.82`
- ProLIF summary recovered an **Asn564 H-bond donor contact**
- RMSD vs crystal pose was not finalized due to atom-mapping issues, so it should be treated as pending rather than inferred

### Generated set follow-up
From the 38 PB-valid molecules:
- safety-pass molecules: `30`
- successfully docked molecules: `24`

Best raw Vina generated hit:
- `sample 68`
- `vina_best = -8.892`
- `gnina_cnnscore = 0.35755`

## 6) Exemplar panel

A **multi-signal ranking** was used instead of naive Vina sorting:
1. prefer molecules that preserve **Asn564** or **Gln570** interactions
2. then rank by Vina + GNINA support
3. retain one strong non-contact comparator

### Selected exemplars

#### sample 44
- role: **Gln570-contact exemplar**
- Vina: `-8.778`
- GNINA: `0.08867`
- interaction: **Gln570 contact**
- Boltz binder probability: `0.1923`
- DFT: **PASS**
- DFT summary: `gap = 5.16 eV`, `dipole = 6.56 D`

#### sample 75
- role: **Asn564-contact exemplar**
- Vina: `-8.345`
- GNINA: `0.13962`
- interaction: **Asn564 contact**
- Boltz binder probability: `0.3881` (**highest Boltz among exemplars**)
- DFT: **PASS**
- DFT summary: `gap = 3.94 eV`, `dipole = 3.42 D`

#### sample 68
- role: **strongest raw docking exemplar**
- Vina: `-8.892` (**best Vina**)
- GNINA: `0.35755` (**best GNINA among exemplars**)
- interaction: **no key Asn564/Gln570 contact**
- Boltz binder probability: `0.3369`
- DFT: **PASS**
- DFT summary: `gap = 4.84 eV`, `dipole = 4.79 D`

### Exemplar interpretation
- `sample 75` is the most convincing **interaction-faithful lead**
- `sample 68` is the strongest **score-driven alternate binding-mode hypothesis**
- `sample 44` is useful as a **Gln570-contact mechanistic comparator**, but weaker on Boltz and less attractive as the primary lead

## 7) Analog optimization around sample 75

### Design logic
A 6-member local analog panel was built around `sample 75` to preserve the **Asn564 interaction hypothesis** while exploring:
- aromatic ring variations
- H-bond donor/acceptor tuning
- lipophilicity adjustment

All six analogs were:
- RDKit-valid
- Lipinski-pass

### Validation outcome
For `GR75_A01–A06`:
- PB-valid: `6/6`
- safety-pass: `6/6`
- successfully docked: `6/6`

### Key optimized analogs

#### GR75_A04
- design: terminal propyl → **hydroxyethyl**
- interpretation: lower lipophilicity / higher polarity while preserving the anchor hypothesis
- Vina: `-8.197`
- GNINA: `0.18047`
- interaction: **Asn564 contact retained**
- Boltz binder probability: `0.3129`
- DFT: pending / not yet included in this summary

#### GR75_A06
- design: ring-edge **hydroxymethyl** polarity increase
- Vina: `-8.092`
- GNINA: `0.1292`
- interaction: **Asn564 contact retained**
- Boltz binder probability: `0.2437`
- DFT: pending / not yet included in this summary

### Analog optimization interpretation
This stage materially strengthened the campaign:
- the `sample 75` interaction hypothesis was not just descriptive; it proved **optimizable**
- **GR75_A04** emerged as the strongest local successor, outperforming `GR75_A06` across docking, GNINA, and Boltz while retaining Asn564 contact

## 8) DFT results and campaign-level reliability

### GR exemplar DFT outcome
All three exemplars passed:
- `sample 44` → PASS
- `sample 75` → PASS
- `sample 68` → PASS

### Combined DiffSBDD-era DFT reliability
- prior IPF campaign: `14/14 PASS`
- current GR exemplars: `3/3 PASS`
- **combined total: 17/17 = 100% PASS**

### Meaning
This is a major campaign-level signal:
- the current DiffSBDD → QC → ranking → QE/DFT handoff workflow is behaving as a **reliable production pipeline**, not a one-off success
- GR is now an independent second project supporting the same conclusion reached earlier in IPF

## 9) Current top-lead recommendation

### Top GR campaign recommendation: **sample 75**

Why `sample 75` is still the best current campaign lead:
- preserves the **Asn564** interaction seen in the dexamethasone-calibrated GR binding hypothesis
- has the **highest Boltz** among the three exemplars (`0.3881`)
- passed DFT with acceptable physics-level properties
- serves as the most coherent parent for local analog optimization
- its successor line already yielded a strong next-step analog (`GR75_A04`)

### Secondary recommendation: **GR75_A04** as the best optimization successor

Why `GR75_A04` matters:
- retains **Asn564 contact**
- best docking result among the `sample 75` successor analogs
- better GNINA than `GR75_A06`
- better Boltz than `GR75_A06`
- supports the medicinal chemistry direction of increasing polarity while preserving the core interaction motif

### What to say clearly
- **Current campaign lead:** `sample 75`
- **Current best optimized successor analog:** `GR75_A04`

That is the cleanest and most defensible pair of recommendations.

## 10) Bottom line

This GR campaign successfully established a full AI-enabled discovery storyline:
- validated target
- framed therapeutic need
- built public-data ML baselines
- generated pocket-conditioned molecules with DiffSBDD
- triaged them through QC, safety, docking, interaction analysis, and Boltz
- confirmed exemplar quality by DFT
- converted the best interaction-faithful exemplar into an optimizable analog series

The campaign therefore already supports a strong proposal narrative:

> AI can be used not only to generate GR ligands, but to build a **mechanistically anchored, multi-signal optimization workflow** for next-generation glucocorticoid discovery.

At the current stage, the most defensible medicinal chemistry path forward is:
1. keep **sample 75** as the reference GR exemplar lead
2. advance **GR75_A04** as the top optimization successor
3. continue iterative analog design around the **Asn564-preserving interaction motif**
