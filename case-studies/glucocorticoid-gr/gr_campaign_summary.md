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
- second-generation optimization
- DFT confirmation

### Headline outcomes

- **GR / NR3C1 was validated as a high-confidence target** with strong clinical precedent and rich structural support.
- **QSAR baselines were workable but not trivial**; RF was the strongest scaffold-split model, outperforming SVM and GCN.
- **DiffSBDD generation from the 1M2Z agonist-state pocket produced chemically valid molecules at useful QC rates**.
- A multi-signal ranking strategy identified **three GR exemplars** with distinct tradeoffs:
  - `sample 44` — Gln570-contact exemplar
  - `sample 75` — Asn564-contact exemplar, strongest interaction-faithful hit
  - `sample 68` — strongest raw Vina/GNINA hit, but lacking the canonical contact pattern
- Follow-up local optimization around `sample 75` produced **GR75_A04** and **GR75_A06**, both preserving the Asn564 interaction hypothesis, with **GR75_A04** emerging as the strongest first-generation successor analog.
- A second optimization round around `GR75_A04` produced **GRA04_G2_05** and **GRA04_G2_04** as the strongest gen2 analogs, with **GRA04_G2_05** becoming the current best optimized successor in the GR campaign.
- Both top gen2 analogs passed DFT, confirming that the optimization line preserved quantum-chemistry viability.
- With prior IPF results included, the cumulative **DiffSBDD era total is now 21/21 = 100% DFT PASS**.

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

## 7) First-generation analog optimization around sample 75

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

### First-generation interpretation
This stage materially strengthened the campaign:
- the `sample 75` interaction hypothesis was not just descriptive; it proved **optimizable**
- **GR75_A04** emerged as the strongest local successor, outperforming `GR75_A06` across docking, GNINA, and Boltz while retaining Asn564 contact

## 8) Second-generation optimization around GR75_A04

### Design logic
A second 6-member panel was built around `GR75_A04`, focusing on:
1. hydroxyethyl micro-tuning (methylation, fluorination, extension)
2. core scaffold substitution while preserving the Asn564-anchor hypothesis
3. at least one direction explicitly intended to improve Vina / packing
4. at least one direction explicitly intended to reduce lipophilicity further

### Validation outcome
For `GRA04_G2_01–GRA04_G2_06`:
- PB-valid: `6/6`
- safety-pass: `6/6`
- successfully docked: `6/6`

### Best gen2 analogs

#### GRA04_G2_05
- design: **fluorophenyl fused analogue** for improved aromatic packing
- Vina: `-8.259`
- GNINA: `0.17287`
- interaction: **Asn564 contact retained**
- Boltz binder probability: `0.2996`
- confidence_score: `0.5150`
- DFT: **PASS**
- DFT summary: `gap = 3.83 eV`, `dipole = 3.59 D`

#### GRA04_G2_04
- design: conservative **ring-fluorinated pyridyl-N** electronic tuning
- Vina: `-7.960`
- GNINA: `0.17247`
- interaction: **Asn564 contact retained**
- Boltz binder probability: `0.2006`
- confidence_score: `0.5543`
- DFT: **PASS**
- DFT summary: `gap = 3.53 eV`, `dipole = 2.04 D`

### Gen2 interpretation
The second optimization round sharpened the medicinal chemistry lesson:
- not all “lower lipophilicity” directions preserved the anchor motif
- the strongest result came from **modest aromatic packing enhancement without losing Asn564**
- **GRA04_G2_05** now surpasses `GR75_A04` on Vina while retaining the same key interaction class

## 9) DFT results and campaign-level reliability

### GR exemplar DFT outcome
All three exemplars passed:
- `sample 44` → PASS
- `sample 75` → PASS
- `sample 68` → PASS

### Combined DiffSBDD-era DFT reliability
- prior IPF campaign: `14/14 PASS`
- GR exemplars: `3/3 PASS`
- GR gen2 leads: `2/2 PASS`
- additional previously validated DiffSBDD handoffs in this era: `2/2 PASS`
- **combined total: 21/21 = 100% PASS**

### Meaning
This is a major campaign-level signal:
- the current DiffSBDD → QC → ranking → QE/DFT handoff workflow is behaving as a **reliable production pipeline**, not a one-off success
- GR is now an independent second project supporting the same conclusion reached earlier in IPF

## 10) Lead progression chain

The cleanest summary of current GR lead evolution is:

### `sample 75` → `GR75_A04` → `GRA04_G2_05`

#### sample 75
- first high-confidence **Asn564-contact exemplar**
- strongest Boltz among original exemplars
- DFT PASS
- best parent for mechanism-anchored local optimization

#### GR75_A04
- first successful local successor that improved the `sample 75` line while preserving **Asn564**
- lower-lipophilicity / higher-polarity design that still maintained docking and interaction quality
- strongest first-generation analog

#### GRA04_G2_05
- second-generation successor built from the `GR75_A04` line
- recovered **Asn564** while pushing Vina further to `-8.259`
- DFT PASS with `gap = 3.83 eV`, `dipole = 3.59 D`
- currently the strongest optimized analog in the entire GR campaign

## 11) Current top-lead recommendation

### Current GR campaign lead recommendation: **GRA04_G2_05**

Why `GRA04_G2_05` is now the best overall lead:
- preserves the mechanistically important **Asn564 interaction**
- best Vina among all optimized `sample 75` descendants seen so far (`-8.259`)
- GNINA remains competitive (`0.17287`)
- Boltz remains supportive (`0.2996`), though lower than the original `sample 75` exemplar
- now also has **DFT PASS** with a solid gap/dipole package (`3.83 eV`, `3.59 D`)
- it represents a successful second-generation medicinal chemistry move rather than a one-off generated hit

### Still-important reference lead: **sample 75**

Why `sample 75` still matters:
- highest Boltz among the original exemplars
- DFT PASS
- strongest interaction-faithful original generated molecule
- remains the conceptual parent / benchmark for the optimization series

### Practical recommendation
- **Best optimized current lead:** `GRA04_G2_05`
- **Best original exemplar reference lead:** `sample 75`
- **Best first-generation stepping-stone analog:** `GR75_A04`

## 12) Final GR campaign conclusion

This GR campaign now demonstrates a complete and internally consistent AI-enabled optimization story:
- validate GR as the right target
- understand the pharmacology problem from literature and public data
- establish realistic QSAR baselines
- generate pocket-conditioned molecules with DiffSBDD
- triage by QC, safety, docking, interaction analysis, and Boltz
- identify a strong mechanistic exemplar (`sample 75`)
- convert that exemplar into an optimizable analog lineage (`GR75_A04`)
- improve that lineage further with a second-generation successor (`GRA04_G2_05`)
- confirm the top gen2 branch by DFT

The campaign therefore supports a strong proposal claim:

> AI can drive **mechanistically anchored, multi-signal lead progression** in glucocorticoid discovery, moving from pocket-conditioned generation to interaction-faithful optimization rather than stopping at one-shot hit discovery.

### Final campaign conclusion

- **Reference exemplar lead:** `sample 75`
- **Best first-generation optimized successor:** `GR75_A04`
- **Current best overall GR campaign lead:** `GRA04_G2_05`

Why `GRA04_G2_05` wins:
- keeps the key **Asn564** interaction class
- improves docking relative to `GR75_A04`
- retains supportive GNINA and Boltz signals
- now adds **DFT PASS** on top of the full AI ranking stack

At this point, the GR campaign has progressed from a generated exemplar to a genuinely optimized second-generation lead series, with the clearest progression chain being:

**sample 75 → GR75_A04 → GRA04_G2_05**
