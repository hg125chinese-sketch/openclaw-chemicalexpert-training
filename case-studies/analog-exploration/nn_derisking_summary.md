# N-N de-risking summary — mol_0021 successor family

Status: **Completed through QE DFT**

This case study archives the N-N de-risking campaign launched after the original `mol_0021` analog family proved that the local chemotype neighborhood was real but remained limited by scaffold-level N-N safety concern.

Goal:
- remove the original N-N safety liability
- preserve hinge-compatible binding mode behavior
- test whether the de-risked family survives the full CE↔QE validation stack

---

## Starting point

Problem identified after analog exploration:
- the `mol_0021` family was structurally strong and QE-valid
- but the original bridge motif triggered scaffold-level `n_n_single_bond` / chronic N-N caution

Main question:
- can the project preserve the validated pose family **without** keeping the original N-N motif?

---

## N-N-free redesign

An 8-member N-N-free analog panel was designed with:
- compact C-N linker replacements
- O-containing linker replacements
- direct C-C linker control
- amide / urea linker replacements

All 8 analogs:
- passed RDKit parsing
- satisfied simple Lipinski-style filters
- cleared the original `n_n_single_bond` SMARTS alert

---

## Validation outcomes

### Structural stack
- all 8 passed PoseBusters and safety_pass in the executed validation stack
- 6/8 retained hinge H-bond after docking

### Multi-seed robustness + GNINA
Three N-N-free analogs became the top validated finalists:
- `NNF_02`
- `NNF_05`
- `NNF_07`

All three were:
- hinge-stable
- strict multi-seed robust
- viable under GNINA rescoring

### Boltz-2
- `NNF_02`: binder_prob **0.154**
- `NNF_05`: binder_prob **0.616**
- `NNF_07`: binder_prob **0.178**

Interpretation:
- `NNF_05` is the strongest parent-like Boltz-supported de-risked continuation
- `NNF_02` is the strongest docking-led but conflict-heavy de-risked winner

### QE DFT
All three N-N-free finalists passed QE DFT:
- `NNF_02`: PASS, gap **4.61 eV**, dipole **3.95 D**, walltime **2.4 h**
- `NNF_05`: PASS, gap **4.85 eV**, dipole **3.65 D**, walltime **1.9 h**
- `NNF_07`: PASS, gap **4.96 eV**, dipole **5.14 D**, walltime **1.4 h**

---

## Main finding

The de-risking campaign succeeded.

Not only because the N-N alert was removed, but because the replacement family preserved the real signals that mattered:
- hinge compatibility
- pose robustness
- QE feasibility

Most important project-level shift:
- parent / parent-like winners were around **~2.9 eV** gap
- N-N-free finalists moved into the **4.6–5.0 eV** range

That makes the de-risked family look not just safer by motif policy, but also stronger by current QE stability proxies.

---

## Final reading

### `NNF_05`
- best parent-like de-risked continuation
- strongest Boltz-supported N-N-free analog

### `NNF_02`
- best docking-led de-risked lead
- valuable because it challenges the assumption that low Boltz must mean weak de-risked chemistry

### `NNF_07`
- best-gap N-N-free analog
- strong balanced fallback / medicinal-chemistry-friendly option

---

## Takeaway

The project now has a credible transition path from:
- a safety-constrained but high-performing parent chemotype

to:
- a **validated N-N-free successor family**.
