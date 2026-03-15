---
name: chem-self-diagnosis
description: Diagnose and recover from failed chemistry workflow scripts before asking for human help. Use when a script, CLI, notebook, or pipeline step fails with stderr/traceback, missing files, broken paths, tool crashes, API timeouts, unrealistic outputs, or other execution faults in DMTA, docking, QC, safety, Boltz, or data-prep workflows.
homepage: https://docs.openclaw.ai
metadata: { "openclaw": { "emoji": "🩺", "requires": { "bins": ["python3", "bash"], "python": ["pandas"] } } }
---

# chem-self-diagnosis — fail better, recover faster

When a chemistry script crashes, do **not** just paste the traceback and stop.

This skill teaches CE to:
- inspect the real failure surface
- classify the failure
- apply known fixes from prior cycle experience
- retry only when the fix is low-risk and auditable
- escalate only when the problem is scientific or genuinely ambiguous

Core rule:
- **diagnose first, retry second, escalate last**

---

## Output contract

Every diagnosis should produce `diagnosis_report.json` with:

```json
{
  "error_type": "ENV_ERROR",
  "traceback_summary": "FileNotFoundError: /opt/conda/envs/chem/bin/python not found",
  "root_cause": "configured interpreter path is invalid on this host",
  "attempted_fixes": [
    "checked whether /opt/conda/envs/chem/bin/python exists",
    "checked fallback interpreter candidates"
  ],
  "fix_result": "recovered",
  "recommendation": "rerun with absolute valid interpreter path"
}
```

Keep it short, specific, and auditable.

---

## Failure classes

Use exactly one primary class.

### 1) `ENV_ERROR`
Environment / dependency / path / permission problem.

Examples:
- `/opt/conda/envs/chem/bin/python` missing
- binary not installed
- bad PATH assumptions
- permission denied
- wrong conda env

### 2) `DATA_ERROR`
Input file missing, malformed, truncated, or semantically incompatible.

Examples:
- expected CSV/SDF/JSON not found
- required column missing
- invalid JSON
- corrupted molecule records

### 3) `RUNTIME_ERROR`
Script logic failed after starting correctly.

Examples:
- uncaught exception in loop
- traceback from bad indexing / None handling
- OOM
- timeout
- one bad molecule crashes whole batch

### 4) `TOOL_ERROR`
External tool / API failed.

Examples:
- Vina `tree.h`
- RDKit kekulization failure during toolkit call path
- ToolUniverse timeout / parameter validation / tool lookup failure
- gnina direct call breaks while wrapper path is required

### 5) `SCIENCE_ERROR`
Execution finished or nearly finished, but the result is scientifically implausible.

Examples:
- Vina score `> 0`
- MACE strain `> 1000 kcal/mol`
- all molecules fail unexpectedly
- all generated molecules collapse to the same chemistry
- geometry or metric clearly outside project reality

Rule:
- `SCIENCE_ERROR` is a **review flag**, not an auto-fix target.

---

## Standard diagnosis flow

When any script fails:

### A. Capture the failure
1. Read **stderr / traceback last 50 lines**.
2. Record the failing command, script, and input/output paths.
3. Ask: did the script fail before work started, mid-loop, or after producing partial outputs?

### B. Classify
Map the failure to one primary class:
- missing path / dependency / env mismatch → `ENV_ERROR`
- missing input / schema mismatch / parse failure → `DATA_ERROR`
- traceback from our own code / timeout / OOM → `RUNTIME_ERROR`
- external binary / API / toolkit crash → `TOOL_ERROR`
- outputs nonsensical despite technical success → `SCIENCE_ERROR`

### C. Apply class-specific checks
Follow the relevant section below.

### D. Decide recovery action
- safe known fix available → apply + retry
- partial-batch failure → isolate + continue
- scientific implausibility → stop and flag for human review

### E. Emit `diagnosis_report.json`
Always write the diagnosis artifact, even when recovery succeeds.

---

## Class-specific workflow

### `ENV_ERROR`

Check in this order:
1. Does the path/binary actually exist?
2. Is the script using the correct absolute path?
3. Is the correct conda env required?
4. Is permission / executable bit the issue?
5. Is an env var required (`PYTHONNOUSERSITE=1`, etc.)?

Typical recovery actions:
- replace implicit binary with known absolute path
- switch to project-approved interpreter path
- add missing env var in the command wrapper
- prefer wrapper script over raw binary

Auto-retry policy:
- allowed if the fix is deterministic and low-risk
- maximum **2 retries**

---

### `DATA_ERROR`

Check in this order:
1. Do all required input files exist?
2. Are they non-empty?
3. Does the schema match expectations (required columns / fields)?
4. Is the file from the expected previous step?
5. Can the missing input be regenerated from the immediately previous step?

Typical recovery actions:
- fix wrong file path
- select the correct prior-step artifact
- validate schema and fail early with a better message
- if upstream artifact is genuinely absent, recommend rerunning the previous step

Do **not** fabricate missing data.

Auto-retry policy:
- only if the correct input already exists and the script pointed to the wrong one
- otherwise stop and recommend upstream regeneration

---

### `RUNTIME_ERROR`

Check in this order:
1. Read traceback carefully; identify the first user-code frame that matters.
2. Decide whether failure is batch-global or molecule-local.
3. If molecule-local, isolate the bad record and avoid killing the whole batch.
4. If timeout/OOM, determine whether scope reduction or checkpointing is possible.

Typical recovery actions:
- add per-molecule `try/except` around fragile operations
- log failure rows instead of aborting the whole job
- checkpoint partial outputs
- reduce batch size / split long loops
- add explicit timeouts around external calls

Auto-retry policy:
- allowed when the patch preserves logic and only improves robustness
- if root logic is ambiguous, stop and ask

---

### `TOOL_ERROR`

Check in this order:
1. Is this a known tool quirk from the issue library below?
2. Are we calling the tool the approved way (wrapper, absolute path, env var)?
3. Is the failure per-input or global?
4. Is there a known fallback mode?

Typical recovery actions:
- apply wrapper / absolute path
- call `load_tools()` once globally for ToolUniverse
- catch per-molecule Vina failures and continue
- switch to fallback layer if API/tooling is flaky but project policy allows it

Auto-retry policy:
- allowed if a known workaround exists
- maximum **2 retries**

---

### `SCIENCE_ERROR`

Check in this order:
1. Is the number outside obvious physical / project bounds?
2. Is this likely due to bad geometry, bad pose, bad parser, or genuine model weirdness?
3. Does the output contradict multiple orthogonal signals?

Response:
- write diagnosis report
- preserve offending outputs
- flag for human review
- recommend the smallest next validation experiment

Auto-retry policy:
- **do not auto-fix**

---

## Known issue library (from Cycle 1–7)

Use these before inventing new theories.

### gnina
Problem:
- direct invocation is fragile in this project

Rule:
- use **`/home/node/.local/bin/gnina-run`**

Diagnosis tag:
- `TOOL_ERROR`

Fix:
- replace raw `gnina` with wrapper and retry

---

### PoseBusters / bust
Problem:
- PATH assumptions break

Rule:
- use **`/home/node/.local/bin/bust`** absolute path

Diagnosis tag:
- `ENV_ERROR` or `TOOL_ERROR`

Fix:
- rewrite command to absolute path

---

### ToolUniverse startup / API calls
Problems:
- `load_tools()` is slow
- tool lookup names can drift
- parameter names can be exacting
- API/network can timeout

Rules:
- call `load_tools()` **once globally**
- reuse the loaded client
- validate exact parameter names before repeated retries
- if project policy allows, fallback to a lower layer when TU fails repeatedly

Diagnosis tag:
- usually `TOOL_ERROR`

---

### Vina `tree.h`
Problem:
- some molecules cause internal Vina failures

Rule:
- do **not** let one molecule crash the whole batch

Diagnosis tag:
- `TOOL_ERROR`

Fix:
- wrap each molecule in `try/except`
- log failure row
- continue batch

---

### DiffSBDD 3D geometry
Problem:
- generated 3D coordinates can show extreme strain and are unsafe for direct QC handoff

Rule:
- do **not** feed raw DiffSBDD coordinates to QE
- use **RDKit ETKDGv3 + MMFF** re-embed before MACE/QE

Diagnosis tag:
- often `SCIENCE_ERROR` if the geometry is absurd, sometimes `TOOL_ERROR` if downstream tools crash on it

---

### Boltz env
Problem:
- user-site pollution can break the run

Rule:
- set **`PYTHONNOUSERSITE=1`**

Diagnosis tag:
- `ENV_ERROR`

Fix:
- add env var and retry

---

### Long docker tasks
Problem:
- foreground exec / short timeout kills legitimate long runs

Rule:
- use detached/background execution plus log files

Diagnosis tag:
- `RUNTIME_ERROR`

Fix:
- relaunch as long-running background job with explicit logs

---

## Auto-recovery rules

### Recover automatically when:
- `ENV_ERROR` + known deterministic fix exists
- `TOOL_ERROR` + known workaround exists
- `RUNTIME_ERROR` is clearly molecule-local and can be isolated safely

### Do not recover automatically when:
- the proposed change alters scientific logic
- the input artifact is genuinely missing and must be regenerated upstream
- the result is scientifically implausible (`SCIENCE_ERROR`)
- the fix would be speculative or destructive

### Retry limits
- maximum **2 retries** per incident
- after second failure, stop and escalate with diagnosis report

---

## Recommended patch patterns

### Pattern 1 — per-molecule shielding
Use when one bad molecule should not kill a batch.

```python
rows = []
for mol_id, smi in items:
    try:
        result = run_one(mol_id, smi)
        rows.append({"mol_id": mol_id, **result, "error": ""})
    except Exception as e:
        rows.append({"mol_id": mol_id, "error": str(e)[:300]})
        continue
```

### Pattern 2 — early input validation
Use before starting long work.

```python
required = [INPUT_CSV, REC_PDBQT, REPORT_JSON]
for p in required:
    if not p.exists():
        raise FileNotFoundError(f"Missing required file: {p}")
```

### Pattern 3 — tool wrapper enforcement
Use project-approved wrappers explicitly.

```python
GNINA_RUN = "/home/node/.local/bin/gnina-run"
BUST = "/home/node/.local/bin/bust"
```

### Pattern 4 — science sanity gate
Use to stop nonsense early.

```python
if vina_score > 0:
    raise ValueError("SCIENCE_ERROR: positive docking score is not credible here")
if strain_kcal is not None and strain_kcal > 1000:
    raise ValueError("SCIENCE_ERROR: extreme geometry strain")
```

---

## Minimal diagnosis procedure (copyable)

1. Read last 50 stderr/traceback lines.
2. Classify into one of the 5 error types.
3. Check the known issue library.
4. Apply the smallest safe fix.
5. Retry at most twice if policy allows.
6. Write `diagnosis_report.json`.
7. If still failing, escalate with:
   - error class
   - root cause
   - what was tried
   - exact next recommendation

---

## What good behavior looks like

Good:
- “This is a `TOOL_ERROR`: gnina was called directly. Replaced with `/home/node/.local/bin/gnina-run`, retried, recovered.”
- “This is a `DATA_ERROR`: Step 5 CSV is missing. It cannot be reconstructed locally; rerun Step 5.”
- “This is a `SCIENCE_ERROR`: MACE strain is 1450 kcal/mol from raw DiffSBDD geometry. Do not auto-send to QE. Re-embed with RDKit first.”

Bad:
- “The script failed, maybe try again later.”
- “Everything crashed.”
- silently changing scientific logic to force success

---

## Bottom line

The job is not just to notice failure.

The job is to:
- classify it correctly
- apply project memory
- recover safely when recovery is justified
- stop early when the science looks wrong
