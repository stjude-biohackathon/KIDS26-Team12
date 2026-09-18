# 21. Blocker Resolution Plan

Status date: **2026-09-17** (was 2026-09-16). This document consolidates every
known blocker into an ordered, owner-assignable plan. Blockers are grouped by
whether they gate the *first real result* (P0), gate *credible claims* (P1),
gate *model quality* (P1.5, new), or gate *deployment and handoff* (P2).

**CHANGED 2026-09-17 — the headline fact is no longer true.** The previous
version of this document stated: *"no biological model has been fit yet. Every
performance statement in this repository is currently a plan, not a finding."*

That is now superseded. LOCO run 01 (array `323078995`) completed all 30 folds
and merged on 2026-09-17 at 06:20. **A real result exists** and is recorded in
`docs/24_RESULTS_LOCO_RUN01.md`.

| Ledger change | Detail |
|---|---|
| B4 → **CLOSED** | 173 GB measured peak vs 240 GB reserved |
| B3 → **CLOSED** | 30/30 folds in 6.7 h wall |
| B7 → **RESOLVED** | provenance gate, commit `696e05e` |
| **C1–C4 added** | new P1.5 section — model-quality defects found by run 01 |
| Critical path | rewritten; old path complete |

**Second update, 2026-09-17 PM** (`docs/25`). Two of the three investigated
defects inverted their priority:

| Defect | Was | Now |
|---|---|---|
| **C3** purity | most serious | **DOWNGRADED** — artefact of C1; partial correlation *rose* 0.612 → 0.621 when controlling for purity |
| **C1** calibration | tractable fix | **ESCALATED** — label-free approach fails LOTO in every configuration (best R² = −0.116) |
| **C2** zero floor | cheap win | clip **ADOPTED** (+0.008 skill, ranking preserved); log1p full run in flight |
| **C1b** few-shot | untested | **VIABLE** — k=10 labels recover 58% of achievable gain |

Performance statements about the *locked CNS partition* and about *pediatric
transfer* remain plans, not findings. The lock is still closed.

---

## P0 — Gates the first real result

### B1. Historical matrix payload is unlocated — **RESOLVED 2026-09-16**

**Symptom.** `docs/14_VERIFICATION_STATUS.md` records that only a 280,257-byte
header fixture is present locally, while the expected PanCanAtlas payload is
41,541,692,788 bytes. `config/published_beta_metadata.tsv` carries the expected
md5 (`a92f50490cf4eca98b0d19e10927de9d`).

**Why it blocks.** No feature matrix can be extracted, so no model can be fit.

**Resolution.**
1. Locate the payload on shared scratch; record its absolute path in
   `docs/03_DATA_ACQUISITION.md`.
2. Verify with `md5sum` against the manifest value. A mismatch means re-download,
   not "proceed with caution."
3. Record path + checksum + verification date in `CHANGELOG.md`.
4. Only then run `scripts/prepare_beta.py`.

**Done when.** A checksum-verified path is recorded and `prepare_beta.py` emits a
`.provenance.json` marked as real (not fixture).

> ### RESOLVED 2026-09-16
>
> All four conditions are met. Verified state:
>
> | Item | Value |
> |---|---|
> | Raw payload | `data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv` |
> | Size on disk | 41,541,692,788 bytes — **exact match** to the manifest |
> | Expected md5 | `a92f50490cf4eca98b0d19e10927de9d` (`config/published_beta_metadata.tsv`) |
> | sha256 sidecar | `1212f48e8f090fd6afad747e625adf7e66c10a337788aeb81a8a8eaac8f962c4` |
> | Extracted matrix | `data/processed/beta.tsv`, 336,480 probes x 7,707 samples |
> | Provenance | `engineering_only: false`, `track: historical-publication` |
> | Allowlist hash | `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc` |
>
> **How the md5 was verified.** `scripts/acquire_tcga.py` refuses to finalise a
> download unless size AND md5 match the manifest (line 46), and writes the
> `.sha256` sidecar only after that check passes (line 50). The sidecar's
> existence is therefore evidence the md5 gate was satisfied at download time,
> not an independent claim. To re-verify from scratch:
>
> ```bash
> cd data/raw/pancanatlas && sha256sum -c \
>   jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.sha256
> ```
>
> **Note on width.** The extracted matrix has 336,480 probes, not the 384,640 of
> the HM450/EPIC bridge. The bridge is the *allowlist*; 336,480 is its
> intersection with probes actually present in the published matrix. Quote
> 336,480 as the modelling feature count and 384,640 as the platform bridge.
>
> **Still open, and distinct from B1.** The matrix is verified, but per-specimen
> array QC (detection p-values, sex concordance, duplicate audit) has NOT gated
> this cohort — every row carries the same blanket `quality_annotation`. The QC
> gate in `train_baseline.R` now warns about exactly this (see B5).

---

### B2. R runtime never executed — **RESOLVED 2026-09-16**

**Symptom.** `docs/19_R_AND_SHINY_AUDIT.md` states no R executable was available;
all leakage and preprocessing guarantees are verified *by code inspection only*.

**Why it blocks.** Static inspection cannot catch runtime errors, and the
leakage guarantees are the scientific foundation of the whole design.

**Resolution.**
1. Run `Rscript scripts/setup.R` on the cluster; commit `renv.lock` only after a
   successful restore.
2. Run `Rscript tests/smoke_model.R` and record pass/fail in
   `docs/14_VERIFICATION_STATUS.md`.
3. Extend the smoke test with an explicit **leakage assertion**: fit
   preprocessing on a train split, confirm that held-out medians, variance
   ranking, and centering constants are byte-identical to the frozen training
   values and are *not* recomputed on test rows.

**Done when.** Smoke tests pass on real infrastructure and the leakage assertion
is part of the test suite rather than a comment.

> ### SINGLE-FOLD DRY RUN 2026-09-16 (LSF job 323075937, noderome117)
>
> The full modelling path executed end to end on real data for the first time.
> Fold 1 (ACC, the smallest development type) completed successfully:
>
> | Quantity | Value |
> |---|---|
> | Train / test | 6,988 / 77 (locked CNS correctly excluded) |
> | `fit_en` | 53.2 min |
> | Wall time | 55.4 min |
> | **Peak memory** | **173 GB** |
> | Selected hyperparameters | alpha 0.1, lambda 1.411, **interior** to its path |
> | Non-zero coefficients | 772 of 5,000 features |
> | Artefacts written | all five per-fold files |
>
> **Two findings that changed the array configuration.**
>
> 1. **Peak memory is 173 GB, not the 75 GB** the preprocessing-only benchmark
>    implied. `fit_en` additionally holds standardised inner-fold matrices and
>    glmnet's working copies on top of the 19 GB beta matrix. The array was
>    revised from `mem=110GB %15` to `mem=230GB %6`. The original setting would
>    have needed ~3.5 TB against 2.4 TB available and caused mass memory kills.
> 2. **A fold costs ~55 min, not ~45 min**, because the benchmark excluded
>    glmnet. Thirty folds at 6-way concurrency is five waves, roughly 5–6 h.
>
> **The lambda anchoring works as intended.** The selected lambda was interior
> to its data-derived path, which is the outcome that indicates the path was
> wide enough. Under the old hardcoded grid the nearest values were 1 and 10.
>
> **Still outstanding for B2:** `tests/smoke_model.R` has not been run on the
> cluster, and the explicit leakage assertion is not yet in the test suite. The
> dry run exercises the same code path but is not a substitute for that test.

> ### SMOKE TESTS + LEAKAGE ASSERTIONS 2026-09-16 (LSF job 323078428, noderome117)
>
> **B2 is now fully closed.** `tests/smoke_model.R` passed on cluster
> infrastructure, and the leakage guarantee is executable rather than a comment.
>
> Check 7 was added to the suite with six assertions:
>
> | Assertion | What it guards |
> |---|---|
> | 7a | Stored centre/scale equal training-only statistics |
> | 7b | Refitting on the same training rows is deterministic |
> | 7c | Held-out data standardised with training constants stays off-centre |
> | 7d | The transform object is not mutated by being applied |
> | 7e | Features align by NAME with unseen probes present and order shuffled |
> | 7f | The lambda path is reproducible from training rows alone |
>
> **7c is the load-bearing one, and it was validated against a negative
> control.** A deliberately leaky implementation that recentres on the data it
> is handed produces column means of ~7e-17, while the correct implementation
> produces ~1.2. The assertion threshold (>0.5) separates them decisively, so
> this is a test that can actually fail — unlike the tautological QC gate
> described in B5.
>
> One adjustment during development: 7a uses `all.equal(tolerance=1e-12)` rather
> than `identical()`, because `fit_preprocess()` uses `matrixStats::colMeans2`
> while the check uses base `colMeans`. They agree to ~1e-16 — a summation-order
> difference, not a difference in what was computed. Using `identical()` would
> have made the test fail for a reason unrelated to leakage.

---

### B3. Preprocessing is ~42 hours for full nested LOCO — **MITIGATED 2026-09-16**

**Symptom.** `R/model.R` lines 46–62: `fit_preprocess()` uses `apply()` across
~336k columns, ~2.8 minutes per call; nested LOCO needs ~900 calls.

**Resolution.**
1. Replace column-wise `apply()` with `matrixStats::colMedians` / `colVars`
   (documented in-file as reducing the estimate to ~17 hours).
2. Restrict the missingness filter and variance ranking to a single pass.
3. Benchmark on one outer fold before launching the full grid.
4. Submit as an LSF array job, one outer fold per task, rather than a single
   long-running job.

**Caution.** Do not reduce cost by screening features once on the full cohort —
that reintroduces exactly the leakage the design prevents. Speed fixes must be
implementation-level only.

> ### MEASURED 2026-09-16 (LSF job 323075601, nodegpu217)
>
> `matrixStats` reductions are in place and verified numerically identical to
> the old `apply()` path on Inf-free input. Benchmarked on the real
> 336,480 x 7,707 matrix rather than extrapolated:
>
> | Quantity | Value |
> |---|---|
> | `fit_preprocess` | 72.8 s median (73.1 / 71.5 / 72.8) |
> | `apply_preprocess` | 8.7 s |
> | Calls per outer fold | 30 x `fit_preprocess`, 58 x `apply_preprocess` |
> | Per outer fold | ~44.8 min |
> | Serial Phase 1 (30 folds) | **~22.4 h** |
> | Peak RSS | 75 GB |
>
> **An earlier extrapolation from 20k/60k-probe test matrices predicted 19.9 s
> per call and 8.7 h serial. It was wrong by 3.7x.** Cost does not scale
> linearly with probe count once the matrix is ~19 GB, because memory bandwidth
> rather than arithmetic becomes the limit. Treat small-matrix extrapolations in
> this project as lower bounds only.
>
> **Mitigation, not elimination.** Item 4 is what actually makes this tractable:
> `scripts/lsf_loco_array.bsub` runs the 30 folds as independent array tasks at
> 15-way concurrency, giving ~1.5–2 h wall time in two waves. A single
> `fit_preprocess` call is still ~73 s and that has not changed.

**Done when.** One outer fold completes in a measured, recorded wall time and the
full run is scheduled within the event window.

---

### B4. Memory ceiling on matrix load — **CLOSED 2026-09-17**

**Symptom.** `scripts/train_baseline.R` lines 33–72: ~21 GB for `fread()`, 60+ GB
peak through the transpose chain.

**Resolution.** Request a large-memory queue; subset to the 384,640-probe
allowlist *during* read rather than after; avoid retaining intermediate copies;
consider `data.table::setDT` in-place transforms. Record peak RSS.

> **Closed by measurement.** All 30 array tasks completed with **173 GB peak RSS**
> against a 240 GB reservation (`-n 4` × `rusage[mem=60GB]` per slot), a 28%
> headroom margin. Zero OOM kills, zero EXIT'd tasks across the full array.
> The earlier 75 GB estimate in the schematic was low; 173 GB is the measured
> figure and is what future reservations should be sized against.
>
> Note this is the *in-loop* peak, not the `train_baseline.R` serial path, which
> remains unmeasured but is superseded by the array for Phase 1 work.

---

## P1 — Gates credible claims

### B5. QC gate is tautological — **RESOLVED 2026-09-16**

**Symptom.** `scripts/train_baseline.R` line 128 checks
`quality_annotation == "published_450K_no_exclusion"`, but `build_master.py`
hardcodes that exact value. The check can never fail.

**Resolution.** Replace with a real gate: assert the provenance sidecar exists,
that its recorded checksum matches the loaded matrix, and that per-sample QC
fields (detection p-value pass rate, sex concordance) are present and within
bounds. A gate that cannot fail is worse than no gate, because it produces false
assurance in the audit trail.

---

### B6. Latent `sample()` bug — **FIXED 2026-09-16**

**Symptom.** `scripts/train_baseline.R` lines 212–217: `sample(ii, k)` where
`ii` has length 1 silently samples from `1:ii` instead of returning `ii`.

**Why it matters.** Currently masked (smallest development group, OV, has n=10)
but will fire on any subset or rarer stratum, producing invalid fold assignment
without an error.

**Resolution.** Guard with `if (length(ii) == 1L) ii else sample(ii, k)`, or use
`resample <- function(x, ...) x[sample.int(length(x), ...)]`. Add a unit test
with a single-element group.

---

### B7. `predict_frozen.R` has no provenance gate — **RESOLVED 2026-09-16**

**Symptom.** Unlike the training script, inference accepts any beta matrix and
will happily score an engineering fixture, then pass it to the Shiny app via
`KIDS26_DEMO_RESULTS` where it renders as "Approved precomputed results."

**Resolution.** Read the `.provenance.json` sidecar; refuse to score matrices
marked as simulated/fixture unless an explicit `--allow-fixture` flag is passed,
and stamp the output file with the provenance class so downstream display cannot
misrepresent it.

> **Done.** The gate now lives in `R/provenance.R` — shared rather than copied,
> so training and inference cannot drift apart — and runs in
> `scripts/predict_frozen.R` *before* the matrix is read, so an unlabelled 28 GB
> file is refused without being loaded.
>
> Three classes are refused, not one: a declared fixture
> (`engineering_only=true`), a matrix with **no sidecar at all**, and a
> non-fixture matrix whose `probe_allowlist_sha256` is absent. The middle case
> matters most — the original gap was reachable by simply dropping an arbitrary
> TSV on disk, which carries no `engineering_only` flag to catch.
>
> `--allow-fixture` permits engineering runs but cannot launder the label: the
> class is suffixed `;OVERRIDDEN_BY_ALLOW_FIXTURE`, a warning is raised,
> `reportable` is forced `FALSE`, estimates and bounds are blanked, and the
> `provenance` column leads with `NOT A SCIENTIFIC RESULT`. Since `app/app.R`
> validates `KIDS26_DEMO_RESULTS` structurally only and prints `provenance`
> verbatim, putting the class *in that column* is what actually closes the path
> to the demo.
>
> Probe-allowlist mismatch between model and scoring matrix is recorded per row
> (`probe_set_match`) and warned about, but is **not** fatal — cross-platform
> transfer to pediatric EPIC arrays is the research goal, and
> `apply_preprocess()` aligns by feature name.
>
> `tests/test_provenance_gate.R` asserts each refusal executably, including
> against the real `data/processed/beta.provenance.json`.

---

### B8. Lambda grid is hardcoded — **RESOLVED 2026-09-16**

**Symptom.** `R/model.R` uses five fixed lambda values spanning four decades, not
anchored to a glmnet-derived `lambda.max`. A boundary optimum is accepted
silently.

**Resolution.** Derive the path from `glmnet`'s own `lambda.max` per inner fold;
warn loudly if the selected lambda sits at either endpoint of the grid.

> **Done.** `lambda_path()` in `R/model.R` computes
> `lambda.max = max|x'y| / (n * alpha)` per alpha from the fold's own training
> rows, log-spaced down to `0.001 * lambda.max`. Verified to reproduce glmnet's
> internal `lambda.max` exactly at alpha 0.1 / 0.5 / 1.
>
> The old fixed grid is confirmed to have been wasteful: at alpha = 1,
> `lambda.max` is ~2.77, so the hardcoded values 10 and 100 were both guaranteed
> intercept-only fits — two of five grid points were dead.
>
> `check_lambda_boundary()` warns when the winner lands on an endpoint, and the
> fitted bundle records `lambda_paths`, `lambda_at_boundary` and
> `lambda_boundary_side` so the condition survives into the audit trail.
> `scripts/loco_merge.R` reports how many folds hit a boundary and warns if a
> majority did.

---

### B9. `HRD_high_probability` is always NA — **RESOLVED 2026-09-16**

**Symptom.** The exploratory threshold (42) lives in
`config/analysis_protocol.json` but is read by no R code.

**Resolution.** Either (a) drop the column until a calibration set justifies it,
or (b) wire it explicitly and label it exploratory in every output. Option (a) is
preferred before any public presentation — an always-NA column invites a reader
to assume a clinical threshold exists.

> **Done — option (a).** The column was removed from `scripts/predict_frozen.R`
> and replaced with a comment stating why, so the absence is deliberate and
> documented rather than an oversight for someone to "fix" by re-adding it.

---

## P2 — Gates deployment and handoff

### B10. Shiny app governance gap

`app/app.R` validates the schema of `KIDS26_DEMO_RESULTS` but enforces no path
allowlist, checksum, or approval token; raw sample identifiers render as-is; the
`results` object is global and shared across sessions.

**Resolution.** Add a path allowlist + checksum check, hash or alias displayed
identifiers, move `results` inside `server()` for per-session isolation, and make
the provenance label derive from the stamped provenance class rather than from
`nzchar(Sys.getenv(...))`.

### B11. PBTP EPIC generation unconfirmed

The exact EPIC generation is unverified, so the 450K/EPIC bridge cannot be
finalized for pediatric transfer. Confirm generation with the data custodian
before freezing; if EPIC v2, the shared-probe intersection must be rebuilt and
its SHA-256 re-recorded.

### B12. Raw IDAT preprocessing bridge unbuilt

Deferred to post-hackathon per `docs/09_POST_HACKATHON_PLAN.md`. Needed only if
raw IDATs enter the pipeline; currently level-3 betas are used.

---

## Suggested execution order

Status as of 2026-09-17: **B1, B2, B5, B6, B7, B8 and B9 are RESOLVED; B3 and
B4 are MITIGATED and closed in practice.** The 30-fold LOCO array (`323078995`)
completed and merged; results are in `docs/24_RESULTS_LOCO_RUN01.md`.

The first real result exists. The remaining work is no longer blocker removal —
it is the **four model-quality defects** that run 01 exposed (C1–C4 below),
which gate opening the CNS lock.

```mermaid
graph TD
    B1["B1 locate + checksum matrix<br/>RESOLVED"] --> B4["B4 memory-safe load<br/>CLOSED: 173 GB measured<br/>vs 240 GB reserved"]
    B2["B2 R runtime + smoke tests<br/>RESOLVED on cluster"] --> B3["B3 matrixStats + LSF array<br/>CLOSED: 30/30 folds in 6.7 h"]
    B4 --> B3
    B5["B5 real QC gate<br/>RESOLVED"] --> FIT
    B6["B6 sample bug<br/>FIXED"] --> FIT
    B3 --> FIT["LOCO run 01 COMPLETE<br/>array 323078995, n=7065<br/>30/30 folds, 0 boundary hits"]
    FIT --> CONF["Confounding controls RUN<br/>tissue + purity + permutation"]
    CONF --> VERDICT{"Tissue confound<br/>disqualifying?"}
    VERDICT -->|"NO — within-tissue r=0.61,<br/>perm p=0.001, 29/30 tissues"| CONTINUE["CONTINUE MODELLING<br/>reframed as within-tissue ranker"]

    CONTINUE --> C1["C1 per-tissue calibration<br/>ESCALATED - label-free<br/>approach FAILED, LOTO R2 negative"]
    CONTINUE --> C2["C2 zero floor - RESOLVED<br/>clip ADOPTED, log1p REJECTED<br/>within-tissue r fell 0.612 to 0.522"]
    CONTINUE --> C3["C3 purity inversion<br/>DOWNGRADED - artefact of C1,<br/>partial cor went UP 0.612 to 0.621"]
    CONTINUE --> C4["C4 OV n=10 disclosure<br/>OPEN - 27k array excluded"]

    C1 --> C1B["C1b few-shot calibration<br/>VIABLE - k=10 recovers 58%<br/>but needs 10 labels per new tissue"]
    C1 --> C1R["C1c within-tissue rank only<br/>UNTESTED - needs same-type<br/>reference cohort at predict time"]

    C1B --> LOCK
    C1R --> LOCK
    C2 --> LOCK
    C3 --> LOCK
    C4 --> LOCK
    LOCK["CNS LOCK - still closed<br/>opens ONCE, after C1-C4"]

    B7["B7 inference provenance<br/>RESOLVED 2026-09-16"] --> B10["B10 app governance<br/>OPEN, P2"]
    LOCK --> B10
    B11["B11 PBTP EPIC<br/>OPEN"] --> TRANSFER["Pediatric transfer"]
    B10 --> TRANSFER

    style B1 fill:#d4edda
    style B2 fill:#d4edda
    style B3 fill:#d4edda
    style B4 fill:#d4edda
    style B5 fill:#d4edda
    style B6 fill:#d4edda
    style B7 fill:#d4edda
    style FIT fill:#d4edda,stroke:#0f5132,stroke-width:3px
    style CONF fill:#d4edda
    style CONTINUE fill:#d1e7dd,stroke:#0f5132,stroke-width:3px
    style C1 fill:#f8d7da,stroke:#b02a37,stroke-width:3px
    style C2 fill:#d1e7dd,stroke:#0f5132
    style C3 fill:#fff3cd,stroke:#d39e00
    style C4 fill:#fff3cd,stroke:#d39e00
    style C1B fill:#d1e7dd,stroke:#0f5132,stroke-width:2px
    style C1R fill:#ffe08a,stroke:#d39e00
    style LOCK fill:#f8d7da,stroke:#b02a37,stroke-width:2px
    style B10 fill:#fff3cd
    style B11 fill:#fff3cd
```

**Critical path is now:** fix C1–C4 → re-run the array → confirm skill improves
→ freeze → open the CNS lock **once** → B10 app governance → pediatric transfer.

---

## P1.5 — Model-quality defects exposed by run 01 (NEW 2026-09-17)

These are not blockers in the original sense (nothing is broken or unverifiable);
they are substantive modelling defects measured in `docs/24_RESULTS_LOCO_RUN01.md`.

### C3. Purity inversion — **DOWNGRADED 2026-09-17, was "most serious"**

`cor(pred, purity)` within tissue = 0.165 vs `cor(label, purity)` = 0.019 — the
model tracks tumour purity ~8× more strongly than the truth does. Skill falls
monotonically with purity: +0.190 (low) → +0.056 (mid) → **−0.041 (high)**. In
the cleanest samples the model **loses to the tissue mean**.

A genuine biological signal should get *stronger* with purity, not weaker. Until
this is explained, no strong biological claim is defensible.

> **Investigated and largely exonerated — `docs/25` §C3.** The decisive test is
> the partial correlation: if purity drove the signal, controlling for it would
> collapse the within-tissue correlation. It does the opposite —
> **0.6124 → 0.6210**, a slight *increase*.
>
> The inversion was an artefact of C1. The per-tissue offset is essentially
> constant across purity tertiles (3.50 / 3.55 / 3.49), so it consumes a fixed
> ~3.5 units of a margin that shrinks as the tissue-mean null becomes easier to
> beat at high purity. Remove the offset and skill is **positive in all three
> tertiles** (0.242 / 0.207 / 0.159). `cor(per-tissue offset, mean purity)` is
> −0.124 — weak, and the wrong sign to explain an inversion.
>
> Purity is a nuisance variable the model partially encodes, **not** the source
> of its signal. A residual gradient survives correction (0.242 → 0.159), so
> this is downgraded rather than closed.
>
> **Not tested:** whether individual selected CpGs are purity-associated. That
> needs a probe-level job against the 28 GB matrix (~4 h). Probe selection is
> stable enough to make it worthwhile — 4,080 of 5,000 probes are shared across
> all folds examined.

### C1. Per-tissue calibration offset — **ESCALATED 2026-09-17**

Mean absolute per-tissue offset is 3.46 HRD units (SD 4.83, range −15.44 to
+8.29). Removing it drops pooled MAE 9.049 → 7.881. This single defect costs
more than the model's entire margin over the tissue-mean null (+0.085).

**Resolution.** Fit a per-tissue offset **inside** the LOCO loop, on training
tissues only. Fitting it on the held-out tissue is leakage and voids the fold.

> **The label-free approach FAILED — `docs/25` §C1.** Learning
> `offset ~ tissue covariates` (mean prediction, prediction SD, mean purity,
> mean OOD distance, n) and evaluating leave-one-tissue-out — which mirrors
> deployment exactly — **every configuration lost to predicting the global mean
> offset.** Best LOTO R² = −0.116; the full four-covariate model reached −0.217.
> The offset is uncorrelated with every available covariate (|r| ≤ 0.267).
>
> The negative result is trustworthy: `tests/test_calibration.R` case 8 recovers
> a synthetic linear offset at LOTO R² > 0.8, and case 9 confirms the machinery
> does not manufacture signal from noise. With 30 tissues as 30 observations,
> there is simply not enough information.
>
> **Few-shot calibration DOES work — `docs/25` §C1b.** ~10 labelled samples from
> the new tissue recover **58% of the achievable gain** (MAE 8.43 → 7.86) and
> beat the fair k-label null in 20 of 29 tissues. k=3 actively *hurts* (offset
> error 4.35 vs true offset SD 4.83). The gain is concentrated almost entirely
> in mis-levelled tissues: `cor(gain, |true offset|) = 0.991`, and only 14 of 29
> tissues improve at all.
>
> **Consequence for the project.** Absolute HRDsum on a brand-new tissue with
> zero labels is **not currently achievable**. Two supported paths remain:
> (a) report within-tissue relative rank, which needs a same-type reference
> cohort, or (b) obtain ~10 labelled samples per new tumour type. Neither solves
> the N-of-1 pediatric case. This is an open scientific problem, not a coding
> task.

### C2. Unbounded predictions against a zero-floored label — **PARTLY RESOLVED**

HRDsum ≥ 0 by construction and 14.3% of samples are exactly 0. The elastic net
emits negative predictions, inflating MAE with a priori impossible values.

**Resolution.** Clip at 0, or model `log1p(HRDsum)` and back-transform.

> **Clipping ADOPTED.** No refit needed: MAE 9.049 → 8.972, skill 0.085 →
> 0.093, Spearman with the raw prediction exactly **1.000**, so no ranking
> changes. 227 samples (3.2%) were negative, most extreme −11.53.
>
> **log1p REJECTED** after the full 30-fold run (arrays `323176856` +
> `323169191`, 30/30, merged 2026-09-17). It improves the absolute metrics
> (MAE 9.049 → 8.827, skill 0.085 → 0.108, zero negative predictions) but
> **degrades within-tissue Pearson 0.6124 → 0.5221**, with correlation falling
> in 20 of 30 tissues. The 4-fold preview understated this because it happened
> to contain both tissues where log1p helps most.
>
> The losses land on BLCA, STAD, LUSC and ESCA — high-HRD tissues where
> discrimination is the clinically useful task. Since C1 forces the project
> toward **within-tissue ranking**, and ranking is measured by correlation,
> log1p optimises the metric we must abandon at the cost of the one we must
> rely on. Retained only as a documented option for the quiet-tumour regime
> (THCA −0.032 → 0.171).
>
> Zero inflation is **not** addressed by log1p (it maps 0 → 0). A genuine
> hurdle / two-part model is recorded as future work.

### C4. Ovarian cohort is n = 10 — **OPEN (disclosure, not a code fix)**

OV is the canonical HRD cancer and the main clinical application. TCGA ovarian
methylation is mostly 27k-array, excluded by the 450k bridge. HRD-high behaviour
is therefore inferred from UCEC/BRCA/STAD.

**Resolution.** State this limitation in every presentation. Optionally rebuild
the bridge to include 27k probes, accepting a much smaller probe intersection.

---

**Superseded 2026-09-17.** The old critical path (dry run → array → merge → read
the pooled tissue-mean null) is **complete**. The pooled null was read: skill
+0.085, within-tissue r = 0.612, permutation p = 0.001. The answer to "is there
a result worth reporting at all" is **yes, reframed** — see
`docs/24_RESULTS_LOCO_RUN01.md` §6. The path forward is C1–C4 above.
