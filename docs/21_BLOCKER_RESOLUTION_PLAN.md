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

**Third update, 2026-09-17 evening — teammate integration + C1 workaround.**

| Ledger change | Detail |
|---|---|
| **B4** → reverted in part | the `select=` allowlist subset in `4ec3493` is wrong-axis and redundant; `rm(beta); gc()` kept |
| **B10** → **PARTLY RESOLVED** | `06cd86c` allowlist/titles reworked — 4 defects fixed; label now derives from the stamped provenance class |
| **B13–B16 added** | process gaps exposed by the integration: no load-path test, no app-governance test, `renv.lock` drift, sentinel representativeness |
| **B13, B14** → **RESOLVED** | `tests/test_load_path.R` (7/7) and `tests/test_app_governance.R` (10/10) |
| **B15** → partly resolved | `matrixStats`/`digest` added to `setup.R`; `renv.lock` snapshot still owed |
| **C1b** → strengthened | EB shrinkage **inverts the k=3 verdict**: harmful → +40% of oracle gain |
| **C1c** → not supported by the cheap screen | percentile remap adds zero ranking information and *increases* tissue R² |
| C1 sentinel | array `323195423`, four folds, `priority` queue, gates pre-registered |

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

> ### CORRECTION 2026-09-17 PM — the `select=` half of commit `4ec3493` was reverted
>
> Commit `4ec3493` ("B4: complete memory optimization") added an allowlist
> subset during `fread()`:
>
> ```r
> allowlist_probes <- c("probe_id", readLines("config/shared_autosomal_probes.txt"))
> beta <- data.table::fread(args[1], ..., select = allowlist_probes)
> ```
>
> This was **reverted**, because it is wrong on two independent counts. The
> resolution text above ("subset to the allowlist *during* read") is what invited
> it, so that instruction is now retracted for this matrix layout.
>
> **1. Wrong axis.** `select=` chooses **columns** by name. `beta.tsv` stores
> probes in **rows** and sample barcodes in columns — the header is
> `probe_id  TCGA-OR-A5J1-01A-…  TCGA-OR-A5J2-01A-…`. Passing probe IDs to
> `select=` matches nothing. Reproduced on a fixture with the same orientation:
>
> ```
> WARNING: Column name 'cg00000029' not found in column name header
>          (case sensitive), skipping.
> ```
>
> This is a **warning, not an error**. `fread()` returns a single `probe_id`
> column, and the next line, `t(as.matrix(beta[,-1,drop=FALSE]))`, transposes a
> zero-column frame. The failure is silent and produces wrong-shaped data rather
> than stopping — the most dangerous class of bug in this repository.
>
> **2. Redundant even if the axis were right.** `prepare_beta.py` already
> intersected the matrix with the bridge at build time.
> `data/processed/beta.provenance.json` records
> `probe_allowlist_sha256 = f359e43b…` and `n_probes = 336480`, and the file on
> disk has 336,480 rows against the allowlist's 384,640. The filtering this
> commit intended had already happened upstream; there was nothing left to
> remove, so it could not have saved memory even in principle.
>
> **What was kept.** The `rm(beta); gc()` in the same commit is a genuine
> improvement — it frees the ~21 GB pre-transpose copy — and remains in place.
> A comment block at the call site records why `select=` must not come back.
>
> **Why this was not caught by a test.** `train_baseline.R` has no unit test that
> exercises the load path against a realistically-oriented fixture. See **B13**.

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

### B10. Shiny app governance gap — **PARTLY RESOLVED 2026-09-17**

`app/app.R` validates the schema of `KIDS26_DEMO_RESULTS` but enforces no path
allowlist, checksum, or approval token; raw sample identifiers render as-is; the
`results` object is global and shared across sessions.

**Resolution.** Add a path allowlist + checksum check, hash or alias displayed
identifiers, move `results` inside `server()` for per-session isolation, and make
the provenance label derive from the stamped provenance class rather than from
`nzchar(Sys.getenv(...))`.

> **Two of four items are now done**, building on commit `06cd86c` (Tarun), which
> introduced the path allowlist and provenance-aware titles. That commit had four
> defects, all fixed rather than reverted — the design intent was right.
>
> | # | Defect in `06cd86c` | Fix |
> |---|---|---|
> | 1 | Allowlist contained only `results/baseline/locked_predictions*.tsv`, but `README.md` documents `results/locked_predictions.tsv`. **The gate rejected the only workflow the repo tells you to run.** | Documented path added to the allowlist |
> | 2 | Raw string comparison, so `./results/…` or an absolute path to the same file was rejected | Compare `normalizePath()` on both sides; also stops a symlink pointing outside the allowlist |
> | 3 | Commit message promised "checksum verification"; no checksum was computed | Implemented against the sidecar `sha256` **when `digest` is available**, and warns explicitly when it is not, rather than silently skipping |
> | 4 | Read a `.provenance.json` sidecar **next to the predictions table**. `predict_frozen.R` does not write one — B7 stamps a `provenance` **column**. Requiring the sidecar made every real run fail, and the two `renderUI` titles were never added to the UI, so they rendered nowhere | Gate now reads the **stamped column**, which is the artefact B7 actually guarantees; a single shared `provenance_label` replaces both dead blocks and **is** wired in via `uiOutput("results_title")` |
>
> **Item 4 of the original resolution is now satisfied**: the displayed label
> derives from the stamped provenance class, not from `nzchar(Sys.getenv(...))`.
> A table stamped `NOT A SCIENTIFIC RESULT;OVERRIDDEN_BY_ALLOW_FIXTURE` by
> `--allow-fixture` is now **refused outright** rather than rendered under an
> "Approved precomputed results" heading — which closes the laundering path B7
> identified from the inference side.
>
> Verified executably — the gate can both accept and fail, unlike the B5 mistake:
>
> | Case | Result |
> |---|---|
> | Path outside the allowlist | REJECTED |
> | Allowlisted path, fixture/override stamp | REJECTED |
> | Allowlisted path, clean provenance | ACCEPTED |
> | `./results/…` equivalent spelling | ACCEPTED (was bug 2) |
> | No env var (default fixture) | Loads unchanged |
>
> **Still open in B10:** displayed identifiers are still raw (no hashing or
> aliasing), and `results` is still a global rather than per-session inside
> `server()`. Both matter only once real patient-derived output is displayed, so
> B10 stays **OPEN at P2** rather than closing. These checks are also not yet in
> an automated test — see **B13**.


### B11. PBTP EPIC generation unconfirmed

The exact EPIC generation is unverified, so the 450K/EPIC bridge cannot be
finalized for pediatric transfer. Confirm generation with the data custodian
before freezing; if EPIC v2, the shared-probe intersection must be rebuilt and
its SHA-256 re-recorded.

### B12. Raw IDAT preprocessing bridge unbuilt

Deferred to post-hackathon per `docs/09_POST_HACKATHON_PLAN.md`. Needed only if
raw IDATs enter the pipeline; currently level-3 betas are used.

---

## New blockers opened 2026-09-17 PM

These were found while integrating the teammate branch and running the C1
sentinel folds. They are recorded here rather than fixed silently, because each
one is a *process* gap that will otherwise recur.

### B13. No test exercises the matrix load path — **RESOLVED 2026-09-17**

**Symptom.** The `select=` defect in B4 shipped, was merged, and survived review.
Nothing in the test suite would have caught it: `tests/smoke_model.R` tests
`R/model.R` functions against in-memory matrices, and `tests/test_core.py` tests
the Python acquisition layer. **No test reads a TSV from disk in the orientation
`train_baseline.R` actually expects.**

**Why it matters.** This is the second silent-wrong-shape risk in this file (the
first being the `storage.mode` coercion warning already noted at line 93). A load
bug does not throw — it produces a matrix of the wrong shape that flows onward
and looks like data.

**Resolution.**
1. Add a fixture TSV: a handful of probe rows × a handful of barcode columns,
   with the real header shape (`probe_id` then TCGA barcodes).
2. Assert post-transpose that `nrow(x) == n_samples`, `ncol(x) == n_probes`,
   `colnames(x)` are probe IDs, and `rownames(x)` are barcodes.
3. Assert the load **errors** — not warns — if requested probes are absent.
4. Run it in CI alongside `tests/smoke_model.R`.

> **Done — `tests/test_load_path.R`, 7/7 passing.**
>
> Case 4 is the regression guard and it reproduces the exact B4 failure on a
> same-orientation fixture. One correction to the original diagnosis is worth
> recording, because it makes the bug *worse* than first described: the
> degenerate transpose is **0 × n_probes**, i.e. **zero samples with the probe
> count intact** — not a zero-probe matrix. A downstream check that only
> validated `ncol(x)` against the expected probe count would therefore have
> passed while holding no data at all.
>
> Item 3 was adjusted. The load path does not itself error on absent probes and
> should not be made to, because `apply_preprocess()` aligns by feature name and
> tolerating unseen probes is required for cross-platform transfer (B7 records
> the same reasoning). Case 6 instead asserts that absent probes are
> **detectable by set difference** rather than silently dropped, and case 5
> shows the correct row-wise subsetting idiom so the guard is not misread as
> "probe subsetting is impossible".
>
> Case 7 adds a source-level check that no production beta load reintroduces
> `fread(select=)`, covering `train_baseline.R`, `loco_one_fold.R` and
> `c1_rank_one_fold.R`. That is what makes this durable against a future merge.

---

### B14. Shiny governance checks are untested — **RESOLVED 2026-09-17**

**Symptom.** The B10 allowlist/provenance gates were verified interactively (see
the table under B10) but that verification lives in this document, not in a test
file. It will rot.

**Resolution.** Add `tests/test_app_governance.R` in the style of
`tests/test_provenance_gate.R`, asserting each row of that table executably,
including the negative controls. Bundle with B13 into one CI entry point.

> **Done — `tests/test_app_governance.R`, 10/10 passing.**
>
> Every REJECT case is paired with an ACCEPT case differing in exactly one
> property, so the suite cannot degenerate into the B5 failure of a gate that
> can only ever fire one way. Cases 6 and 7 are load-bearing: a table stamped
> `NOT A SCIENTIFIC RESULT;OVERRIDDEN_BY_ALLOW_FIXTURE` (precisely what
> `predict_frozen.R` writes under `--allow-fixture`) is refused rather than
> rendered under an approved heading, which is what actually closes the
> laundering path B7 identified from the inference side.
>
> Cases 4, 5 and 10 are regression guards tied to specific defects in commit
> `06cd86c`: the documented README path must stay in the allowlist, an
> equivalent path spelling must still resolve, and the displayed label must
> derive from the stamped provenance class rather than from
> `nzchar(Sys.getenv(...))`.
>
> Accepts are checked by **which table loaded**, not merely by the absence of an
> error, so a silent fallback to the synthetic fixture cannot masquerade as a
> pass.

---

### B15. `renv.lock` does not cover the packages now in use — **PARTLY RESOLVED 2026-09-17**

**Symptom.** `scripts/setup.R` installs `glmnet`, `data.table`, `jsonlite`,
`shiny`. Code on `main` now also uses **`matrixStats`** (B3's speedup, load-
bearing for the 22 h → tractable result) and, after the B10 fix, optionally
**`digest`** for checksum verification. `app/app.R` calls `jsonlite::` without
`library(jsonlite)`.

**Why it matters.** B2 closed on the premise that the environment is
reproducible. A fresh `renv::restore()` on another machine can currently produce
an environment where the array silently falls back or the app fails at runtime.

**Resolution.**
1. Add `matrixStats` and `digest` to `scripts/setup.R`.
2. Re-run `renv::snapshot(prompt = FALSE)` after `tests/smoke_model.R` passes and
   commit the updated `renv.lock`.
3. Add an explicit `requireNamespace()` guard wherever a package is used but not
   attached (done for `jsonlite` in `app/app.R`; audit the rest).

> **Items 1 and 3 done.** `scripts/setup.R` now installs `matrixStats` and
> `digest`, with a comment recording why each is load-bearing and how its
> absence fails — `matrixStats` errors outright, while a missing `digest`
> degrades a security check to a warning, which is the quieter and more
> dangerous of the two. `app/app.R` guards `jsonlite` with `requireNamespace()`.
> The setup script's closing comment now lists all four R suites, not just
> `smoke_model.R`.
>
> **Item 2 is NOT done and this blocker stays open because of it.** Running
> `renv::snapshot()` would rewrite `renv.lock` from whatever is installed on the
> current interactive node, which is not verified to be the environment the
> cluster jobs used. Snapshotting from an unvalidated session would replace a
> known-incomplete lockfile with a confidently-wrong one. The snapshot should be
> taken on the analysis host after a clean `renv::restore()` and a full suite
> run.

**Done when.** A clean `renv::restore()` followed by the full test suite passes
on a machine that has never run this project.

---

### B16. Sentinel/preview runs are not representative by construction — **OPEN, process**

**Symptom.** The C2 four-fold log1p preview pointed the *opposite way* from the
full 30-fold run, because the four folds chosen happened to include both tissues
where log1p helps most (THCA, PCPG). The preview suggested correlation was
near-neutral; the full run showed it falling in 20 of 30 tissues.

**Why it matters.** The C1 rank sentinel (array `323195423`) is the same shape of
experiment and carries the same risk. It is mitigated — the four folds were
chosen *a priori* to span the failure modes, and THCA is designated a **negative
control** rather than a win condition — but mitigation is not immunity.

**Resolution.** Any sentinel result must state, before the numbers are read:
which folds were chosen, why, and what result would falsify the hypothesis. A
sentinel may **stop** a full run; it may not on its own **authorise** one unless
the pre-registered gates in `docs/26_C1_N_OF_1_WORKAROUND.md` §7 are met.

**Done when.** The C1 sentinel is reported against its pre-registered gates, with
THCA explicitly excluded from the success criteria.

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

    C1 --> C1B["C1b few-shot calibration<br/>STRENGTHENED - EB shrinkage<br/>makes k=3 helpful (+40% oracle)<br/>but still needs labels"]
    C1 --> C1R["C1c relative-target model<br/>SENTINEL RUNNING - array 323195423<br/>percentile remap FAILED screen<br/>tissue R2 rose 0.562 to 0.593"]

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

> ### C1 WORKAROUND INVESTIGATION 2026-09-17 PM — `docs/26_C1_N_OF_1_WORKAROUND.md`
>
> Three follow-ups were run without touching PBTP or the CNS lock.
>
> **C1b upgraded — shrinkage makes small k safe.** Empirical-Bayes shrinkage of
> the offset toward a cross-tissue prior (τ² = 22.3 against σ² = 108.6, so
> w ≈ k/(k+4.87)) **inverts the k=3 verdict recorded above**. The unshrunk mean
> is harmful at k=3 (−22% of oracle gain); shrunk it recovers **+40%**. At k=1 it
> converts a −214% catastrophe into +22%. Against a *fair hierarchical null*
> given the same prior and the same k labels, the model wins in 25–26 of 29
> tissues at every k.
>
> Honest limits: only 14/29 tissues are helped at k=3, and the entire benefit
> sits in the six tissues with |offset| ≥ 5 (6/6 helped) versus 0/11 for
> |offset| < 2. Uniform shrinkage also **over-shrinks** genuinely large offsets,
> where the unshrunk mean still wins. It remains a **labelled-panel** method and
> does not touch the zero-label case.
>
> **C1c is NOT supported by the cheap screen, and this is the important negative
> result.** Mapping the raw prediction to a within-tissue percentile with a
> frozen source-only monotone map beats its null on MAE (0.226 vs 0.245) and
> gives AUROC 0.729 — but it adds **exactly zero** ranking information (macro
> Spearman identical to the raw prediction to four decimals; a monotone map
> cannot reorder anything), is badly miscalibrated within tissue (macro slope
> 2.61), and its output is **more** explained by tissue identity (R² = 0.593)
> than the raw prediction it came from (0.562). A monotone remap of a
> tissue-confounded score is still a tissue-confounded score.
>
> **The honest test of C1c is the relative-target model**, which learns from a
> within-tissue target instead of remapping afterwards. Implemented in
> `R/rank_model.R` with 12 leakage/invariance tests passing; validated on
> synthetic data only (transfers at ρ ≈ 0.96, manufactures nothing under a pure
> tissue intercept or a noise negative control). **Sentinel array `323195423`
> (BRCA, KICH, THCA, UCEC) submitted to the `priority` queue 2026-09-17 19:37.**
> Gates are pre-registered in `docs/26_C1_N_OF_1_WORKAROUND.md` §7; see **B16**.


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
