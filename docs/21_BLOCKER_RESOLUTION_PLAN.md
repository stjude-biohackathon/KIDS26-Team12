# 21. Blocker Resolution Plan

Status date: **2026-09-18** (was 2026-09-17, was 2026-09-16). This document consolidates every
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
| **C1c** → **REFUTED for V1/V2** | sentinel `323195423` failed 3 of 7 gates; tissue R² **57.6%** vs 52.7% for the absolute model; V3 untested |
| C1 sentinel | array `323195423`, 4/4 complete, gates scored in `docs/26` §7 — **full array NOT authorised** |
| **CNS lock** | **stays closed** — C1 is more firmly unresolved, not less |

**Fourth update, 2026-09-18 — the lock was opened.** Two source-only
experiments ran (one passed, one failed), a candidate was frozen under a
pre-registered rule, the shipping inference path was audited and repaired, and
the locked CNS cohort was scored **exactly once** under LSF job `323264626`.
Every line below is evidenced in `results/` and in `docs/27`, `docs/28`,
`docs/29`.

| Ledger change | Detail |
|---|---|
| **CNS lock** | **OPENED 2026-09-18, CONSUMED, CANNOT BE REUSED** — job `323264626`, 642 samples (135 GBM + 507 LGG), 2 ledger rows in `results/LOCKED_EVALUATION_LEDGER.tsv`, one per candidate |
| **C1** → **TESTED AND CONFIRMED** | no longer "unresolved pending investigation": **characterised**, **mitigated ~42% inside the source distribution** by V3-abs, and **still decisive outside it** on CNS |
| **C1c** → **V3-abs PASSED and SHIPPED** | 7/7 pre-registered sentinel gates (`docs/27` §5) and 7/7 selection criteria S1–S7 (`results/v3_loco_full/selection_scorecard.tsv`); frozen as **Candidate B**, sha256 `df9f7e82…` |
| **C1c** → **V4 lineage-penalised FAILED** | the one pre-registered alternative (`docs/27` §6) failed **G3, the mechanism gate**: tissue R² of predictions rose 0.5219 → **0.5361** against a 0.4919 ceiling. Branch terminated as pre-registered. Passed G1, G2, G4, G5, G6, G7 |
| **B17 added** | the frozen OOD/reportability rule **did not flag CNS as out of distribution** — 8 of 642, all LGG. A defect, not a reassurance |
| **B18 added** | the exploratory `HRDsum ≥ 42` cutoff is **vacuous in CNS**: 0 GBM and exactly 1 LGG sample reach it. No binary claim is possible in this lineage |
| **B19 added → RESOLVED** | commit `adf420b`: five defects on the **shipping inference path**, including one that meant it **could not score anything at all** |
| **B20 added** | the CNS `bsub` exited non-zero on a **script-contract mismatch**, after scoring had already succeeded. Process gap, not a data problem |
| C2 | carried forward — but see **B19** item 1: the adopted clip was never actually applied on the shipping path until `adf420b` |
| C3 | carried forward — purity populated for 632 of 642 CNS rows; residual diagnostic reported, still downgraded, still not closed |
| C4 | carried forward — OV n = 10 is now **joined by a second cohort-composition problem**, **B18** |
| **B10, B11, B12, B15, B16** | **remain OPEN**, unchanged by today |
| **Deployment gate** | **stays shut** (`docs/27` §2A). Nothing today produces a validated clinical HRD assay |

**The scientific headline, stated so it cannot be rounded up.** V3-abs removes
**42.0%** (95% CI 37.9–46.5%) of the excess lineage structure Candidate A
carries beyond the target's own, `ΔR²_excess = 0.0913` [0.0824, 0.0996],
p < 0.0001 (`results/tables/tableA2_A_vs_V3_full30.md`) — and that reduction,
measured on the **source** tissues, **did not translate into better CNS
calibration**. Reducing lineage imprinting inside the training distribution was
**not sufficient** for transfer outside it. That is the most informative single
sentence the project produced today, and it is a negative one.

The pre-declared verdict is **outcome A of the three fixed in `docs/28` §7**:
rank transfers weakly, absolute cross-lineage calibration does not
(`results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md` §5).

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

> **Still OPEN 2026-09-18, and it earned its keep twice today.** B16 is the rule
> that forced `docs/27` and `docs/28` to be written *before* their inputs
> existed. It worked: the V3-abs sentinel passed 7/7 and was then required to
> pass a *separately pre-registered* 30-fold rule (S1–S7) before it could
> displace run 01; the V4 sentinel failed G3 and the branch was terminated
> rather than re-tuned. The blocker stays open because it is a standing process
> rule, not a task — and because the sentinel-vs-full divergence it warns about
> **recurred**: V3-abs's sentinel tissue R² was 0.2991 on four tissues but
> 0.4672 across all thirty. The direction held; the magnitude did not.

---

## New blockers opened 2026-09-18

Found during the pre-unlock audit and the locked CNS evaluation. B1–B16 were
taken; these take the next free numbers.

### B17. The OOD reportability rule did not flag CNS — **OPEN, P1, and this is a defect**

**Symptom.** The frozen `ood_score` / `reportable` flags were carried into the
CNS evaluation exactly as pre-declared (`docs/28` §6). They flagged **8 of 642
samples — 1.2%, all of them LGG, none of them GBM**
(`results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md` §3). GBM
ood_score ranges 0.81 / 1.07 / 1.43 (min/median/max), LGG 0.76 / 1.22 / 1.55.

**Why this is a defect and not a reassurance.** CNS is, on the project's own
primary interpretive frame, **a poor relation rather than a typical held-out
tissue**: GBM Pearson sits at the **13.3rd percentile** of the 30 source LOCO
tissues (28th of 32), LGG at the **16.7th percentile** (26th of 32), and GBM
mean bias at the **10th percentile** (29th of 32). The oracle tissue-mean null
beats the model in every group (skill −1.009 GBM, −0.048 LGG, −0.233 pooled).
A rule whose entire purpose is to warn a user "this sample is unlike anything
the model was trained on" declared **98.8% of a genuinely unseen lineage
reportable**, including the lineage where the model is worst.

The flag is therefore **not a safety net**. It is currently evidence only that
CNS methylation is within the numeric envelope of the source feature
distribution — which is a statement about the features, not about whether the
prediction can be trusted.

**Why it matters for the actual goal.** The pediatric transfer (B11) is a
*further* domain shift on top of this one. If the rule cannot flag adult CNS, a
green "reportable" on a PBTP sample carries no information, and would be read
as assurance by exactly the audience least able to check it.

**Resolution.**
1. Do **not** silently re-tune the threshold — that would be fitting the OOD
   rule to the CNS labels we have now consumed, and is the B5 failure mode in
   a new costume.
2. Characterise what the score actually measures against **source** data only:
   its per-tissue distribution across the 30 LOCO folds, and whether it
   correlates with per-tissue MAE or |bias| at all. If it does not correlate
   with error on held-out source tissues, it is not an error-prediction rule
   and must not be labelled as one.
3. Report reportability alongside a **calibration** statement, not instead of
   one. "Reportable" must not be allowed to mean "calibrated".
4. Until 2 is done, every surface that displays `reportable` (the Shiny app —
   see B10) states that the flag was tested on CNS and **did not fire**.

**Done when.** Either the score is shown to predict held-out-tissue error on
source data, or the flag is relabelled to say only what it can support.

---

### B18. The `HRDsum ≥ 42` threshold is vacuous in CNS — **OPEN, disclosure + design**

**Symptom.** The exploratory cutoff has **0 GBM samples** and **exactly 1 LGG
sample** at or above it
(`results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md` §4). The
confusion matrices are degenerate: GBM 0/0/0/135, LGG 0/0/1/506. GBM's AUC is
`NA` because there is no positive class at all.

**Why it matters.** Any binary "HRD-high" claim in this lineage is
**arithmetically impossible**, not merely underpowered. The reported LGG AUC of
0.966 rests on a **single positive sample** and must never be quoted as
evidence of anything; it is a one-sample order statistic. This is the same
class of problem as **C4** (OV n = 10) — the cohort cannot support the claim —
but sharper, because here the limiting count is one rather than ten.

It also retroactively vindicates **B9**, which removed the always-`NA`
`HRD_high_probability` column: had that column shipped, it would now be
rendering a binary call in a lineage where the threshold separates nothing.

**Resolution.**
1. State in every CNS-facing figure, caption and slide that the 42 cutoff is
   **inapplicable to this lineage**, with the counts (0 and 1) given.
2. Keep the continuous within-lineage ranking as the only reported endpoint
   for CNS. It is the one the data can support.
3. Do not derive a CNS-specific cutoff from these 642 samples. The lock is
   consumed; a threshold fitted on them would be fitted on the test set.

**Done when.** No presented artefact implies a binary HRD call in CNS.

---

### B19. The shipping inference path was broken and unguarded — **RESOLVED 2026-09-18, commit `adf420b`**

**Symptom.** A pre-unlock audit of `scripts/predict_frozen.R` — the script that
was about to consume the one-shot lock — found five defects, of which the
second is the serious one:

| # | Defect | Consequence |
|---|---|---|
| 1 | The adopted **C2 clipping was never applied** on the shipping path | Inference emitted **negative HRDsum**, an a priori impossible value, while `docs/21` C2 recorded clipping as ADOPTED. The doc and the code disagreed |
| 2 | `predict_frozen.R` **never loaded `glmnet`**, so `predict()` had no S3 method | **The shipping inference path could not score anything at all.** Pre-existing, not introduced by the audit; confirmed against stashed pre-change code |
| 3 | No locked-partition guard and **no record of scoring** | The lock could have been consumed silently, or repeatedly, with no artefact proving it |
| 4 | **Two sources of truth** for the lock: a hardcoded `c("GBM","LGG")` in five files versus the `partition` column | The two could drift and nothing would notice |
| 5 | `loco_merge.R`'s completeness gate only globbed `predictions_*.tsv` | A fold could be missing its metrics or nullpanel, or run a *different* transform or feature filter, and still merge |

**Why #2 is worth dwelling on.** It is the third instance in this repository of
the failure mode B13 was opened for: a defect on a path that no test exercised.
B7's provenance gate, B10's governance gate and C2's clipping rule were all
*specified* on this path, and all three were guarding a script that would have
errored before reaching them.

> **Done — `adf420b`, all 7 suites pass, each fix paired with a REJECT and an
> ACCEPT test.**
>
> 1. Clipping now goes through the bundle's own `target_transform` via
>    `get_transform()` — the same object `loco_one_fold.R` back-transforms with,
>    so the two cannot drift — and the pre-clip value is retained in
>    `predicted_reference_HRDsum_raw` so the clip is auditable rather than
>    invisible.
> 2. `glmnet` is loaded and the requirement guarded.
> 3. Scoring `locked_CNS` rows requires an explicit **`--locked-evaluation`**
>    opt-in and appends to an **append-only ledger**,
>    `results/LOCKED_EVALUATION_LEDGER.tsv`, recording commit, artifact sha256,
>    matrix sha256 and output sha256. Repeated scoring becomes *visible* rather
>    than silent — it does not become impossible, and should not be mistaken
>    for a lock that enforces itself.
> 4. `assert_partition_matches_cns()` requires **two-way agreement** between the
>    hardcoded set and the `partition` column, and is called from every fold
>    script.
> 5. `loco_merge.R` now checks metrics / nullpanel / prediction counts agree and
>    that all folds share **one `target_transform` AND one `feature_rank`**, and
>    ships the tissue-identity R² — the lineage-confound number, previously
>    computed ad hoc in one probe script — on **every** run.
>
> The ledger was used as intended: both rows' artifact checksums were verified
> with `sha256sum -c` **before** scoring, against the values recorded in
> `docs/29` §1.

---

### B20. The CNS job exited non-zero on a script-contract mismatch — **OPEN, process**

**Symptom.** The CNS `bsub` (job `323264626`) returned a non-zero exit status.
The cause is a contract mismatch between two scripts, not a data fault:
`predict_frozen.R` scores the **whole 7,707-sample matrix** — which is correct,
it is a scoring script and the partition guard is what restricts what may be
*reported* — while `analyze_cns.R` requires **CNS-only** input.

**What actually happened, stated plainly because the audit trail depends on
it.** Scoring had **already succeeded** and **both ledger rows were written**
(`results/LOCKED_EVALUATION_LEDGER.tsv`, 20:05:08Z and 20:06:47Z, both
`n_samples_scored = 7707`, `n_locked_CNS_scored = 642`) before the failure. The
analysis was then completed by filtering to the 642 locked rows
(`cns_only_v3.tsv`, `cns_only_A.tsv`). **Nothing was rescored. The lock was
consumed exactly once.**

**Why it is recorded anyway.** A non-zero exit on the single most consequential
job in the project is exactly the signal that must never be ambiguous. Had the
ledger not existed (**B19** item 3), there would be no machine-checkable way to
distinguish "failed before scoring" from "failed after scoring" — and the
temptation to re-run would have been strong and unanswerable.

**Resolution.** Make the contract explicit: either `predict_frozen.R` gains an
option to emit a partition-filtered output alongside the full one, or
`analyze_cns.R` filters its own input and says so. Whichever is chosen, the
wrapper must distinguish "scoring failed" from "post-scoring analysis failed"
in its exit status.

**Done when.** A one-shot evaluation job either succeeds end to end, or fails
with an exit status that identifies which stage failed.

---

## Suggested execution order

Status as of 2026-09-17: **B1, B2, B5, B6, B7, B8 and B9 are RESOLVED; B3 and
B4 are MITIGATED and closed in practice.** The 30-fold LOCO array (`323078995`)
completed and merged; results are in `docs/24_RESULTS_LOCO_RUN01.md`.

The first real result exists. The remaining work is no longer blocker removal —
it is the **four model-quality defects** that run 01 exposed (C1–C4 below),
which gate opening the CNS lock.

**Revised 2026-09-18.** The gate above was **wrong**, and `docs/27` §2 explains
why: conditioning the CNS unlock on *resolving* C1 is circular, because the CNS
evaluation is the only experiment capable of measuring zero-shot calibration in
an untouched lineage. C1 is a **finding to be tested on CNS**, not a
precondition for testing. The unlock was re-gated on **Gate B** (freeze
discipline: pre-registered candidate, checksummed artifact, predeclared
metrics, auditable pre-unlock record) and executed on 2026-09-18.

**The lock is now consumed.** The critical path below is retained as a
historical record of what was believed on 2026-09-17; the current path is
stated after the diagram.

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

    CONTINUE --> C1["C1 per-tissue calibration<br/>TESTED ON CNS - CONFIRMED<br/>dominant failure mode"]
    CONTINUE --> C2["C2 zero floor - RESOLVED<br/>clip ADOPTED, log1p REJECTED<br/>within-tissue r fell 0.612 to 0.522"]
    CONTINUE --> C3["C3 purity inversion<br/>DOWNGRADED - artefact of C1,<br/>partial cor went UP 0.612 to 0.621"]
    CONTINUE --> C4["C4 OV n=10 disclosure<br/>OPEN - 27k array excluded"]

    C1 --> C1B["C1b few-shot calibration<br/>STRENGTHENED - EB shrinkage<br/>makes k=3 helpful (+40% oracle)<br/>but still needs labels"]
    C1 --> C1R["C1c relative-target model<br/>V1/V2 REFUTED - sentinel failed<br/>3 of 7 gates: macro rho -0.219,<br/>BRCA -0.405, tissue R2 57.6%"]
    C1R --> V3["V3-abs within-tissue variance filter<br/>PASSED 7/7 gates + 7/7 S1-S7<br/>removes 42.0% of excess lineage R2<br/>SHIPPED as Candidate B"]
    C1R --> V4["V4 lineage-penalized selection<br/>FAILED G3, the mechanism gate<br/>tissue R2 ROSE 0.5219 to 0.5361<br/>branch TERMINATED per docs/27 s6"]

    V3 --> LOCK
    C1B --> LOCK
    C2 --> LOCK
    C3 --> LOCK
    C4 --> LOCK
    LOCK["CNS LOCK - OPENED 2026-09-18<br/>job 323264626, CONSUMED<br/>642 samples, cannot be reused"]

    LOCK --> CNS["CNS RESULT: rank transfers weakly,<br/>calibration does NOT<br/>oracle null BEATS the model<br/>skill -1.009 GBM / -0.048 LGG"]
    CNS --> B17["B17 OOD rule did NOT flag CNS<br/>OPEN - 8 of 642, all LGG"]
    CNS --> B18["B18 HRDsum>=42 vacuous in CNS<br/>OPEN - 0 GBM, 1 LGG above cutoff"]
    CNS --> GATEA["DEPLOYMENT GATE<br/>STAYS SHUT"]

    B19["B19 shipping path audit<br/>RESOLVED adf420b - could not<br/>score at all; no clip; no ledger"] --> LOCK
    LOCK --> B20["B20 non-zero exit after<br/>successful scoring<br/>OPEN, process"]

    B7["B7 inference provenance<br/>RESOLVED 2026-09-16"] --> B10["B10 app governance<br/>OPEN, P2"]
    GATEA --> B10
    B11["B11 PBTP EPIC<br/>OPEN"] --> TRANSFER["Pediatric transfer"]
    B17 --> TRANSFER
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
    style V3 fill:#d1e7dd,stroke:#0f5132,stroke-width:3px
    style V4 fill:#f8d7da,stroke:#b02a37,stroke-width:2px
    style LOCK fill:#d1e7dd,stroke:#0f5132,stroke-width:4px
    style CNS fill:#ffe08a,stroke:#d39e00,stroke-width:3px
    style GATEA fill:#f8d7da,stroke:#b02a37,stroke-width:3px
    style B17 fill:#f8d7da,stroke:#b02a37,stroke-width:2px
    style B18 fill:#fff3cd,stroke:#d39e00
    style B19 fill:#d4edda
    style B20 fill:#fff3cd
    style B10 fill:#fff3cd
    style B11 fill:#fff3cd
```

**Critical path as of 2026-09-17 (superseded):** fix C1–C4 → re-run the array →
confirm skill improves → freeze → open the CNS lock **once** → B10 app
governance → pediatric transfer.

**Critical path as of 2026-09-18.** Everything up to and including the unlock is
done, and the unlock cannot be repeated. What remains:

1. **B17** — characterise the OOD score against source data before any
   pediatric sample is scored. This is now the highest-value open item, because
   B11's transfer is a *larger* shift than the one the flag just failed to see.
2. **B18 / C4** — cohort-composition disclosures in every presented artefact.
3. **B20** — make the one-shot job's exit status unambiguous by stage.
4. **B10** — app governance, with the B17 caveat surfaced, before any public
   demo.
5. **B15** — `renv.lock` snapshot on the analysis host.
6. **B11 / B12** — pediatric readiness. Unchanged by today, and now carrying the
   CNS result as its prior: absolute HRDsum on an unseen lineage is not
   supported.

The next-generation hypothesis — that reducing lineage imprinting *within* the
source distribution is necessary but not sufficient for cross-lineage
calibration — is **documented, not implemented**, per `docs/28` §7 and
`docs/29` §4 item 5.

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
>
> **Carried forward 2026-09-18.** Purity was populated for **632 of 642** CNS
> rows, so the pre-declared secondary residual diagnostic (`docs/28` §6) was
> computable and is reported per group in
> `results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md` §1. Nothing
> in the CNS result changes C3's status: it stays **DOWNGRADED, not closed**,
> and the probe-level purity job remains the untested piece.

### C1. Per-tissue calibration offset — **ESCALATED 2026-09-17; TESTED ON CNS AND CONFIRMED 2026-09-18**

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
> synthetic data (transfers at ρ ≈ 0.96, manufactures nothing under a pure
> tissue intercept or a noise negative control).

> ### C1c SENTINEL RAN AND FAILED 2026-09-17 ~21:00 — array `323195423`
>
> Four folds, `priority` queue, variant V2. All completed, zero errors, all
> lambdas interior. **Three of seven pre-registered gates failed, including both
> load-bearing ones. The full 30-fold rank array is NOT authorised.**
>
> | Fold | Rank ρ | Run 01 ρ | Δ |
> |---|---:|---:|---:|
> | UCEC | 0.550 | 0.703 | −0.154 |
> | BRCA | 0.261 | 0.666 | **−0.405** |
> | KICH | 0.443 | 0.543 | −0.101 |
> | THCA † | 0.125 | −0.029 | +0.154 |
>
> † negative control; its near-constant target makes rank evaluation degenerate,
> so its "gain" is the artefact the gate exists to catch.
>
> Macro ρ excluding THCA: **0.418 vs 0.637, Δ = −0.219**. Tissue identity
> explains **57.6%** of the relative score versus **52.7%** for the absolute
> prediction on the same samples, against a **0.001** floor set by the true
> relative target.
>
> **The model did not manufacture signal** — every fold beats its within-tissue
> permutation null at p ≤ 0.004. It simply learned *less* than the absolute
> model while remaining *more* tissue-confounded.
>
> **Diagnosed mechanism.** The 5,000-probe filter ranks by **pooled** variance
> and runs *before* the target is consulted, so it selects lineage-discriminating
> probes. Changing the target cannot undo a feature set already chosen for
> between-tissue variance.
>
> **C1c is therefore REFUTED for the pooled-feature variants (V1/V2).** One
> pointwise option survives and is untested: **V3**, which ranks probes by
> training-only *within-tissue* variance and attacks the diagnosed mechanism
> directly. Implementation is written and benchmarked (~2.2 min for 336,480
> probes). This is the highest-value next C1 experiment.
>
> **Pediatric consequence — see `docs/26` §7b.** Both routes to an n-of-1
> pediatric answer are now closed. Absolute HRDsum on a new tumour type carries a
> 95% predictive offset interval of **−9.7 to +8.1 units** (74% of the label's
> own IQR); among the 1,171 patients within ±10 units of the threshold, a
> tissue-level offset flips the HRD-high call for a median of **9.1%** of them
> (up to 79.2%). The relative-rank escape hatch is majority tissue identity.
> Only the **labelled** path survives: ~10 PBTP labels, which we do not have.
>
> **The CNS lock stays closed.** It opens once, after C1–C4; C1 is now more
> firmly unresolved, not less.

> ### C1c: V3-abs PASSED, AND IS THE SHIPPED MODEL — 2026-09-18
>
> **This supersedes the "V3 untested" line above, and it is a different
> experiment from the V3 defined in `docs/26` §4.** docs/26-V3 bundled *two*
> changes — V2's relative target (already refuted) *plus* the within-tissue
> filter. **V3-abs** changes **one** thing: the run-01 absolute model, entirely
> unchanged, except that the 5,000-probe filter ranks by pooled **within-tissue**
> variance on training-fold rows only, instead of pooled **total** variance.
> `docs/27` §3 records that distinction in writing before either result existed.
>
> **Sentinel: 7 of 7 pre-registered gates passed** (`docs/27` §5, thresholds
> frozen at 10:27 before the job returned; scorecard in
> `results/v3_sentinel/gate_scorecard.tsv`). **Full 30-fold LOCO: 7 of 7
> selection criteria S1–S7 passed** (`docs/28` §4, frozen at 11:58 while
> `results/v3_loco_full/` was empty; scorecard in
> `results/v3_loco_full/selection_scorecard.tsv`), selecting **Candidate B**.
>
> Full 30-fold head-to-head, macro estimator under adopted C2 clipping, from
> `results/tables/tableA2_A_vs_V3_full30.md`:
>
> | Statistic | A (run 01) | V3-abs | Δ |
> |---|---:|---:|---:|
> | Macro within-tissue Pearson | 0.5197 | 0.5036 | −0.0161 |
> | Macro within-tissue Spearman | 0.4747 | 0.4472 | −0.0275 |
> | Pooled MAE | 8.9717 | 8.4456 | −0.5261 |
> | Mean absolute tissue bias (C1) | 3.4628 | 2.9737 | −0.4891 |
> | **Tissue R² of predictions** | **0.5584** | **0.4672** | **−0.0913** |
> | Tissue R² of truth | 0.3414 | 0.3414 | 0, identical by construction |
> | Macro skill vs oracle tissue-mean null | 0.4192 | 0.6304 | +0.2112 |
> | Tissues with positive within-tissue r | 29 of 30 | 29 of 30 | 0 |
>
> **Headline.** `ΔR²_excess = 0.0913`, 95% CI [0.0824, 0.0996], bootstrap
> p < 0.0001 — **42.0% (CI 37.9–46.5%) of the excess lineage structure removed**
> beyond the 0.3414 the target itself carries. Since R²_truth is identical for
> both candidates it cancels exactly, so this equals R²_pred(A) − R²_pred(B).
> **This is the first intervention in the project to move C1 at all.**
>
> **And the honest limit, which must be stated in the same breath.** The paired
> per-tissue **ranking change is not significant**: Pearson mean −0.0161, 95% CI
> [−0.0410, 0.0051], Wilcoxon p = 0.428, **12 of 30** tissues improved; Spearman
> mean −0.0275, p = 0.477, **14 of 30** improved. So V3-abs is **a targeted
> reduction in lineage imprinting at no measurable ranking cost** — it is **not**
> a general accuracy improvement, and must never be described as one.
>
> Sentinel-vs-full divergence, per **B16**: tissue R² was 0.2991 on the four
> sentinel tissues but 0.4672 across all thirty. The direction replicated; the
> magnitude did not. This is why S1–S7 existed as a separate 30-fold gate.

> ### C1c: V4 LINEAGE-PENALISED FAILED THE MECHANISM GATE — 2026-09-18
>
> The **one** alternative source-only experiment pre-registered in `docs/27` §6 —
> supervised within-tissue HRD meta-association feature selection with a lineage
> penalty, `score_j = |z_j^meta| · c_j − λ · ℓ_j`, λ chosen in the inner folds —
> was run on the same four sentinel tissues against the same gates
> (`results/v4_sentinel/gate_scorecard.tsv`).
>
> | Gate | Threshold | Baseline | Observed | Verdict |
> |---|---|---:|---:|---|
> | G1 macro within-r ex-THCA | ≥ 0.7096 | 0.7596 | 0.7186 | PASS |
> | G2 BRCA within-r | ≥ 0.5823 | 0.6323 | 0.6125 | PASS |
> | **G3 tissue R² of predictions** | **≤ 0.4919** | 0.5219 | **0.5361** | **FAIL** |
> | G4 pooled MAE | ≤ 12.143 | 11.039 | 8.6034 | PASS |
> | G5 mean abs tissue bias | ≤ 5.674 | 5.1581 | 2.8776 | PASS |
> | G6 non-THCA not degraded | ≥ 2 of 3 | 3 of 3 | 2 of 3 | PASS |
> | G7 provenance recorded | all TRUE | — | TRUE | PASS |
>
> **Six of seven gates passed. The one that failed is the only one that was
> declared non-negotiable.** G3 is the mechanism gate: V4 exists solely to
> penalise lineage, and **tissue R² of the predictions went UP**, 0.5219 →
> 0.5361, against a ceiling of 0.4919 — *above the baseline it was supposed to
> beat*. **The experiment designed to penalise lineage produced more of it.**
>
> The λ values are the tell: selected by inner folds, λ came out **0 for BRCA**,
> **0.5 for UCEC**, **2 for KICH and THCA**. Where the inner folds had the most
> data to speak with, they chose **no penalty at all** — the penalty was not
> being selected *for*, and a supervised within-tissue association score
> evidently still concentrates on lineage-discriminating probes.
>
> **Branch terminated**, per `docs/27` §6 and the §5 decision rule: no
> re-tuning against the same four tissues, no second variant. Recorded here
> rather than deleted, because a pre-registered mechanism that fails its own
> mechanism gate is a genuinely informative negative result — it narrows the
> space of explanations for C1.
>
> Note the contrast with V3-abs, which is the scientifically interesting part:
> an **unsupervised** change to the *variance* statistic cut lineage imprinting
> by 42%, while a **supervised** penalty aimed directly at lineage did not.

> ### C1 CONFIRMED ON CNS — 2026-09-18, lock consumed, job `323264626`
>
> C1 was never a precondition for the unlock; it was **the hypothesis the unlock
> tested** (`docs/27` §2). It is now tested, on 642 genuinely untouched samples
> (135 GBM + 507 LGG), and it is **confirmed as the dominant failure mode**.
>
> The status of C1 is therefore no longer "unresolved pending investigation". It
> is: **characterised** (a lineage-dependent additive offset arising from a
> lineage-imprinted feature basis), **mitigated by ~42% within the source
> distribution** (V3-abs, above), and **still decisive out of distribution**.
>
> Primary candidate B (V3-abs), from
> `results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md`:
>
> | Group | n | Pearson [95% CI] | Spearman | MAE | bias | calib slope |
> |---|---:|---|---:|---:|---:|---:|
> | GBM | 135 | 0.320 [0.173, 0.461] | 0.279 | 8.86 | +7.99 | 0.358 |
> | LGG | 507 | 0.329 [0.234, 0.421] | 0.297 | 5.13 | +1.62 | 0.578 |
> | CNS pooled | 642 | 0.261 [0.183, 0.335] | 0.242 | 5.92 | +2.96 | 0.358 |
>
> **The oracle tissue-mean null beats the model everywhere**: skill −1.009
> (GBM), −0.048 (LGG), −0.233 (pooled). In plain terms, **"use the GBM average"
> beats the model's absolute predictions.** The null is an *oracle* — it uses the
> CNS labels and is unavailable at deployment — and it is reported precisely
> because it does not flatter.
>
> **CNS is a poor relation, not a typical held-out tissue.** Placed in the
> empirical distribution of the 30 source LOCO tissues (`docs/28` §6's primary
> interpretive frame): GBM Pearson at the **13.3rd percentile**, 28th of 32; LGG
> Pearson **16.7th**, 26th of 32; GBM mean |bias| **10th percentile**, 29th of
> 32. One counterweight, reported because omitting it would be selective: LGG
> MAE sits at the **83.3rd percentile**, 6th of 32 — but that is largely because
> LGG's own HRDsum is low and tightly distributed, which makes small absolute
> errors easy and is not evidence of calibration.
>
> **The scientific point, and the reason this is worth more than a pass/fail.**
> V3-abs's 42% reduction in excess lineage structure was measured on the
> **source** tissues, and it **did not translate into better CNS calibration**.
> Reducing lineage imprinting *within* the training distribution was **not
> sufficient** for transfer *outside* it. Whether it was *necessary* is untested
> and now untestable on this cohort — the lock is consumed.
>
> Sensitivity candidate A (run 01 pooled) was scored in the same session, as
> `docs/28` §5 requires and permits, and is reported because reporting only the
> better one is prohibited: GBM r = 0.293 [0.142, 0.446], ρ = 0.310, MAE 9.79,
> bias +9.11, slope 0.327; LGG r = 0.371 [0.278, 0.460], ρ = 0.311, MAE 5.95,
> bias +3.69, slope 0.692. **The CIs overlap heavily and neither candidate is
> clearly better on CNS.** A has the better LGG Pearson; B has the better MAE,
> bias and slope in both lineages. Under the pre-registered rule this changes
> nothing — the designation was final before unlock (`docs/29` §4 item 3) — and
> it is a useful check on over-reading the source-side selection: a 42%
> reduction in lineage imprinting bought no distinguishable rank advantage in a
> new lineage.
>
> **Verdict: outcome A of the three pre-declared in `docs/28` §7** — *rank
> transfers weakly, absolute cross-lineage calibration does not*. **The
> deployment gate (`docs/27` §2A) stays shut.**
>
> Two defects surfaced by this evaluation are logged as **B17** (the OOD
> reportability rule did not flag CNS) and **B18** (the 42 threshold is vacuous
> in CNS).
>
> Consequences for C1b, unchanged and now more pointed: the labelled few-shot
> path remains the only supported route to absolute HRDsum in a new lineage, and
> we still have no PBTP labels.



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

> **Carried forward 2026-09-18 — with a correction that matters.** Clipping was
> recorded as ADOPTED on 2026-09-17, but **`predict_frozen.R` never applied it**
> until commit `adf420b` (see **B19** item 1). The adopted rule was true of the
> evaluation path and false of the shipping path for a day. The CNS evaluation
> ran under the corrected path, stamped
> `C2_ADOPTED:pmax(inverse_clip(raw),0)` in
> `results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md`, with the
> pre-clip value retained in `predicted_reference_HRDsum_raw`. C2 stays
> **PARTLY RESOLVED**: the clip is in place everywhere it should be, the hurdle
> model is still future work.

### C4. Ovarian cohort is n = 10 — **OPEN (disclosure, not a code fix)**

OV is the canonical HRD cancer and the main clinical application. TCGA ovarian
methylation is mostly 27k-array, excluded by the 450k bridge. HRD-high behaviour
is therefore inferred from UCEC/BRCA/STAD.

**Resolution.** State this limitation in every presentation. Optionally rebuild
the bridge to include 27k probes, accepting a much smaller probe intersection.

> **Carried forward 2026-09-18, and no longer alone.** C4 is a
> cohort-composition limit: the tissue where the claim matters most is nearly
> absent from development. **B18** is the same class of problem on the
> evaluation side — CNS has **0 GBM and 1 LGG** sample at the exploratory
> HRDsum ≥ 42 cutoff, so the high-HRD regime is nearly absent from the external
> test too. Taken together: the project has developed and externally tested a
> model largely **in the low-HRD regime**, and the disclosure should say that
> rather than listing two separate footnotes.

---

**Superseded 2026-09-17.** The old critical path (dry run → array → merge → read
the pooled tissue-mean null) is **complete**. The pooled null was read: skill
+0.085, within-tissue r = 0.612, permutation p = 0.001. The answer to "is there
a result worth reporting at all" is **yes, reframed** — see
`docs/24_RESULTS_LOCO_RUN01.md` §6. The path forward is C1–C4 above.

**Superseded again 2026-09-18.** The C1–C4 path is also complete, in the sense
that matters: the lock was opened once and C1 was tested rather than assumed.
The project's answer to "does this transfer to an unseen lineage" is **rank
weakly yes, calibration no** — a complete, honest, pre-declared answer
(`docs/28` §7, outcome A), not a project failure. What remains open is **B17**
(the OOD rule that should have warned us and did not), **B18** and **C4**
(cohort composition), **B20** (job-status hygiene), **B10** (app governance),
**B15** (`renv.lock`), and **B11/B12** (pediatric readiness).

**Nothing on this page has been deleted.** The refuted V1/V2 sentinel, the
failed label-free C1 attack, the rejected log1p transform and the terminated V4
branch are all retained above. Each one narrowed the hypothesis space, and the
V3-abs result is only interpretable against them.
