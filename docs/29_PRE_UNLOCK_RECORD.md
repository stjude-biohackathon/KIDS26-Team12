# 29 — Pre-unlock audit record

**Written and committed 2026-09-18, before any CNS HRDsum label was read.**

`results/LOCKED_EVALUATION_LEDGER.tsv` does not exist at the time of this
commit. Its absence is the machine-checkable proof that the locked cohort has
never been scored. If you are auditing this later, confirm this file's commit
precedes the first ledger row.

---

## 1. Designation

Per `docs/28` §4, the S1–S7 selection rule — frozen at commit `28bb913`, 11:58
CDT, while `results/v3_loco_full/` was empty — was scored at commit `5ee4bac`
and **passed 7 of 7**, selecting Candidate B.

| Role | Candidate | Artifact | sha256 | feature_rank |
|---|---|---|---|---|
| **PRIMARY** | **B — V3-abs** | `results/frozen_v3_2026-09-18/frozen_nonCNS.rds` | `df9f7e82b0b371b81ecca6a1d99a1a50128a95e78daa63307d5f3b0633a4aa72` | `within_tissue` |
| Sensitivity | A — run 01 | `results/frozen_2026-09-18/frozen_nonCNS.rds` | `0e975cd6bc7caededffecf34332a3afd250e01bfea54fccab82266021f325763` | `pooled` |

Designation confirmed by the project lead before unlock, on the record, with the
counter-argument for A (independent replication in `pipeline_glmnet/`, and a
ranking comparison that was a statistical tie) explicitly considered and
declined in favour of following the pre-registered rule.

Both artifacts were frozen **before** this designation: A at 17:06:38 UTC, B at
17:53:28 UTC, both from git commit `28bb913`. Neither has been modified since,
and neither can be modified after unlock.

## 2. Rationale for the primary choice

`ΔR²_excess = 0.0913`, 95% CI [0.0824, 0.0996], bootstrap p < 0.0001 — V3-abs
carries 42.0% (CI 37.9–46.5%) less excess lineage structure than A beyond the
0.3414 the target itself carries. Macro skill against the oracle tissue-mean
null rises 0.419 → 0.630. The paired per-tissue ranking change is
indistinguishable from zero (Pearson p = 0.43, Spearman p = 0.48).

Since C1 — lineage-dependent calibration offset — is the failure mode the CNS
test is designed to probe, the candidate that demonstrably carries less lineage
structure is the scientifically appropriate one to send into an unseen lineage.

## 3. Metrics, fixed before unlock

Exactly as pre-declared in `docs/28` §6: GBM and LGG reported separately first,
pooled CNS second. Ranking (Pearson, Spearman, bootstrap CIs, seed 3391);
absolute (MAE, RMSE, median AE); calibration (bias, intercept, slope, range
compression); the within-CNS tissue-mean null **labelled as an oracle**; OOD
reportability; purity if available. Placement of GBM and LGG within the
30-tissue source LOCO distribution is the primary interpretive frame.

HRDsum ≥ 42 is exploratory only. The score is not a probability.

## 4. Binding commitments

1. **No refitting, recalibration, preprocessing change, CpG change,
   alpha/lambda change, or post-hoc clipping change after unlock.**
2. **Both candidates' results will be reported**, whatever they show. Reporting
   only the better one is prohibited by `docs/28` §5 clause 3.
3. **The CNS result will not be used to select a model.** The designation above
   is final and was made before any label was visible.
4. **No third model** may be built after unlock.
5. If transfer is poor, that is the finding. The next-generation hypothesis gets
   documented, not implemented today.

## 5. Confirmation

No CNS outcome label was inspected during candidate selection. Every source-side
decision — the V3-abs gates (`docs/27` §5), the S1–S7 selection rule
(`docs/28` §4), and the V4 termination — was made on development data only, and
each was committed to git before the results it would judge existed.

Technical properties of the CNS feature matrix (642 rows: 135 GBM, 507 LGG;
probe dimensions; orientation; provenance) were verified, which `docs/28` and
the handoff both permit. No outcome column was read.

**Git commit at unlock:** recorded in the ledger row written by
`scripts/predict_frozen.R` at scoring time.
