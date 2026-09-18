# Table A - Source-domain model comparison

Development LOCO, 7,065 samples, 30 non-CNS cancer types. Scored 2026-09-18 14:00 CDT.
Per-tissue statistics computed under the adopted C2 clipping, then macro-averaged.
This macro estimator is not interchangeable with the pooled-within estimator (0.612) in docs/24.

| model | feature_strategy | n_features | macro_pearson | macro_spearman | pooled_MAE | mean_abs_tissue_bias | tissue_R2_predicted | tissue_R2_observed | n_tissues_positive_r | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A: run 01 baseline | Pooled total variance (unsupervised) | 5000 | 0.5197 | 0.4747 | 8.9717 | 3.4628 | 0.5584 | 0.3414 | 29 | Incumbent; independently replicated in pipeline_glmnet/ |
| B: V3-abs | Pooled within-tissue variance (unsupervised) | 5000 | 0.5036 | 0.4472 | 8.4456 | 2.9737 | 0.4672 | 0.3414 | 29 | Selection rule S1-S7: 7 of 7 passed |

**Primary candidate: B (V3-abs, within_tissue)** (S1-S7: 7 of 7 passed).
