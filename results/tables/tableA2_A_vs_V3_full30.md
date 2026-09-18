# Table A2 - Candidate A vs V3-abs, full 30-fold LOCO

Development LOCO, 7,065 samples, 30 non-CNS cancer types. Generated 2026-09-18 14:23 CDT, seed 5813, 10,000 bootstrap replicates.
Per-tissue statistics under adopted C2 clipping, macro-averaged; same estimator for both candidates.

| metric | model_A | v3_abs | difference | ci_95 |
| --- | --- | --- | --- | --- |
| Macro Pearson (within-tissue) | 0.5197 | 0.5036 | -0.0161 | - |
| Macro Spearman (within-tissue) | 0.4747 | 0.4472 | -0.0275 | - |
| Pooled MAE | 8.9717 | 8.4456 | -0.5261 | - |
| Macro skill vs oracle tissue-mean null | 0.4192 | 0.6304 | 0.2112 | - |
| Mean absolute tissue bias | 3.4628 | 2.9737 | -0.4891 | - |
| Tissue R2 of predictions | 0.5584 | 0.4672 | -0.0913 | A [0.5461, 0.5736] | B [0.4539, 0.4837] |
| Tissue R2 of truth | 0.3414 | 0.3414 | 0 (identical by construction) | - |
| Excess tissue R2 (pred - truth) | 0.2171 | 0.1258 | -0.0913 | A [0.2018, 0.2310] | B [0.1104, 0.1407] |
| Fraction of excess removed | - | 42.0% | 0.0913 removed | [37.9, 46.5] |
| Tissues with positive r | 29 of 30 | 29 of 30 | +0 | - |
| Median per-tissue delta Pearson | - | -0.0022 | 12 of 30 improved | [-0.0410, 0.0051] |
| Median per-tissue delta Spearman | - | -0.0095 | 14 of 30 improved | [-0.0777, 0.0078] |
| Paired difference, Pearson | - | mean -0.0161 | p = 0.428 (Wilcoxon) | [-0.0410, 0.0051] |
| Paired difference, Spearman | - | mean -0.0275 | p = 0.477 (Wilcoxon) | [-0.0777, 0.0078] |

## Headline quantity

delta_R2_excess = 0.0913, 95% CI [0.0824, 0.0996], bootstrap p = 0.0000

R2_truth (0.3414) is identical for both candidates and cancels exactly, so this equals R2_pred_A - R2_pred_B.
V3-abs removes 42.0% of the excess lineage structure that Candidate A carries beyond the target's own.

## Paired ranking change

Pearson:  mean -0.0161, 95% CI [-0.0410, 0.0051], Wilcoxon p = 0.428, 12 of 30 tissues improved
Spearman: mean -0.0275, 95% CI [-0.0777, 0.0078], Wilcoxon p = 0.477, 14 of 30 tissues improved
