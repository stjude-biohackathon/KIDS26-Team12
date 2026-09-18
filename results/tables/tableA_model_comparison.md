# Table A. Model comparison (STUB -- awaiting candidate selection)

| Model | Description | Pooled_Pearson | Pooled_Spearman | Macro_within_tissue_Pearson | Macro_within_tissue_Spearman | Macro_MAE | Skill_vs_tissue_mean_null | Tissue_R2_of_prediction | Mean_abs_tissue_offset | Selected |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
*(stub -- no rows yet; columns are defined above)*


**STATUS: STUB.** Column headers are frozen; rows are added by the main agent
after candidate selection completes. No model has been selected at the time this
file was generated, and no row may be added from a run that saw the held-out
cancer type or the locked CNS cohort.

Metric definitions match `results/loco_run01/`: `Macro_*` statistics are computed
inside each held-out cancer type and then averaged unweighted over cancer types;
`Skill_vs_tissue_mean_null` is 1 - MAE_model / MAE_tissue_mean_null;
`Tissue_R2_of_prediction` is the one-way ANOVA R2 of cancer type on the clipped
prediction. All predictions must have the adopted C2 rule `pmax(pred, 0)` applied.
