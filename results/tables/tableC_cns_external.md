# Table C. Locked CNS external test (STUB -- cohort sealed, not yet scored)

| Cancer | N | Mean_observed_HRD | Mean_predicted_HRD | Pearson | Spearman | MAE | RMSE | Bias_mean_signed_error | Skill_vs_tissue_mean_null | OOD_abstention | Prereg_gate_outcome |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
*(stub -- no rows yet; columns are defined above)*


**STATUS: STUB. THE CNS COHORT IS LOCKED.** 642 GBM + LGG specimens are held sealed and
will be scored exactly once, after the candidate model is frozen. No GBM or LGG HRDsum
value has been read, predicted, plotted or summarised by any script in the figure/table
layer -- `R/figures.R::assert_development()` aborts if either type appears.

Column headers are frozen now so that the analysis is pre-committed. The main agent fills
the rows after the single external run. Results must be reported against the
pre-registered gates in `docs/27_V3_PREREGISTRATION.md`, whatever they show; no reselection,
refitting or gate revision is permitted after the cohort is opened.
