# Run once on the analysis host; review package versions and test before snapshot.
if(!requireNamespace("renv",quietly=TRUE)) install.packages("renv",repos="https://cloud.r-project.org")
renv::init(bare=TRUE)

# B15: this list must cover every package the code actually loads, including
# ones reached via pkg::fn without library(). Two were missing and are now here:
#
#   matrixStats  R/model.R::fit_preprocess uses colMedians/colVars/colMeans2.
#                This is load-bearing, not cosmetic - it is what took nested
#                LOCO preprocessing from ~42 h to tractable (B3). On a restore
#                without it the code errors rather than silently slowing down,
#                but either way the environment is not reproducible.
#   digest       app/app.R uses digest::digest(algo="sha256") to verify a
#                predictions table against its provenance sidecar (B10). If it
#                is absent the app warns that the checksum was NOT verified and
#                continues, so its absence degrades a security check quietly -
#                exactly the failure mode worth pinning against.
#
# jsonlite is used in app/app.R via jsonlite:: and is guarded there by an
# explicit requireNamespace() check rather than being assumed present.
# The interactive demo also loads shinydashboard, plotly, DT and VizModules.
renv::install(c("glmnet","data.table","jsonlite","shiny","matrixStats","digest","shinydashboard","plotly","DT","VizModules"))

# Add Bioconductor methylation callers only for the intensity branch, with reviewed versions.
# After the R test suites pass:
#   Rscript tests/smoke_model.R
#   Rscript tests/test_provenance_gate.R
#   Rscript tests/test_calibration.R
#   Rscript tests/test_c1_rank_model.R
#   Rscript tests/test_app_hrd_scores.R
# then: renv::snapshot(prompt=FALSE)
writeLines(capture.output(sessionInfo()),"renv_setup_session.txt")
