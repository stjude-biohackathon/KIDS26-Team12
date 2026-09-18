# Run once on the analysis host; review package versions and test before snapshot.
if(!requireNamespace("renv",quietly=TRUE)) install.packages("renv",repos="https://cloud.r-project.org")
renv::init(bare=TRUE)
renv::install(c("glmnet","data.table","jsonlite","shiny","ggplot2","bioc::dittoViz"))
# Add Bioconductor methylation callers only for the intensity branch, with reviewed versions.
# After Rscript tests/smoke_model.R passes: renv::snapshot(prompt=FALSE)
writeLines(capture.output(sessionInfo()),"renv_setup_session.txt")
