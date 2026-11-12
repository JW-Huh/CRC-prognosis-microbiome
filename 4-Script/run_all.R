################################################################################
##  run_all.R  -  reproduce every table and figure from the raw inputs
################################################################################
##  - run from the repository root:  Rscript 4-Script/run_all.R
##  - total run time ~ 15 min (Monte-Carlo CV in 02 and the bootstrap in 03 dominate)
##  - each script is self-contained and can also be run on its own in order
################################################################################
scripts <- c("01_preprocessing.R",
             "02_screening_prognostic_bacteria.R",
             "03_microbial_hazard_score.R",
             "04_community_overview_figures.R",
             "05_picrust2_pathways.R",
             "06_diet_exploratory_2022.R")
for (s in scripts) {
  cat("\n\n==========", s, "==========\n")
  t0 <- Sys.time()
  source(file.path("4-Script", s), echo = FALSE)
  cat("---- finished in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}
