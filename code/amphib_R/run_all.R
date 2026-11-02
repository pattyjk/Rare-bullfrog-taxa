## Run the full pipeline from the project root: Rscript run_all.R
## Steps 01 and 04 are the slow ones (~1 min each); their outputs are already in data/derived/, so you can
## start from any later script. Set RERUN_PREP <- FALSE to skip 01 and 04 and reuse data/derived/.
RERUN_PREP <- TRUE
steps <- c("R/01_prepare_meta.R", "R/02_meta_genus_abundance.R", "R/03_bullfrog_genus_abundance.R", "R/04_amphibac_match.R",
           "R/05_sOTU_function_variance.R", "R/06_cultured_share.R", "R/07_abundance_class_function.R", "R/08_source_host.R")
if (!RERUN_PREP) steps <- setdiff(steps, c("R/01_prepare_meta.R", "R/04_amphibac_match.R"))
for (s in steps) { message("\n===== ", s, " ====="); t0 <- Sys.time(); source(s, local = new.env()); message("done in ", round(difftime(Sys.time(), t0, units = "secs")), " s") }
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
