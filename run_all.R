# ============================================================
# run_all.R
# Entry point do pacote de reprodução
# ============================================================

options(stringsAsFactors = FALSE)
source(file.path("R", "00_config.R"), chdir = FALSE)

scripts <- c(
  "01_build_static_exposure.R",
  "02_static_results.R",
  "03_dynamic_results.R",
  "04_boundary_network_results.R",
  "05_appendix_robustness.R",
  "06_export_paper_outputs.R"
)

RUN_STARTED_AT <- Sys.time()
start_time <- RUN_STARTED_AT
for (s in scripts) {
  message("\n============================================================")
  message("Executando: R/", s)
  message("============================================================")
  source(file.path(REPRO_ROOT, "R", s), chdir = FALSE)
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
writeLines(
  c(
    paste0("Reprodução concluída em ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Tempo total (min): ", round(as.numeric(elapsed), 2)),
    paste0("R version: ", R.version.string)
  ),
  file.path(LOG_DIR, "run_all.log")
)
message("\nReprodução concluída. Ver outputs/tables e outputs/figures.")
