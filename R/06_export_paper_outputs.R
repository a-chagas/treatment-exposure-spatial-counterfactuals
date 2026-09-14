# ============================================================
# 06_export_paper_outputs.R
# Consolida e verifica os arquivos efetivamente usados no manuscrito
# ============================================================

message("[06] Consolidando outputs do paper...")

# Effective-support table added in Revision 2.
if (!file.exists(DYNAMIC_RESULTS_RDS)) stop("Dynamic results object not found.")
dyn <- readRDS(DYNAMIC_RESULTS_RDS)
eff <- data.table::as.data.table(dyn$effective_support)
pooled_eff <- data.table::as.data.table(dyn$pooled_effective_support)

eff[, stack_order := match(as.character(stack_label), c(
  "Central ULEZ, Apr. 2019", "Inner London expansion, Oct. 2021", "London-wide expansion, Aug. 2023"
))]
data.table::setorder(eff, stack_order)
eff[, paper_stack_label := data.table::fcase(
  as.character(stack_label) == "Central ULEZ, Apr. 2019", "2019 Central ULEZ",
  as.character(stack_label) == "Inner London expansion, Oct. 2021", "2021 Inner London",
  as.character(stack_label) == "London-wide expansion, Aug. 2023", "2023 London-wide",
  default = NA_character_
)]

if (nrow(eff) != 3L || anyNA(eff$paper_stack_label) || nrow(pooled_eff) != 1L) {
  stop("Unexpected effective-support structure for the 3 km paper table.")
}

eff_lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Effective identifying support for the currently untreated exposure component, 0--3 km benchmark}",
  "\\label{tab:ulez_effective_support_3km}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabular}{lrrrr}", "\\toprule",
  " & Nominal exposed & Raw target & Nonzero residualized & \\\\",
  "Stack & receivers & stations & target stations & $N_{\\mathrm{eff}}$ \\\\",
  "\\midrule",
  paste0(
    eff$paper_stack_label, " & ", eff$nominal_exposed_receivers, " & ",
    eff$raw_target_stations, " & ", eff$nonzero_residualized_target_stations, " & ",
    formatC(eff$n_eff, 2, format = "f"), " \\\\"
  ),
  "\\midrule",
  paste0(
    "Pooled across stacks & ",
    sum(eff$nominal_exposed_receivers), " & ",
    pooled_eff$raw_nonzero_regressor_stations, " & ",
    pooled_eff$nonzero_residualized_target_stations, " & ",
    formatC(pooled_eff$n_eff, 2, format = "f"), " \\\\"
  ),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
  paste0(
    "\\item \\emph{Notes:} ``Nominal exposed receivers'' counts currently untreated stations with positive entering-cohort exposure. ",
    "``Raw target stations'' counts stations for which the post-treatment target regressor is nonzero. ",
    "``Nonzero residualized target stations'' counts stations with nonzero target residuals after partialling out the fixed effects and the other decomposition components; it is \\emph{not} a count of raw treatment-exposure support. ",
    "Residualization can assign nonzero target residuals to stations whose raw target regressor is zero. ",
    "Support should therefore be assessed using the raw counts together with the distribution of residualized sum of squares, $N_{\\mathrm{eff}}$, the largest $q_i$ shares, and Frisch--Waugh--Lovell coefficient contributions. ",
    "In the 2023 stack, the largest station share is ", formatC(eff$max_share[eff$stack_order == 3L], 3, format = "f"),
    "; for the pooled model it is ", formatC(pooled_eff$max_share, 3, format = "f"), "."
  ),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
)
writeLines(eff_lines, file.path(TABLE_DIR, "table_ulez_effective_support_3km.tex"))

# Expected paper outputs. If one is missing, the public patch should fail visibly.
expected_tables <- c(
  "table_ulez_protocol_implementation.tex",
  "table_ulez_support_protocol.tex",
  "table_ulez_design_diagnostics.tex",
  "table_ulez_no_restricted_unrestricted.tex",
  "table_ulez_across_matrices_protocol.tex",
  "table_ulez_dynamic_support_3km.tex",
  "table_ulez_descriptive_statistics.tex",
  "table_ulez_benchmark_implied_effects.tex",
  "table_ulez_spatial_inference.tex",
  "table_ulez_exposure_bins.tex",
  "table_ulez_conventional_robustness.tex",
  "table_ulez_effective_support_3km.tex",
  "table_ulez_boundary_mapping_correspondence.tex",
  "table_ulez_boundary_dynamic_support.tex",
  "table_ulez_boundary_effects.tex",
  "table_ulez_predetermined_network.tex",
  "table_ulez_announcement_2019.tex",
  "table_ulez_exposure_timing_2021.tex"
)

expected_figures <- c(
  "fig_ulez_benchmark_pretrend_cells_with_D1S0.png",
  "fig_ulez_benchmark_restricted_unrestricted.png",
  "fig_ulez_across_matrices_protocol.png",
  "fig_12_01a_dynamic_static_decomp_implied_3km.png",
  "fig_12_02a_dynamic_static_decomp_event_study_3km.png",
  "Figure_S1.pdf"
)

missing_tables <- expected_tables[!file.exists(file.path(TABLE_DIR, expected_tables))]
missing_figures <- expected_figures[!file.exists(file.path(FIGURE_DIR, expected_figures))]

actual_tables <- list.files(TABLE_DIR, all.files = FALSE, no.. = TRUE)
actual_figures <- list.files(FIGURE_DIR, all.files = FALSE, no.. = TRUE)
unexpected_tables <- setdiff(actual_tables, expected_tables)
unexpected_figures <- setdiff(actual_figures, expected_figures)

if (length(missing_tables) || length(missing_figures) ||
    length(unexpected_tables) || length(unexpected_figures)) {
  stop(
    "Reprodução incompleta ou diretório de outputs não canônico.\n",
    if (length(missing_tables)) paste0("Tabelas ausentes:\n", paste(missing_tables, collapse = "\n"), "\n") else "",
    if (length(missing_figures)) paste0("Figuras ausentes:\n", paste(missing_figures, collapse = "\n"), "\n") else "",
    if (length(unexpected_tables)) paste0("Tabelas inesperadas:\n", paste(unexpected_tables, collapse = "\n"), "\n") else "",
    if (length(unexpected_figures)) paste0("Figuras inesperadas:\n", paste(unexpected_figures, collapse = "\n")) else ""
  )
}

if (exists("RUN_STARTED_AT", inherits = TRUE)) {
  expected_paths <- c(
    file.path(TABLE_DIR, expected_tables),
    file.path(FIGURE_DIR, expected_figures)
  )
  mtimes <- file.info(expected_paths)$mtime
  stale <- is.na(mtimes) | mtimes < (RUN_STARTED_AT - 2)
  if (any(stale)) {
    stop(
      "Outputs antigos foram encontrados durante a execução atual:\n",
      paste(expected_paths[stale], collapse = "\n")
    )
  }
}

manifest <- data.table::rbindlist(list(
  data.table::data.table(type = "table", file = expected_tables),
  data.table::data.table(type = "figure", file = expected_figures)
))
manifest[, path := ifelse(type == "table", file.path("outputs", "tables", file), file.path("outputs", "figures", file))]
data.table::fwrite(manifest, file.path(REPRO_ROOT, "outputs", "manifest.csv"))

message("[06] ", nrow(manifest), " paper outputs verified.")
