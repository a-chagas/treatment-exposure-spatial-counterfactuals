# ============================================================
# 05_appendix_robustness.R
# Descritivas, implied effects, inferência espacial e forma funcional
# ============================================================

message("[05] Gerando tabelas de robustez do apêndice...")
if (!file.exists(STATIC_PANEL_FILE)) stop("Execute 01_build_static_exposure.R primeiro.")

PANEL_FILE <- STATIC_PANEL_FILE
TAB <- file.path(DERIVED_DIR, "appendix_core_tables")
RDS <- file.path(DERIVED_DIR, "appendix_models")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Read and prepare data --------------------------------

panel <- readr::read_csv(PANEL_FILE, show_col_types = FALSE) |>
  janitor::clean_names()

required_vars <- c("code", "month", "no2", "d_active", "s_bench_0_3km")
missing_vars <- setdiff(required_vars, names(panel))
if (length(missing_vars) > 0) {
  stop("Missing required variable(s): ", paste(missing_vars, collapse = ", "))
}

panel <- panel |>
  mutate(
    code = stringr::str_to_upper(as.character(code)),
    month = as.Date(month),
    D_active = as.integer(d_active),
    log_y = ifelse("log_y" %in% names(panel) & !is.na(log_y), log_y, log(no2)),
    S_bench = s_bench_0_3km,
    DS_bench = D_active * S_bench,
    US_bench = (1 - D_active) * S_bench
  ) |>
  filter(
    is.finite(no2),
    is.finite(log_y),
    is.finite(D_active),
    is.finite(S_bench),
    !is.na(code),
    !is.na(month)
  )

panel_central <- panel |>
  filter(month <= as.Date("2021-09-01"))

if (nrow(panel_central) == 0) {
  stop("Central-ULEZ-only subsample is empty. Check the month variable.")
}

# ---- 3. Formatting helpers ----------------------------------

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", formatC(x, digits = 0, format = "f", big.mark = ","))
}

stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

effect_pct <- function(x) 100 * (exp(x) - 1)

safe_nobs <- function(model_obj) {
  out <- tryCatch(
    as.numeric(stats::nobs(model_obj)),
    error = function(e) {
      tryCatch(length(stats::residuals(model_obj)), error = function(e2) NA_real_)
    }
  )
  if (length(out) == 0) NA_real_ else out[1]
}

safe_within_r2 <- function(model_obj) {
  out <- tryCatch(
    as.numeric(fixest::fitstat(model_obj, "wr2")),
    error = function(e) NA_real_
  )
  if (length(out) == 0) NA_real_ else out[1]
}

coef_table_safe <- function(model_obj) {
  tryCatch(
    fixest::coeftable(model_obj),
    error = function(e) {
      sm <- summary(model_obj)
      if (!is.null(sm$coeftable)) sm$coeftable else stop(e)
    }
  )
}

mean_pos_exposure <- function(data, component = c("treated", "untreated", "all")) {
  component <- match.arg(component)
  x <- switch(
    component,
    treated = data$S_bench[data$D_active == 1 & data$S_bench > 0],
    untreated = data$S_bench[data$D_active == 0 & data$S_bench > 0],
    all = data$S_bench[data$S_bench > 0]
  )
  if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
}

summarize_vector <- function(x, label) {
  tibble::tibble(
    variable = label,
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    p25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
    p50 = as.numeric(stats::quantile(x, 0.50, na.rm = TRUE, names = FALSE)),
    p75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
    max = max(x, na.rm = TRUE)
  )
}

effect_cell <- function(model_obj, term, model_name, scale = c("log", "level"),
                        exposure_scale = 1, digits = 2) {
  scale <- match.arg(scale)
  ct <- coef_table_safe(model_obj)

  if (!term %in% rownames(ct)) {
    return(tibble::tibble(
      model = model_name,
      term = term,
      beta = NA_real_,
      se_beta = NA_real_,
      p_value = NA_real_,
      exposure_scale = exposure_scale,
      implied_effect_linear = NA_real_,
      implied_se_linear = NA_real_,
      reported_effect = NA_real_,
      reported_se = NA_real_,
      cell = ""
    ))
  }

  beta <- as.numeric(ct[term, 1])
  se_beta <- as.numeric(ct[term, 2])
  p_value <- as.numeric(ct[term, 4])

  implied_effect_linear <- beta * exposure_scale
  implied_se_linear <- se_beta * exposure_scale

  if (scale == "log") {
    reported_effect <- effect_pct(implied_effect_linear)
    reported_se <- 100 * exp(implied_effect_linear) * implied_se_linear
  } else {
    reported_effect <- implied_effect_linear
    reported_se <- implied_se_linear
  }

  tibble::tibble(
    model = model_name,
    term = term,
    beta = beta,
    se_beta = se_beta,
    p_value = p_value,
    exposure_scale = exposure_scale,
    implied_effect_linear = implied_effect_linear,
    implied_se_linear = implied_se_linear,
    reported_effect = reported_effect,
    reported_se = reported_se,
    cell = paste0(fmt_num(reported_effect, digits), stars(p_value), " (", fmt_num(reported_se, digits), ")")
  )
}

# ---- 4. Appendix descriptive statistics ----------------------

desc_stats <- dplyr::bind_rows(
  summarize_vector(panel$no2, "$NO_2$"),
  summarize_vector(panel$log_y, "$\\log(NO_2)$"),
  summarize_vector(panel$D_active, "$D_{it}$"),
  summarize_vector(panel$S_bench, "$S_{it}$"),
  summarize_vector(panel$DS_bench, "$D_{it}S_{it}$"),
  summarize_vector(panel$US_bench, "$(1-D_{it})S_{it}$"),
  summarize_vector(panel$S_bench[panel$S_bench > 0], "$S_{it}\\mid S_{it}>0$"),
  summarize_vector(panel$S_bench[panel$D_active == 1 & panel$S_bench > 0],
                   "$S_{it}\\mid D_{it}=1,S_{it}>0$"),
  summarize_vector(panel$S_bench[panel$D_active == 0 & panel$S_bench > 0],
                   "$S_{it}\\mid D_{it}=0,S_{it}>0$")
)

readr::write_csv(desc_stats, file.path(TAB, "table_ulez_descriptive_statistics.csv"))

desc_lines <- c(
"\\begin{table}[!htbp]\\centering",
"\\caption{Descriptive statistics for the ULEZ empirical illustration}",
"\\label{tab:ulez_descriptive_statistics}",
"\\begin{threeparttable}",
"\\small",
"\\setlength{\\tabcolsep}{3pt}",
"\\renewcommand{\\arraystretch}{1.08}",
"\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}Xrrrrrrrr@{}}",
"\\toprule",
"Variable & $N$ & Mean & SD & Min & $p25$ & $p50$ & $p75$ & Max \\\\",
"\\midrule",
paste0(
  desc_stats$variable, " & ",
  fmt_int(desc_stats$n), " & ",
  fmt_num(desc_stats$mean, 3), " & ",
  fmt_num(desc_stats$sd, 3), " & ",
  fmt_num(desc_stats$min, 3), " & ",
  fmt_num(desc_stats$p25, 3), " & ",
  fmt_num(desc_stats$p50, 3), " & ",
  fmt_num(desc_stats$p75, 3), " & ",
  fmt_num(desc_stats$max, 3), " \\\\"
),
"\\bottomrule",
"\\end{tabularx}",
"\\begin{tablenotes}[flushleft]",
"\\footnotesize",
"\\item Notes: $NO_2$ is measured in $\\mu$g/m$^3$. $S_{it}$ is the benchmark row-standardized 0--3 km exposure index. The last three rows report exposure intensity among exposed observations, including the treated exposed and currently untreated exposed subsamples used to compute $\\bar S^t$ and $\\bar S^u$. The empirical specifications include station and month fixed effects and no additional time-varying covariates. Source: London Air Quality Network/Imperial College London monitoring data and Greater London Authority/Transport for London ULEZ boundary files. Author's calculations.",
"\\end{tablenotes}",
"\\end{threeparttable}",
"\\end{table}"
)

writeLines(desc_lines, file.path(TAB, "table_ulez_descriptive_statistics.tex"))

# ---- 5. Benchmark model comparison: slopes and implied effects

m_noexp <- fixest::feols(
  log_y ~ D_active | code + month,
  data = panel,
  cluster = ~code,
  warn = FALSE,
  notes = FALSE
)

m_restricted <- fixest::feols(
  log_y ~ D_active + S_bench | code + month,
  data = panel,
  cluster = ~code,
  warn = FALSE,
  notes = FALSE
)

m_unrestricted <- fixest::feols(
  log_y ~ D_active + DS_bench + US_bench | code + month,
  data = panel,
  cluster = ~code,
  warn = FALSE,
  notes = FALSE
)

sbar_all <- mean_pos_exposure(panel, "all")
sbar_treated <- mean_pos_exposure(panel, "treated")
sbar_untreated <- mean_pos_exposure(panel, "untreated")

benchmark_effects_long <- dplyr::bind_rows(
  effect_cell(m_noexp, "D_active", "No-exposure DiD", "log", 1) |>
    mutate(component = "$\\beta^d$", exposure_label = "", exposure_mean = 1),

  effect_cell(m_restricted, "D_active", "Restricted exposure", "log", 1) |>
    mutate(component = "$\\beta^d$", exposure_label = "", exposure_mean = 1),

  effect_cell(m_restricted, "S_bench", "Restricted exposure", "log", sbar_all) |>
    mutate(component = "$\\beta^r\\bar S$", exposure_label = "$\\bar S\\mid S>0$", exposure_mean = sbar_all),

  effect_cell(m_unrestricted, "D_active", "Unrestricted exposure", "log", 1) |>
    mutate(component = "$\\beta^d$", exposure_label = "", exposure_mean = 1),

  effect_cell(m_unrestricted, "DS_bench", "Unrestricted exposure", "log", sbar_treated) |>
    mutate(component = "$\\beta^t\\bar S^t$", exposure_label = "$\\bar S^t$", exposure_mean = sbar_treated),

  effect_cell(m_unrestricted, "US_bench", "Unrestricted exposure", "log", sbar_untreated) |>
    mutate(component = "$\\beta^u\\bar S^u$", exposure_label = "$\\bar S^u$", exposure_mean = sbar_untreated)
)

readr::write_csv(benchmark_effects_long, file.path(TAB, "table_ulez_benchmark_implied_effects.csv"))

benchmark_wide <- benchmark_effects_long |>
  select(model, component, cell) |>
  tidyr::pivot_wider(names_from = component, values_from = cell) |>
  arrange(match(model, c("No-exposure DiD", "Restricted exposure", "Unrestricted exposure")))

for (nm in c("$\\beta^d$", "$\\beta^r\\bar S$", "$\\beta^t\\bar S^t$", "$\\beta^u\\bar S^u$")) {
  if (!nm %in% names(benchmark_wide)) benchmark_wide[[nm]] <- ""
  benchmark_wide[[nm]][is.na(benchmark_wide[[nm]])] <- ""
}

benchmark_lines <- c(
"\\begin{table}[!htbp]\\centering",
"\\caption{Benchmark ULEZ estimates and implied exposure effects}",
"\\label{tab:ulez_benchmark_implied_effects}",
"\\begin{threeparttable}",
"\\small",
"\\setlength{\\tabcolsep}{3pt}",
"\\renewcommand{\\arraystretch}{1.08}",
"\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}Xrrrr@{}}",
"\\toprule",
"Specification & $\\beta^d$ & $\\beta^r\\bar S$ & $\\beta^t\\bar S^t$ & $\\beta^u\\bar S^u$ \\\\",
"\\midrule",
paste0(
  benchmark_wide$model, " & ",
  benchmark_wide$`$\\beta^d$`, " & ",
  benchmark_wide$`$\\beta^r\\bar S$`, " & ",
  benchmark_wide$`$\\beta^t\\bar S^t$`, " & ",
  benchmark_wide$`$\\beta^u\\bar S^u$`, " \\\\"
),
"\\bottomrule",
"\\end{tabularx}",
"\\begin{tablenotes}[flushleft]",
"\\footnotesize",
paste0(
"\\item Notes: Entries are percentage effects, $100[\\exp(\\widehat{\\mathrm{effect}})-1]$, with standard errors in parentheses. ",
"The restricted exposure effect is evaluated as $\\hat\\beta^r\\bar S$, where $\\bar S$ is the mean positive benchmark exposure. ",
"The treated-exposure and currently untreated exposure effects are evaluated as $\\hat\\beta^t\\bar S^t$ and $\\hat\\beta^u\\bar S^u$, where $\\bar S^t$ and $\\bar S^u$ are the mean positive exposures among treated exposed and currently untreated exposed observations, respectively. ",
"Estimated exposure means are $\\bar S=", fmt_num(sbar_all, 3), "$, $\\bar S^t=", fmt_num(sbar_treated, 3), "$, and $\\bar S^u=", fmt_num(sbar_untreated, 3), "$. ",
"Significance: * $p<0.10$, ** $p<0.05$, *** $p<0.01$. Source: London Air Quality Network/Imperial College London monitoring data and Greater London Authority/Transport for London ULEZ boundary files. Author's calculations."
),
"\\end{tablenotes}",
"\\end{threeparttable}",
"\\end{table}"
)

writeLines(benchmark_lines, file.path(TAB, "table_ulez_benchmark_implied_effects.tex"))

# ---- 6. Conventional robustness checks -----------------------

m_level <- fixest::feols(
  no2 ~ D_active + DS_bench + US_bench | code + month,
  data = panel,
  cluster = ~code,
  warn = FALSE,
  notes = FALSE
)

m_central <- fixest::feols(
  log_y ~ D_active + DS_bench + US_bench | code + month,
  data = panel_central,
  cluster = ~code,
  warn = FALSE,
  notes = FALSE
)

m_twoway <- fixest::feols(
  log_y ~ D_active + DS_bench + US_bench | code + month,
  data = panel,
  cluster = ~code + month,
  warn = FALSE,
  notes = FALSE
)

robust_models <- list(
  baseline_log = m_unrestricted,
  level_no2 = m_level,
  central_only = m_central,
  twoway_cluster = m_twoway
)

saveRDS(
  list(
    no_exposure = m_noexp,
    restricted = m_restricted,
    unrestricted = m_unrestricted,
    level_no2 = m_level,
    central_only = m_central,
    twoway_cluster = m_twoway
  ),
  file.path(RDS, "ulez_appendix_models_v6.rds")
)

robust_specs <- tibble::tribble(
  ~model_key, ~model_label, ~outcome_scale, ~sample, ~inference, ~scale,
  "baseline_log", "Baseline", "Percent", "Full panel", "Station cluster", "log",
  "level_no2", "Level outcome", "$\\mu$g/m$^3$", "Full panel", "Station cluster", "level",
  "central_only", "Central ULEZ only", "Percent", "2016--Sep. 2021", "Station cluster", "log",
  "twoway_cluster", "Two-way cluster", "Percent", "Full panel", "Station and month cluster", "log"
)

# Keep data objects outside the tibble to avoid nested list-column issues.
robust_data <- list(
  baseline_log = panel,
  level_no2 = panel,
  central_only = panel_central,
  twoway_cluster = panel
)

robust_results <- purrr::map_dfr(seq_len(nrow(robust_specs)), function(k) {
  key <- robust_specs$model_key[k]
  mod <- robust_models[[key]]
  scale <- robust_specs$scale[k]
  dat <- robust_data[[key]]

  if (is.null(dat) || !is.data.frame(dat)) {
    stop("Robustness data object is not a data frame for key: ", key)
  }

  s_treated <- mean_pos_exposure(dat, "treated")
  s_untreated <- mean_pos_exposure(dat, "untreated")

  bind_rows(
    effect_cell(mod, "D_active", robust_specs$model_label[k], scale, exposure_scale = 1),
    effect_cell(mod, "DS_bench", robust_specs$model_label[k], scale, exposure_scale = s_treated),
    effect_cell(mod, "US_bench", robust_specs$model_label[k], scale, exposure_scale = s_untreated)
  ) |>
    mutate(
      model_key = key,
      outcome_scale = robust_specs$outcome_scale[k],
      sample = robust_specs$sample[k],
      inference = robust_specs$inference[k],
      sbar_treated = s_treated,
      sbar_untreated = s_untreated,
      nobs = safe_nobs(mod),
      within_r2 = safe_within_r2(mod),
      .before = model
    )
})

# Diagnostic: stop early if implied exposure scales are missing.
if (any(is.na(robust_results$exposure_scale[robust_results$term %in% c("DS_bench", "US_bench")]))) {
  warning("Some exposure scales are NA in robustness results. Check treatment-exposure support.")
}

readr::write_csv(robust_results, file.path(TAB, "table_ulez_conventional_robustness.csv"))

robust_wide <- robust_results |>
  mutate(term_label = dplyr::case_when(
    term == "D_active" ~ "beta_d",
    term == "DS_bench" ~ "beta_t_sbar_t",
    term == "US_bench" ~ "beta_u_sbar_u",
    TRUE ~ term
  )) |>
  select(model, outcome_scale, sample, inference, sbar_treated, sbar_untreated, nobs, within_r2, term_label, cell) |>
  tidyr::pivot_wider(names_from = term_label, values_from = cell) |>
  arrange(match(model, robust_specs$model_label))

for (nm in c("beta_d", "beta_t_sbar_t", "beta_u_sbar_u")) {
  if (!nm %in% names(robust_wide)) robust_wide[[nm]] <- ""
}

robust_lines <- c(
"\\begin{table}[!htbp]\\centering",
"\\caption{Conventional robustness checks for implied ULEZ exposure effects}",
"\\label{tab:ulez_conventional_robustness}",
"\\begin{threeparttable}",
"\\scriptsize",
"\\setlength{\\tabcolsep}{2pt}",
"\\renewcommand{\\arraystretch}{1.08}",
"\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}X>{\\RaggedRight\\arraybackslash}p{0.11\\textwidth}>{\\RaggedRight\\arraybackslash}p{0.13\\textwidth}>{\\RaggedRight\\arraybackslash}p{0.14\\textwidth}rrrrr@{}}",
"\\toprule",
"Specification & Scale & Sample & Inference & $\\bar S^t$ & $\\bar S^u$ & $\\beta^d$ & $\\beta^t\\bar S^t$ & $\\beta^u\\bar S^u$ \\\\",
"\\midrule",
paste0(
  robust_wide$model, " & ",
  robust_wide$outcome_scale, " & ",
  robust_wide$sample, " & ",
  robust_wide$inference, " & ",
  fmt_num(robust_wide$sbar_treated, 3), " & ",
  fmt_num(robust_wide$sbar_untreated, 3), " & ",
  robust_wide$beta_d, " & ",
  robust_wide$beta_t_sbar_t, " & ",
  robust_wide$beta_u_sbar_u, " \\\\"
),
"\\bottomrule",
"\\end{tabularx}",
"\\begin{tablenotes}[flushleft]",
"\\footnotesize",
"\\item Notes: The table reports implied effects from the unrestricted model with the benchmark 0--3 km exposure mapping. $\\bar S^t$ is the mean positive exposure among treated exposed observations, and $\\bar S^u$ is the mean positive exposure among untreated exposed observations in the relevant estimation sample. Exposure effects are evaluated as $\\hat\\beta^t\\bar S^t$ and $\\hat\\beta^u\\bar S^u$, rather than as unit-exposure slopes. For log-outcome specifications, entries are percentage effects, $100[\\exp(\\widehat{\\mathrm{effect}})-1]$, with standard errors in parentheses computed by the delta method treating $\\bar S$ as fixed. For the level-outcome specification, entries are in $\\mu$g/m$^3$. The central-ULEZ-only specification uses observations through September 2021 to exclude the Inner London and London-wide expansions. Significance: * $p<0.10$, ** $p<0.05$, *** $p<0.01$. Source: London Air Quality Network/Imperial College London monitoring data and Greater London Authority/Transport for London ULEZ boundary files. Author's calculations.",
"\\end{tablenotes}",
"\\end{threeparttable}",
"\\end{table}"
)

writeLines(robust_lines, file.path(TAB, "table_ulez_conventional_robustness.tex"))

# ---- 7. Appendix section snippet -----------------------------

appendix_lines <- c(
"\\section{Additional diagnostics for the ULEZ illustration}",
"\\label{app:ulez_additional_diagnostics}",
"",
"This appendix reports descriptive statistics and conventional robustness checks for the ULEZ empirical illustration.",
"The purpose is to document the empirical design without turning the main illustration into a standalone policy evaluation of London's ULEZ.",
"Changing the exposure matrix is treated in the main text as a mechanism check because it changes the exposure definition and potentially the estimand.",
"The checks below instead hold the benchmark 0--3 km exposure mapping fixed and vary the outcome scale, sample window, and inference procedure.",
"",
"\\input{table_ulez_descriptive_statistics.tex}",
"",
"\\input{table_ulez_benchmark_implied_effects.tex}",
"",
"\\input{table_ulez_conventional_robustness.tex}",
"",
"The benchmark and robustness tables report implied exposure effects, not only exposure slopes.",
"For each specification, the treated-exposure effect is evaluated as $\\hat\\beta^t\\bar S^t$, where $\\bar S^t$ is the mean positive exposure among treated exposed observations.",
"The currently untreated exposure effect is evaluated as $\\hat\\beta^u\\bar S^u$, where $\\bar S^u$ is the mean positive exposure among untreated exposed observations.",
"This convention matches the definition $S_{it}(\\mathbf{W})=\\sum_jw_{ij}D_{jt}$ and avoids interpreting exposure slopes as if all exposed observations had unit exposure.",
"",
"The central-ULEZ-only subsample should be interpreted with caution for the direct and treated-exposure components.",
"Restricting the sample to the first ULEZ episode reduces the amount of identifying variation for $\\beta^d$ and $\\beta^t$, which leads to substantially imprecise estimates for these coefficients.",
"By contrast, the currently untreated exposure component remains negative and more precisely estimated in this subsample.",
"This pattern is consistent with the support diagnostics in the main text: the empirical illustration provides stronger support for the currently untreated exposure component than for a precise decomposition between the direct-treatment and treated-exposure components in restricted windows."
)

writeLines(appendix_lines, file.path(TAB, "appendix_ulez_descriptives_robustness.tex"))

# ---- 8. Console report ---------------------------------------

report_lines <- c(
  "ULEZ appendix descriptives and implied-effect outputs",
  "====================================================",
  paste0("Panel observations: ", nrow(panel)),
  paste0("Panel stations: ", dplyr::n_distinct(panel$code)),
  paste0("Panel months: ", dplyr::n_distinct(panel$month)),
  paste0("Full-panel mean positive exposure: ", fmt_num(sbar_all, 4)),
  paste0("Full-panel mean positive treated exposure: ", fmt_num(sbar_treated, 4)),
  paste0("Full-panel mean positive untreated exposure: ", fmt_num(sbar_untreated, 4)),
  paste0("Central-only observations: ", nrow(panel_central)),
  paste0("Central-only mean positive treated exposure: ", fmt_num(mean_pos_exposure(panel_central, "treated"), 4)),
  paste0("Central-only mean positive untreated exposure: ", fmt_num(mean_pos_exposure(panel_central, "untreated"), 4)),
  "",
  "Generated files:",
  file.path(TAB, "table_ulez_descriptive_statistics.tex"),
  file.path(TAB, "table_ulez_descriptive_statistics.csv"),
  file.path(TAB, "table_ulez_benchmark_implied_effects.tex"),
  file.path(TAB, "table_ulez_benchmark_implied_effects.csv"),
  file.path(TAB, "table_ulez_conventional_robustness.tex"),
  file.path(TAB, "table_ulez_conventional_robustness.csv"),
  file.path(TAB, "appendix_ulez_descriptives_robustness.tex")
)

writeLines(report_lines, file.path(TAB, "appendix_ulez_descriptives_robustness_report.txt"))
cat(paste(report_lines, collapse = "\n"), "\n")

# Copiar para outputs/ apenas as três tabelas que aparecem no manuscrito.
for (ff in c(
  "table_ulez_descriptive_statistics.tex",
  "table_ulez_benchmark_implied_effects.tex",
  "table_ulez_conventional_robustness.tex"
)) {
  ok <- file.copy(file.path(TAB, ff), file.path(TABLE_DIR, ff), overwrite = TRUE)
  if (!ok) stop("Falha ao exportar: ", ff)
}

# ============================================================
# 8. Spatial inference diagnostics
# ============================================================
# O ponto estimado permanece fixo; apenas a matriz de covariância muda.
# Coordinate preparation mirrors 11-build_ulez_spatial_inference_table_v6.R.
coord_names_present <- any(c("latitude", "lat") %in% names(panel)) &&
  any(c("longitude", "lon") %in% names(panel))

if (!coord_names_present) {
  if (!file.exists(ASSIGN_FILE)) {
    stop("Coordinate variables not found in panel and assignment file is missing: ", ASSIGN_FILE)
  }

  station_coord <- readr::read_csv(ASSIGN_FILE, show_col_types = FALSE) |>
    janitor::clean_names() |>
    dplyr::mutate(
      code = stringr::str_to_upper(as.character(code)),
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude)
    ) |>
    dplyr::select(code, longitude, latitude) |>
    dplyr::filter(is.finite(longitude), is.finite(latitude)) |>
    dplyr::distinct(code, .keep_all = TRUE)

  panel <- panel |>
    dplyr::left_join(station_coord, by = "code")
}

if ("lat" %in% names(panel) && !"latitude" %in% names(panel)) {
  panel <- panel |> dplyr::rename(latitude = lat)
}
if ("lon" %in% names(panel) && !"longitude" %in% names(panel)) {
  panel <- panel |> dplyr::rename(longitude = lon)
}

if (!all(c("latitude", "longitude") %in% names(panel))) {
  stop("Spatial inference requires latitude and longitude.")
}
if (any(!is.finite(panel$latitude)) || any(!is.finite(panel$longitude))) {
  stop("Spatial inference found non-finite station coordinates after joining the assignment file.")
}

m_inf <- fixest::feols(
  log_y ~ D_active + DS_bench + US_bench | code + month,
  data = panel,
  cluster = ~ code,
  warn = FALSE,
  notes = FALSE
)

extract_with_vcov <- function(model, vc, inference_label) {
  ct <- as.data.frame(fixest::coeftable(model, vcov = vc))
  keep <- intersect(c("D_active", "DS_bench", "US_bench"), rownames(ct))
  out <- tibble::tibble(
    inference = inference_label,
    term = keep,
    estimate = ct[keep, 1],
    std_error = ct[keep, 2],
    statistic = ct[keep, 3],
    p_value = ct[keep, 4]
  ) |>
    dplyr::mutate(
      parameter = dplyr::recode(
        term,
        D_active = "$\\beta^d$",
        DS_bench = "$\\beta^t$",
        US_bench = "$\\beta^u$"
      ),
      effect_pct = pct_effect(estimate),
      se_effect_pct = 100 * exp(estimate) * std_error,
      cell = paste0(fmt_num(effect_pct, 2), stars(p_value), " (", fmt_num(se_effect_pct, 2), ")")
    )
  out
}

spatial_rows <- list()
spatial_rows[[1]] <- extract_with_vcov(m_inf, ~ code, "Station cluster")
spatial_rows[[2]] <- extract_with_vcov(m_inf, ~ code + month, "Station x month cluster")

# Conley-type HAC. Coordinates are supplied explicitly and the cutoff is in km.
for (cc in c(3, 5, 10)) {
  vc_req <- fixest::vcov_conley(
    lat = "latitude",
    lon = "longitude",
    cutoff = cc,
    distance = "triangular"
  )
  spatial_rows[[length(spatial_rows) + 1L]] <- extract_with_vcov(
    m_inf, vc_req, paste0("Conley ", cc, " km")
  )
}

# Spatial block bootstrap. Eight geographic station blocks are formed by k-means
# on scaled longitude/latitude coordinates. For prospective reproducibility, the
# k-means step and the bootstrap resampling use separate, explicit RNG seeds.
# Cluster labels are canonicalized by centroid location so arbitrary k-means label
# permutations cannot change the finite bootstrap draw sequence.
station_blocks <- panel |>
  dplyr::distinct(code, longitude, latitude) |>
  dplyr::filter(is.finite(longitude), is.finite(latitude))

if (nrow(station_blocks) < 8L) {
  stop("Spatial block bootstrap requires at least 8 stations with valid coordinates.")
}

KMEANS_SEED <- 20260615L
BOOTSTRAP_SEED <- 20260615L
B_BOOT <- 499L

set.seed(KMEANS_SEED)
km8 <- stats::kmeans(
  scale(station_blocks |> dplyr::select(longitude, latitude)),
  centers = 8,
  nstart = 50,
  algorithm = "Hartigan-Wong"
)

# Canonicalize cluster labels by centroid longitude and then latitude.
centers_raw <- vapply(seq_len(8L), function(k) {
  idx <- km8$cluster == k
  c(
    longitude = mean(station_blocks$longitude[idx]),
    latitude = mean(station_blocks$latitude[idx])
  )
}, numeric(2)) |>
  t() |>
  as.data.frame()
centers_raw$old_cluster <- seq_len(8L)
centers_raw <- centers_raw |>
  dplyr::arrange(longitude, latitude) |>
  dplyr::mutate(spatial_block = dplyr::row_number())
cluster_key <- stats::setNames(centers_raw$spatial_block, centers_raw$old_cluster)

block_map <- station_blocks |>
  dplyr::mutate(spatial_block = as.integer(cluster_key[as.character(km8$cluster)]))

panel_boot <- panel |>
  dplyr::left_join(block_map |> dplyr::select(code, spatial_block), by = "code")

boot_terms <- c("D_active", "DS_bench", "US_bench")
boot_coef <- matrix(NA_real_, nrow = B_BOOT, ncol = length(boot_terms),
                    dimnames = list(NULL, boot_terms))
blocks <- sort(unique(panel_boot$spatial_block))
if (!identical(blocks, seq_len(8L))) {
  stop("Spatial block bootstrap did not produce canonical block labels 1:8.")
}

# Use an independent bootstrap RNG stream so resampling does not depend on the
# number of random values consumed internally by k-means. The target is 499
# valid bootstrap replications. Rare resamples in which one or more target
# coefficients are not identified are discarded and deterministically replaced
# by the next draw from the same RNG stream.
set.seed(BOOTSTRAP_SEED)
MAX_BOOT_ATTEMPTS <- 10L * B_BOOT
n_valid_boot <- 0L
n_attempted_boot <- 0L

while (n_valid_boot < B_BOOT && n_attempted_boot < MAX_BOOT_ATTEMPTS) {
  n_attempted_boot <- n_attempted_boot + 1L
  sampled <- sample(blocks, length(blocks), replace = TRUE)
  pieces <- lapply(seq_along(sampled), function(j) {
    x <- panel_boot[panel_boot$spatial_block == sampled[j], , drop = FALSE]
    # IDs recebem sufixo para permitir a repetição de um mesmo bloco no bootstrap.
    x$code_boot <- paste0(x$code, "__draw", j)
    x
  })
  boot_dat <- dplyr::bind_rows(pieces)
  mb <- tryCatch(
    fixest::feols(
      log_y ~ D_active + DS_bench + US_bench | code_boot + month,
      data = boot_dat,
      warn = FALSE, notes = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(mb)) {
    b <- stats::coef(mb)
    if (all(boot_terms %in% names(b))) {
      b_target <- as.numeric(b[boot_terms])
      if (all(is.finite(b_target))) {
        n_valid_boot <- n_valid_boot + 1L
        boot_coef[n_valid_boot, ] <- b_target
      }
    }
  }
}

if (n_valid_boot != B_BOOT) {
  stop(
    "Spatial block bootstrap obtained only ", n_valid_boot,
    " valid replications after ", n_attempted_boot,
    " attempts; expected ", B_BOOT, "."
  )
}

boot_se <- apply(boot_coef, 2, stats::sd)

base_b <- stats::coef(m_inf)[c("D_active", "DS_bench", "US_bench")]
boot_stat <- base_b / boot_se
boot_p <- 2 * stats::pnorm(-abs(boot_stat))
boot_row <- tibble::tibble(
  inference = "Spatial block bootstrap (8 blocks)",
  term = names(base_b), estimate = as.numeric(base_b),
  std_error = as.numeric(boot_se), statistic = as.numeric(boot_stat), p_value = as.numeric(boot_p)
) |>
  dplyr::mutate(
    parameter = dplyr::recode(term, D_active = "$\\beta^d$", DS_bench = "$\\beta^t$", US_bench = "$\\beta^u$"),
    effect_pct = pct_effect(estimate),
    se_effect_pct = 100 * exp(estimate) * std_error,
    cell = paste0(fmt_num(effect_pct, 2), stars(p_value), " (", fmt_num(se_effect_pct, 2), ")")
  )
spatial_inference <- dplyr::bind_rows(spatial_rows, list(boot_row))

# Tabela compacta.
inf_order <- c("Station cluster", "Station x month cluster", "Conley 3 km", "Conley 5 km", "Conley 10 km", "Spatial block bootstrap (8 blocks)")
inf_wide <- spatial_inference |>
  dplyr::mutate(
    inference = factor(inference, levels = inf_order),
    parameter = factor(parameter, levels = c("$\\beta^d$", "$\\beta^t$", "$\\beta^u$"))
  ) |>
  dplyr::select(inference, parameter, cell) |>
  tidyr::pivot_wider(names_from = parameter, values_from = cell) |>
  dplyr::arrange(inference)
writeLines(c(
  "\\begin{table}[!htbp]", "\\centering", "\\begin{threeparttable}",
  "\\caption{Inference diagnostics for the benchmark unrestricted exposure model}",
  "\\label{tab:ulez_spatial_inference}", "\\begin{tabular}{lccc}", "\\toprule",
  "Inference procedure & $\\beta^d$ & $\\beta^t$ & $\\beta^u$ \\\\", "\\midrule",
  paste0(as.character(inf_wide$inference), " & ", inf_wide[["$\\beta^d$"]], " & ", inf_wide[["$\\beta^t$"]], " & ", inf_wide[["$\\beta^u$"]], " \\\\"),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: Entries are percentage effects, $100[\\exp(\\widehat\\beta)-1]$, with standard errors in parentheses computed by the delta method. $\\beta^d$ is the direct-treatment parameter, $\\beta^t$ the treated-exposure parameter, and $\\beta^u$ the currently untreated exposure parameter. Conley rows hold the 0--3 km exposure mapping fixed and vary only the covariance estimator. The spatial block bootstrap uses eight geographic station blocks.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_spatial_inference.tex"))

# ============================================================
# 9. Exposure-response bins
# ============================================================
positive_s <- panel$S_bench[panel$S_bench > EXPOSURE_EPS & is.finite(panel$S_bench)]
if (length(positive_s) < 20L) stop("Too few positive-exposure observations to construct bins.")

# Original manuscript rule (12-build_ulez_exposure_bins_linearity_table.R):
# start from positive-exposure quartiles; if ties make the five quartile cutpoints
# non-unique, fall back to four equal-width bins over the positive-exposure range.
q_raw <- stats::quantile(
  positive_s,
  probs = c(0, 0.25, 0.50, 0.75, 1),
  na.rm = TRUE,
  type = 7
)
q_breaks <- as.numeric(q_raw)
if (length(unique(q_breaks)) < 5L) {
  q_breaks <- seq(
    min(positive_s, na.rm = TRUE),
    max(positive_s, na.rm = TRUE),
    length.out = 5L
  )
}

panel_bins <- panel |>
  dplyr::mutate(
    exposure_bin = dplyr::case_when(
      S_bench <= EXPOSURE_EPS ~ "S0",
      S_bench > EXPOSURE_EPS ~ as.character(cut(
        S_bench,
        breaks = q_breaks,
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE,
        right = TRUE
      )),
      TRUE ~ NA_character_
    ),
    exposure_bin = factor(exposure_bin, levels = c("S0", paste0("Q", 1:4)))
  )
for (qq in paste0("Q", 1:4)) {
  panel_bins[[paste0("T_", qq)]] <- as.integer(panel_bins$D_active == 1 & panel_bins$exposure_bin == qq)
  panel_bins[[paste0("U_", qq)]] <- as.integer(panel_bins$D_active == 0 & panel_bins$exposure_bin == qq)
}

bin_terms <- c(paste0("T_Q", 1:4), paste0("U_Q", 1:4))
if (anyNA(panel_bins[, bin_terms])) {
  stop("Exposure-bin indicators contain NA values; zero exposure must be retained as the S0 reference category.")
}

bin_rhs <- paste(c("D_active", bin_terms), collapse = " + ")
mod_bins <- fixest::feols(
  stats::as.formula(paste0("log_y ~ ", bin_rhs, " | code + month")),
  data = panel_bins, cluster = ~ code, warn = FALSE, notes = FALSE
)
if (as.numeric(stats::nobs(mod_bins)) != nrow(panel_bins)) {
  stop(
    "Exposure-bin model used ", as.numeric(stats::nobs(mod_bins)),
    " observations; expected ", nrow(panel_bins),
    ". Check missing values in the bin indicators."
  )
}
ct_bin <- as.data.frame(fixest::coeftable(mod_bins))
make_bin_cell <- function(term) {
  if (!term %in% rownames(ct_bin)) return(NA_character_)
  est <- ct_bin[term, 1]; se <- ct_bin[term, 2]; p <- ct_bin[term, 4]
  paste0(fmt_num(pct_effect(est), 2), stars(p), " (", fmt_num(100 * exp(est) * se, 2), ")")
}
bin_support <- panel_bins |>
  dplyr::filter(!is.na(exposure_bin)) |>
  dplyr::count(exposure_bin, D_active, name = "n") |>
  tidyr::pivot_wider(names_from = D_active, values_from = n, values_fill = 0, names_prefix = "D")
bin_table <- tibble::tibble(
  bin_internal = paste0("Q", 1:4),
  treated = vapply(paste0("T_Q", 1:4), make_bin_cell, character(1)),
  untreated = vapply(paste0("U_Q", 1:4), make_bin_cell, character(1))
) |>
  dplyr::left_join(bin_support, by = c("bin_internal" = "exposure_bin")) |>
  dplyr::rename(n_untreated = D0, n_treated = D1) |>
  dplyr::mutate(bin = paste0("B", dplyr::row_number()), .before = 1) |>
  dplyr::select(-bin_internal)

writeLines(c(
  "\\begin{table}[!htbp]", "\\centering", "\\begin{threeparttable}",
  "\\caption{Diagnostic for the linear exposure-response restriction}",
  "\\label{tab:ulez_exposure_bins}", "\\begin{tabular}{lrcrc}", "\\toprule",
  "Exposure bin & Treated obs. & Treated exposure & Untreated obs. & Currently untreated exposure \\\\", "\\midrule",
  paste0(
    bin_table$bin, " & ", bin_table$n_treated, " & ",
    ifelse(is.na(bin_table$treated), "", bin_table$treated), " & ",
    bin_table$n_untreated, " & ",
    ifelse(is.na(bin_table$untreated), "", bin_table$untreated), " \\\\"
  ),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: The continuous 0--3 km exposure index is replaced by four indicators among observations with positive exposure. The implementation first attempts quartile cutpoints; because tied positive-exposure quantiles make those cutpoints non-unique in this sample, the reported diagnostic uses the prespecified fallback of four equal-width bins over the positive-exposure range. Treated-exposure coefficients are associated with $D_{it}1\\{S_{it}\\in B_k\\}$; currently untreated exposure coefficients are associated with $(1-D_{it})1\\{S_{it}\\in B_k\\}$. Entries are percentage effects, $100[\\exp(\\widehat\\theta)-1]$, with station-clustered standard errors in parentheses. The model includes monitoring-station and month fixed effects.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_exposure_bins.tex"))

saveRDS(
  list(spatial_inference = spatial_inference, exposure_bins = bin_table, exposure_bin_model = mod_bins),
  file.path(DERIVED_DIR, "appendix_additional_results.rds")
)
