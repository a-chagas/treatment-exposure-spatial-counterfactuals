# ============================================================
# 02_static_results.R
# Resultados estáticos, suporte, pre-trends e comparação de W
# ============================================================

message("[02] Estimando especificações estáticas...")
if (!file.exists(STATIC_PANEL_FILE)) stop("Execute 01_build_static_exposure.R primeiro.")

panel <- readr::read_csv(STATIC_PANEL_FILE, show_col_types = FALSE) |>
  janitor::clean_names() |>
  dplyr::mutate(
    code = stringr::str_to_upper(as.character(code)),
    month = as.Date(month),
    d_active = as.integer(d_active),
    s_bench = s_bench_0_3km,
    ds_bench = d_active * s_bench,
    us_bench = (1 - d_active) * s_bench,
    s_pos_bench = as.integer(s_bench > EXPOSURE_EPS)
  ) |>
  dplyr::filter(is.finite(log_y), !is.na(code), !is.na(month), !is.na(d_active))

matrix_specs <- tibble::tribble(
  ~matrix_name,       ~s_col,                ~matrix_label,
  "bench_0_2km",     "s_bench_0_2km",      "Threshold 0--2 km",
  "bench_0_3km",     "s_bench_0_3km",      "Threshold 0--3 km",
  "bench_0_5km",     "s_bench_0_5km",      "Threshold 0--5 km",
  "near_0_2km",      "s_near_0_2km",       "Ring 0--2 km",
  "near_2_5km",      "s_near_2_5km",       "Ring 2--5 km",
  "medium_5_10km",   "s_medium_5_10km",    "Ring 5--10 km",
  "far_10_20km",     "s_far_10_20km",      "Ring 10--20 km",
  "far_20_40km",     "s_far_20_40km",      "Ring 20--40 km"
)

extract_coef <- function(mod, term) {
  ct <- as.data.frame(fixest::coeftable(mod))
  if (!term %in% rownames(ct)) {
    return(tibble::tibble(estimate = NA_real_, std_error = NA_real_, p_value = NA_real_))
  }
  tibble::tibble(
    estimate = ct[term, 1],
    std_error = ct[term, 2],
    p_value = ct[term, 4]
  )
}

support_one <- function(s_col, matrix_name, matrix_label) {
  s <- panel[[s_col]]
  tibble::tibble(
    matrix_name = matrix_name,
    matrix_label = matrix_label,
    d0_s0 = sum(panel$d_active == 0 & s <= EXPOSURE_EPS, na.rm = TRUE),
    d0_s1 = sum(panel$d_active == 0 & s > EXPOSURE_EPS, na.rm = TRUE),
    d1_s0 = sum(panel$d_active == 1 & s <= EXPOSURE_EPS, na.rm = TRUE),
    d1_s1 = sum(panel$d_active == 1 & s > EXPOSURE_EPS, na.rm = TRUE)
  )
}

fit_one_matrix <- function(s_col, matrix_name, matrix_label) {
  dat <- panel |>
    dplyr::mutate(
      s_i = .data[[s_col]],
      ds_i = d_active * s_i,
      us_i = (1 - d_active) * s_i
    )

  m_u <- fixest::feols(
    log_y ~ d_active + ds_i + us_i | code + month,
    data = dat, cluster = ~ code, warn = FALSE, notes = FALSE
  )
  m_r <- fixest::feols(
    log_y ~ d_active + s_i | code + month,
    data = dat, cluster = ~ code, warn = FALSE, notes = FALSE
  )

  cd <- extract_coef(m_u, "d_active")
  ct <- extract_coef(m_u, "ds_i")
  cu <- extract_coef(m_u, "us_i")
  cr <- extract_coef(m_r, "s_i")
  sup <- support_one(s_col, matrix_name, matrix_label)

  b <- stats::coef(m_u)
  V <- stats::vcov(m_u)
  diff <- unname(b["ds_i"] - b["us_i"])
  se_diff <- sqrt(max(unname(V["ds_i", "ds_i"] + V["us_i", "us_i"] - 2 * V["ds_i", "us_i"]), 0))
  t_diff <- diff / se_diff
  p_equal <- 2 * stats::pt(abs(t_diff), df = data.table::uniqueN(panel$code) - 1L, lower.tail = FALSE)

  row <- dplyr::bind_cols(
    sup,
    tibble::tibble(
      beta_d = cd$estimate, beta_d_se = cd$std_error, beta_d_p_value = cd$p_value,
      beta_t = ct$estimate, beta_t_se = ct$std_error, beta_t_p_value = ct$p_value,
      beta_u = cu$estimate, beta_u_se = cu$std_error, beta_u_p_value = cu$p_value,
      beta_r = cr$estimate, beta_r_se = cr$std_error, beta_r_p_value = cr$p_value,
      p_value_equal_bt_bu = p_equal,
      bic_unrestricted = BIC(m_u),
      bic_restricted = BIC(m_r),
      delta_bic_df_minus_unrestricted = BIC(m_r) - BIC(m_u)
    )
  ) |>
    dplyr::mutate(
      beta_d_effect_pct = pct_effect(beta_d),
      beta_d_conf_low_pct = pct_effect(beta_d - 1.96 * beta_d_se),
      beta_d_conf_high_pct = pct_effect(beta_d + 1.96 * beta_d_se),
      beta_t_effect_pct = pct_effect(beta_t),
      beta_t_conf_low_pct = pct_effect(beta_t - 1.96 * beta_t_se),
      beta_t_conf_high_pct = pct_effect(beta_t + 1.96 * beta_t_se),
      beta_u_effect_pct = pct_effect(beta_u),
      beta_u_conf_low_pct = pct_effect(beta_u - 1.96 * beta_u_se),
      beta_u_conf_high_pct = pct_effect(beta_u + 1.96 * beta_u_se),
      beta_r_effect_pct = pct_effect(beta_r),
      beta_r_conf_low_pct = pct_effect(beta_r - 1.96 * beta_r_se),
      beta_r_conf_high_pct = pct_effect(beta_r + 1.96 * beta_r_se)
    )

  list(row = row, unrestricted = m_u, restricted = m_r)
}

fits <- lapply(seq_len(nrow(matrix_specs)), function(i) {
  fit_one_matrix(matrix_specs$s_col[i], matrix_specs$matrix_name[i], matrix_specs$matrix_label[i])
})
compare <- dplyr::bind_rows(lapply(fits, `[[`, "row"))
support <- compare |>
  dplyr::select(matrix_name, matrix_label, d0_s0, d0_s1, d1_s0, d1_s1)

# Modelos benchmark sem exposição / restrito / irrestrito.
m_no <- fixest::feols(log_y ~ d_active | code + month, data = panel, cluster = ~ code,
                      warn = FALSE, notes = FALSE)
m_res <- fixest::feols(log_y ~ d_active + s_bench | code + month, data = panel, cluster = ~ code,
                       warn = FALSE, notes = FALSE)
m_unres <- fixest::feols(log_y ~ d_active + ds_bench + us_bench | code + month, data = panel,
                         cluster = ~ code, warn = FALSE, notes = FALSE)

extract_model <- function(mod, model_name) {
  ct <- as.data.frame(fixest::coeftable(mod))
  tibble::tibble(
    model_name = model_name,
    term = rownames(ct),
    estimate = ct[, 1],
    std_error = ct[, 2],
    p_value = ct[, 4]
  ) |>
    dplyr::mutate(
      effect_pct = pct_effect(estimate),
      se_effect_approx_pct = 100 * exp(estimate) * std_error,
      cell = paste0(fmt_num(effect_pct, 2), stars(p_value), " (", fmt_num(se_effect_approx_pct, 2), ")")
    )
}
coef_all <- dplyr::bind_rows(
  extract_model(m_no, "No exposure"),
  extract_model(m_res, "Restricted exposure"),
  extract_model(m_unres, "Unrestricted exposure")
)
gof_all <- tibble::tibble(
  model_name = c("No exposure", "Restricted exposure", "Unrestricted exposure"),
  nobs = c(stats::nobs(m_no), stats::nobs(m_res), stats::nobs(m_unres)),
  bic = c(BIC(m_no), BIC(m_res), BIC(m_unres)),
  within_r2 = c(
    as.numeric(fixest::fitstat(m_no, "wr2")),
    as.numeric(fixest::fitstat(m_res, "wr2")),
    as.numeric(fixest::fitstat(m_unres, "wr2"))
  )
)

# ------------------------------------------------------------
# Tabelas do texto principal
# ------------------------------------------------------------
protocol_rows <- tibble::tribble(
  ~step, ~protocol_item, ~implementation,
  "1", "Treatment, unit, and timing", "Monthly monitoring-station panel; treatment equals being inside an active ULEZ boundary. Central ULEZ starts in April 2019, Inner London expansion in October 2021, and London-wide expansion in August 2023.",
  "2", "Propagation mechanism", "Vehicular-emissions regulation may affect nearby untreated stations through local pollution diffusion and traffic reallocation.",
  "3", "Choose $\\mathbf{W}$ from the mechanism", "Benchmark $\\mathbf{W}$ is a binary 0--3 km threshold. Alternatives include shorter thresholds and distance rings: 0--2, 2--5, 5--10, 10--20, and 20--40 km.",
  "4", "Construct treatment-exposure cells", "For each matrix, construct $S_{it}(\\mathbf{W})=\\sum_j w_{ij}D_{jt}$ and classify station-months into $D0/S0$, $D0/S1$, $D1/S0$, and $D1/S1$.",
  "5", "Report support and overlap", paste0(
    "Support diagnostics are reported before estimates. For the 0--3 km benchmark, support is $D0/S0=",
    compare$d0_s0[compare$matrix_name == "bench_0_3km"], "$; $D0/S1=",
    compare$d0_s1[compare$matrix_name == "bench_0_3km"], "$; $D1/S0=",
    compare$d1_s0[compare$matrix_name == "bench_0_3km"], "$; and $D1/S1=",
    compare$d1_s1[compare$matrix_name == "bench_0_3km"], "$ ."
  ),
  "6", "Pre-trends by treatment-exposure cell", "A benchmark-cell pre-trend figure and a pre-period slope test compare $D0/S1$, $D1/S0$, and $D1/S1$ with the $D0/S0$ reference group before April 2019.",
  "7", "Estimate direct-treatment, treated-exposure, and currently untreated exposure components", "The unrestricted model estimates $\\beta^d$, $\\beta^t$, and $\\beta^u$ separately. Restricted and no-exposure models are reported only for comparison.",
  "8", "Inference robust to dependence", "Baseline inference clusters standard errors by monitoring station. Serial/spatial alternatives are reported as appendix robustness checks.",
  "9", "Challenge $\\mathbf{W}$", "Alternative threshold and ring matrices are interpreted as competing exposure hypotheses, not as repeated estimates of one invariant estimand.",
  "10", "Total effects after decomposition", "The illustration emphasizes decomposed effects. The restricted common-exposure coefficient is interpreted only after $\\beta^t$ and $\\beta^u$ are reported."
)
writeLines(c(
  "\\begin{table}[!htbp]\\centering",
  "\\caption{Implementation of the exposure-mapping protocol in the London ULEZ illustration}",
  "\\label{tab:ulez_protocol_implementation}",
  "\\begin{threeparttable}", "\\small", "\\setlength{\\tabcolsep}{3pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "\\begin{tabularx}{\\textwidth}{@{}c>{\\RaggedRight\\arraybackslash}p{0.25\\textwidth}>{\\RaggedRight\\arraybackslash}X@{}}",
  "\\toprule", "Step & Protocol item & Implementation in the ULEZ illustration \\\\", "\\midrule",
  paste0(protocol_rows$step, " & ", protocol_rows$protocol_item, " & ", protocol_rows$implementation, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: The table maps the ten-step protocol into the empirical illustration. The exercise illustrates exposure mapping rather than providing a complete policy evaluation of the ULEZ.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_protocol_implementation.tex"))

# Benchmark model comparison.
lookup <- coef_all |>
  dplyr::mutate(term_label = dplyr::case_when(
    term == "d_active" ~ "beta_d",
    term == "s_bench" ~ "beta_r",
    term == "ds_bench" ~ "beta_t",
    term == "us_bench" ~ "beta_u",
    TRUE ~ term
  ))
expected_lookup <- tibble::tribble(
  ~model_name,               ~term_label,
  "No exposure",             "beta_d",
  "Restricted exposure",     "beta_d",
  "Restricted exposure",     "beta_r",
  "Unrestricted exposure",   "beta_d",
  "Unrestricted exposure",   "beta_t",
  "Unrestricted exposure",   "beta_u"
)

missing_lookup <- expected_lookup |>
  dplyr::anti_join(
    lookup |>
      dplyr::distinct(model_name, term_label),
    by = c("model_name", "term_label")
  )

if (nrow(missing_lookup) > 0L) {
  stop(
    "Missing coefficients required for ",
    "table_ulez_no_restricted_unrestricted.tex: ",
    paste(
      paste0(missing_lookup$model_name, " / ", missing_lookup$term_label),
      collapse = "; "
    )
  )
}
get_cell <- function(model, target_term) {
  x <- lookup |> dplyr::filter(.data$model_name == .env$model, .data$term_label == .env$target_term) |> dplyr::pull(.data$cell)
  if (!length(x)) "" else x[1]
}
get_gof <- function(model, var) {
  x <- gof_all |> dplyr::filter(model_name == model) |> dplyr::pull(.data[[var]])
  if (!length(x)) return("")
  if (var == "nobs") return(as.character(x[1]))
  fmt_num(x[1], ifelse(var == "bic", 1, 3))
}
bench_compare <- compare |> dplyr::filter(matrix_name == "bench_0_3km") |> dplyr::slice_head(n = 1)
writeLines(c(
  "\\begin{table}[!htbp]\\centering",
  "\\caption{Benchmark comparison of no-exposure, restricted-exposure, and unrestricted-exposure models}",
  "\\label{tab:ulez_no_restricted_unrestricted}", "\\begin{threeparttable}", "\\small",
  "\\setlength{\\tabcolsep}{3pt}", "\\renewcommand{\\arraystretch}{1.10}",
  "\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}Xccc@{}}",
  "\\toprule", " & No exposure & Restricted exposure & Unrestricted exposure \\\\", "\\midrule",
  paste0("$\\beta^d D_{it}$ & ", get_cell("No exposure", "beta_d"), " & ", get_cell("Restricted exposure", "beta_d"), " & ", get_cell("Unrestricted exposure", "beta_d"), " \\\\"),
  paste0("$\\beta^r S_{it}$ &  & ", get_cell("Restricted exposure", "beta_r"), " &  \\\\"),
  paste0("$\\beta^t D_{it}S_{it}$ &  &  & ", get_cell("Unrestricted exposure", "beta_t"), " \\\\"),
  paste0("$\\beta^u(1-D_{it})S_{it}$ &  &  & ", get_cell("Unrestricted exposure", "beta_u"), " \\\\"),
  "\\midrule",
  paste0("$H_0:\\beta^t=\\beta^u$ p-value &  &  & ", fmt_p(bench_compare$p_value_equal_bt_bu), " \\\\"),
  paste0("$\\Delta$BIC: restricted $-$ unrestricted &  & ", fmt_num(bench_compare$delta_bic_df_minus_unrestricted, 1), " &  \\\\"),
  "\\midrule",
  paste0("Observations & ", get_gof("No exposure", "nobs"), " & ", get_gof("Restricted exposure", "nobs"), " & ", get_gof("Unrestricted exposure", "nobs"), " \\\\"),
  paste0("Within $R^2$ & ", get_gof("No exposure", "within_r2"), " & ", get_gof("Restricted exposure", "within_r2"), " & ", get_gof("Unrestricted exposure", "within_r2"), " \\\\"),
  paste0("BIC & ", get_gof("No exposure", "bic"), " & ", get_gof("Restricted exposure", "bic"), " & ", get_gof("Unrestricted exposure", "bic"), " \\\\"),
  "Station FE & Yes & Yes & Yes \\\\", "Month FE & Yes & Yes & Yes \\\\",
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: Entries are percentage effects, $100[\\exp(\\hat\\beta)-1]$, with cluster-robust standard errors by monitoring station in parentheses. The benchmark exposure mapping is the binary 0--3 km threshold matrix. The restricted model imposes $\\beta^t=\\beta^u=\\beta^r$. Significance: * $p<0.10$, ** $p<0.05$, *** $p<0.01$.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_no_restricted_unrestricted.tex"))

# Support table.
support_main <- support |>
  dplyr::filter(matrix_name %in% c("bench_0_3km", "near_2_5km", "far_10_20km")) |>
  dplyr::mutate(matrix_name = factor(matrix_name, levels = c("bench_0_3km", "near_2_5km", "far_10_20km"))) |>
  dplyr::arrange(matrix_name)
panel_summary <- panel |>
  dplyr::summarise(
    observations = dplyr::n(), stations = dplyr::n_distinct(code), months = dplyr::n_distinct(month),
    min_month = min(month), max_month = max(month), mean_no2 = mean(no2, na.rm = TRUE),
    treated_share = mean(d_active == 1)
  )
writeLines(c(
  "\\begin{table}[!htbp]\\centering", "\\caption{Exposure-support diagnostics for main ULEZ exposure mappings}",
  "\\label{tab:ulez_support_protocol}", "\\begin{threeparttable}", "\\begin{tabular}{lrrrr}",
  "\\toprule", "Exposure mapping & $D=0,S=0$ & $D=0,S>0$ & $D=1,S=0$ & $D=1,S>0$ \\\\", "\\midrule",
  paste0(support_main$matrix_label, " & ", support_main$d0_s0, " & ", support_main$d0_s1, " & ", support_main$d1_s0, " & ", support_main$d1_s1, " \\\\"),
  "\\midrule",
  paste0("\\multicolumn{5}{l}{Panel: ", panel_summary$observations, " station-months; ", panel_summary$stations, " stations; ", panel_summary$months, " months; ", panel_summary$min_month, "--", panel_summary$max_month, ".} \\\\"),
  paste0("\\multicolumn{5}{l}{Mean NO$_2$: ", fmt_num(panel_summary$mean_no2, 2), " $\\mu$g/m$^3$; treated observation share: ", fmt_num(100 * panel_summary$treated_share, 1), "\\%.} \\\\"),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: $D$ denotes direct treatment and $S>0$ denotes positive exposure under the corresponding matrix. The cell $D=0,S>0$ supports the currently untreated exposure parameter $\\beta^u$; $D=1,S=0$ helps separate $\\beta^d$ from the treated-exposure parameter $\\beta^t$.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_support_protocol.tex"))

# Residualized-design diagnostics.
resid_fe <- function(varname) {
  fml <- stats::as.formula(paste0(varname, " ~ 1 | code + month"))
  as.numeric(stats::residuals(fixest::feols(fml, data = panel, warn = FALSE, notes = FALSE)))
}
resid_design <- tibble::tibble(D = resid_fe("d_active"), DS = resid_fe("ds_bench"), US = resid_fe("us_bench"))
vif_one <- function(y, xs) {
  fml <- stats::as.formula(paste0(y, " ~ ", paste(xs, collapse = " + ")))
  1 / (1 - summary(stats::lm(fml, data = resid_design))$r.squared)
}
bench_support <- support |> dplyr::filter(matrix_name == "bench_0_3km") |> dplyr::slice(1)
design_diag <- tibble::tibble(
  diagnostic = c(
    "$D0/S0$", "$D0/S1$", "$D1/S0$", "$D1/S1$",
    "$\\mathrm{corr}(\\widetilde{D}_{it},\\widetilde{D_{it}S_{it}})$",
    "$\\mathrm{corr}(\\widetilde{D}_{it},\\widetilde{(1-D_{it})S_{it}})$",
    "$\\mathrm{corr}(\\widetilde{D_{it}S_{it}},\\widetilde{(1-D_{it})S_{it}})$",
    "$\\mathrm{VIF}(\\widetilde{D}_{it})$", "$\\mathrm{VIF}(\\widetilde{D_{it}S_{it}})$",
    "$\\mathrm{VIF}(\\widetilde{(1-D_{it})S_{it}})$"
  ),
  value = c(
    bench_support$d0_s0, bench_support$d0_s1, bench_support$d1_s0, bench_support$d1_s1,
    fmt_num(stats::cor(resid_design$D, resid_design$DS), 3),
    fmt_num(stats::cor(resid_design$D, resid_design$US), 3),
    fmt_num(stats::cor(resid_design$DS, resid_design$US), 3),
    fmt_num(vif_one("D", c("DS", "US")), 2),
    fmt_num(vif_one("DS", c("D", "US")), 2),
    fmt_num(vif_one("US", c("D", "DS")), 2)
  ),
  interpretation = c(
    "Untreated and unexposed", "Untreated and exposed", "Treated and unexposed", "Treated and exposed",
    "Residualized direct-treatment versus treated-exposure regressor",
    "Residualized direct-treatment versus currently untreated exposure regressor",
    "Residualized treated-exposure versus currently untreated exposure regressor",
    "Collinearity diagnostic for the direct-treatment regressor",
    "Collinearity diagnostic for the treated-exposure regressor",
    "Collinearity diagnostic for the currently untreated exposure regressor"
  )
)
writeLines(c(
  "\\begin{table}[!htbp]\\centering", "\\caption{Benchmark support and residualized-design diagnostics}",
  "\\label{tab:ulez_design_diagnostics}", "\\begin{threeparttable}", "\\small", "\\setlength{\\tabcolsep}{3pt}",
  "\\renewcommand{\\arraystretch}{1.10}",
  "\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}Xr>{\\RaggedRight\\arraybackslash}X@{}}",
  "\\toprule", "Diagnostic & Value & Interpretation \\\\", "\\midrule",
  paste0(design_diag$diagnostic, " & ", design_diag$value, " & ", design_diag$interpretation, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: Residualized variables are obtained after partialling out monitoring-station and month fixed effects. VIF denotes the variance-inflation factor computed from the residualized design matrix of the benchmark unrestricted model.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_design_diagnostics.tex"))

# Across-matrices table.
matrix_compact <- compare |>
  dplyr::mutate(matrix_name = factor(matrix_name, levels = matrix_specs$matrix_name)) |>
  dplyr::arrange(matrix_name) |>
  dplyr::transmute(
    matrix_label, d0_s1, d1_s0,
    beta_d = paste0(fmt_num(beta_d_effect_pct, 2), stars(beta_d_p_value)),
    beta_t = paste0(fmt_num(beta_t_effect_pct, 2), stars(beta_t_p_value)),
    beta_u = paste0(fmt_num(beta_u_effect_pct, 2), stars(beta_u_p_value)),
    beta_r = paste0(fmt_num(beta_r_effect_pct, 2), stars(beta_r_p_value)),
    p_equal = fmt_p(p_value_equal_bt_bu)
  )
writeLines(c(
  "\\begin{table}[!htbp]\\centering",
  "\\caption{Direct-treatment, treated-exposure, currently untreated exposure, and restricted exposure effects across matrices}",
  "\\label{tab:ulez_across_matrices_protocol}", "\\begin{threeparttable}", "\\scriptsize",
  "\\setlength{\\tabcolsep}{2pt}", "\\renewcommand{\\arraystretch}{1.08}",
  "\\begin{tabularx}{\\textwidth}{@{}>{\\RaggedRight\\arraybackslash}Xrrrrrrr@{}}",
  "\\toprule", "Mapping & $D0/S1$ & $D1/S0$ & $\\beta^d$ & $\\beta^t$ & $\\beta^u$ & $\\beta^r$ & $p_{tu}$ \\\\", "\\midrule",
  paste0(matrix_compact$matrix_label, " & ", matrix_compact$d0_s1, " & ", matrix_compact$d1_s0, " & ", matrix_compact$beta_d, " & ", matrix_compact$beta_t, " & ", matrix_compact$beta_u, " & ", matrix_compact$beta_r, " & ", matrix_compact$p_equal, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  "\\item Notes: Effects are expressed in percent. $\\beta^r$ is the common exposure coefficient from the restricted model. The equality test is computed in the unrestricted model; $p_{tu}$ denotes the p-value for $H_0:\\beta^t=\\beta^u$. Significance: * $p<0.10$, ** $p<0.05$, *** $p<0.01$.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_across_matrices_protocol.tex"))

# ------------------------------------------------------------
# Pre-trend gráfico: célula definida pelo primeiro episódio
# de tratamento ou exposição
# ------------------------------------------------------------

central_start <- as.Date("2019-04-01")
pre_start     <- as.Date("2016-01-01")
pre_end       <- as.Date("2019-03-01")
post_end      <- as.Date("2020-03-01")

# A célula fixa de cada estação é definida pelo primeiro mês
# em que tratamento ou exposição se tornam positivos.
# Estações nunca tratadas e nunca expostas permanecem D0/S0.
station_cells <- panel |>
  dplyr::arrange(code, month) |>
  dplyr::group_by(code) |>
  dplyr::summarise(
    cell = {
      idx <- which(d_active == 1L | s_pos_bench == 1L)
      
      if (!length(idx)) {
        "D0/S0"
      } else {
        paste0(
          "D", d_active[idx[1]],
          "/S", s_pos_bench[idx[1]]
        )
      }
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    cell = factor(
      cell,
      levels = c("D0/S0", "D0/S1", "D1/S0", "D1/S1")
    ),
    cell_label = dplyr::recode(
      as.character(cell),
      "D0/S0" = "D0/S0: untreated, unexposed",
      "D0/S1" = "D0/S1: untreated, exposed",
      "D1/S0" = "D1/S0: treated, unexposed",
      "D1/S1" = "D1/S1: treated, exposed"
    ),
    cell_label = factor(
      cell_label,
      levels = c(
        "D0/S0: untreated, unexposed",
        "D0/S1: untreated, exposed",
        "D1/S0: treated, unexposed",
        "D1/S1: treated, exposed"
      )
    )
  )

# Checagem de reprodução do suporte usado na figura.
expected_pretrend_support <- c(
  `D0/S0` = 27L,
  `D0/S1` = 18L,
  `D1/S0` = 17L,
  `D1/S1` = 67L
)

observed_pretrend_support <- station_cells |>
  dplyr::count(cell) |>
  dplyr::mutate(cell = as.character(cell)) |>
  tibble::deframe()

if (!identical(
  as.integer(observed_pretrend_support[names(expected_pretrend_support)]),
  as.integer(expected_pretrend_support)
)) {
  stop(
    "Unexpected pre-trend cell support. Observed: ",
    paste(
      names(observed_pretrend_support),
      observed_pretrend_support,
      collapse = "; "
    )
  )
}

pretrend_panel <- panel |>
  dplyr::left_join(
    station_cells |>
      dplyr::select(code, cell, cell_label),
    by = "code"
  ) |>
  dplyr::filter(
    !is.na(cell),
    month >= pre_start,
    month <= post_end
  ) |>
  dplyr::mutate(
    event_month =
      (as.integer(format(month, "%Y")) - 2019L) * 12L +
      (as.integer(format(month, "%m")) - 4L),
    
    month_index =
      (as.integer(format(month, "%Y")) - 2016L) * 12L +
      (as.integer(format(month, "%m")) - 1L)
  )

pretrend_mod <- fixest::feols(
  log_y ~ cell * month_index,
  data = pretrend_panel |>
    dplyr::filter(month <= pre_end),
  cluster = ~ code,
  warn = FALSE,
  notes = FALSE
)

event_means <- pretrend_panel |>
  dplyr::group_by(
    event_month,
    month,
    cell,
    cell_label
  ) |>
  dplyr::summarise(
    mean_log_no2 = mean(log_y, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::filter(
    event_month >= -36L,
    event_month <= 11L
  ) |>
  dplyr::arrange(event_month, cell)

# Figura principal.
p_pre <- ggplot2::ggplot(
  event_means,
  ggplot2::aes(
    x = event_month,
    y = mean_log_no2,
    linetype = cell_label,
    shape = cell_label,
    group = cell_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::labs(
    x = "Months relative to central ULEZ start",
    y = "Mean log(NO2)",
    linetype = "Future cell",
    shape = "Future cell"
  ) +
  ggplot2::guides(
    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

ggplot2::ggsave(
  file.path(
    FIGURE_DIR,
    "fig_ulez_benchmark_pretrend_cells_with_D1S0.png"
  ),
  p_pre,
  width = 14.17,
  height = 8.66,
  dpi = 300,
  bg = "white"
)


# Figura suplementar com painéis separados.
p_pre_facets <- ggplot2::ggplot(
  event_means,
  ggplot2::aes(
    x = event_month,
    y = mean_log_no2,
    group = cell_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.4) +
  ggplot2::facet_wrap(
    ~ cell_label,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::labs(
    x = "Months relative to central ULEZ start",
    y = "Mean log(NO2)"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  file.path(FIGURE_DIR, "Figure_S1.pdf"),
  p_pre_facets,
  width = 6.292,
  height = 4.722
)

# Restricted versus unrestricted benchmark figure.
plot_bench <- dplyr::bind_rows(
  coef_all |> dplyr::filter(model_name == "Restricted exposure", term == "s_bench") |>
    dplyr::mutate(parameter = "Restricted: beta_r"),
  coef_all |> dplyr::filter(model_name == "Unrestricted exposure", term %in% c("ds_bench", "us_bench")) |>
    dplyr::mutate(parameter = dplyr::if_else(term == "ds_bench", "Unrestricted: beta_t", "Unrestricted: beta_u"))
) |>
  dplyr::mutate(
    low = pct_effect(estimate - 1.96 * std_error), high = pct_effect(estimate + 1.96 * std_error),
    parameter = factor(parameter, levels = c("Restricted: beta_r", "Unrestricted: beta_t", "Unrestricted: beta_u"))
  )
p_bench <- ggplot2::ggplot(plot_bench, ggplot2::aes(effect_pct, parameter)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = low, xmax = high), height = 0.15) + ggplot2::geom_point() +
  ggplot2::labs(
    x = "Effect, percent", y = NULL,
    title = "Benchmark restricted and unrestricted exposure effects",
    subtitle = "Binary 0--3 km exposure mapping"
  ) + ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(FIGURE_DIR, "fig_ulez_benchmark_restricted_unrestricted.png"),
                p_bench, width = 7.5, height = 3.8, dpi = 300, bg = "white")

# Figure 3: label corrigido para a terminologia final do paper.
plot_matrix <- compare |>
  dplyr::mutate(matrix_label = factor(matrix_label, levels = matrix_specs$matrix_label)) |>
  dplyr::select(
    matrix_label,
    beta_t_effect_pct, beta_t_conf_low_pct, beta_t_conf_high_pct,
    beta_u_effect_pct, beta_u_conf_low_pct, beta_u_conf_high_pct,
    beta_r_effect_pct, beta_r_conf_low_pct, beta_r_conf_high_pct
  ) |>
  tidyr::pivot_longer(cols = -matrix_label, names_to = "name", values_to = "value") |>
  dplyr::mutate(
    parameter = dplyr::case_when(
      stringr::str_starts(name, "beta_t") ~ "beta_t: treated exposure",
      stringr::str_starts(name, "beta_u") ~ "beta_u: currently untreated exposure",
      stringr::str_starts(name, "beta_r") ~ "beta_r: restricted common exposure"
    ),
    stat = dplyr::case_when(
      stringr::str_detect(name, "conf_low") ~ "low",
      stringr::str_detect(name, "conf_high") ~ "high",
      TRUE ~ "estimate"
    )
  ) |>
  dplyr::select(matrix_label, parameter, stat, value) |>
  tidyr::pivot_wider(names_from = stat, values_from = value)
p_matrix <- ggplot2::ggplot(plot_matrix, ggplot2::aes(estimate, matrix_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = low, xmax = high), height = 0.15) + ggplot2::geom_point() +
  ggplot2::facet_wrap(~ parameter, scales = "free_x") +
  ggplot2::labs(
    x = "Effect, percent", y = NULL,
    title = "Exposure effects across admissible mappings",
    subtitle = "Alternative matrices are interpreted as alternative exposure hypotheses"
  ) + ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(FIGURE_DIR, "fig_ulez_across_matrices_protocol.png"),
                p_matrix, width = 11, height = 6, dpi = 300, bg = "white")

# Objetos internos para os apêndices; não são outputs do paper.
saveRDS(
  list(
    panel = panel, compare = compare, support = support,
    benchmark_models = list(no_exposure = m_no, restricted = m_res, unrestricted = m_unres),
    pretrend_model = pretrend_mod, event_means = event_means
  ),
  STATIC_RESULTS_RDS
)
message("[02] Resultados estáticos concluídos.")
