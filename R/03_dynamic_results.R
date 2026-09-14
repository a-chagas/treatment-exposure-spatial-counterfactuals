# ============================================================
# 03_dynamic_results.R
# Decomposição dinâmica / stacked usada no paper
# ============================================================
#
# A implementação preserva a regra de risk set da versão validada durante
# a revisão: unidades de coortes futuras permanecem no conjunto atualmente
# não tratado apenas enquanto t < C_i. A matriz de fontes é construída no
# universo completo do arquivo de assignment; receivers são as estações com
# outcome observado.
#
# Os diagnósticos auxiliares permanecem calculáveis dentro do core, mas apenas
# os objetos efetivamente usados no paper são exportados em outputs/.
# ============================================================

message("[03] Estimando decomposição dinâmica...")

PANEL_FILE_CANDIDATES <- c(RAW_PANEL_FILE)
ASSIGN_FILE_CANDIDATES <- c(ASSIGN_FILE)

OUT_DIR <- file.path(DERIVED_DIR, "dynamic_core")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# O paper reporta o benchmark de 3 km. Os demais raios pertenciam à trilha
# de auditoria da revisão e não são necessários ao patch público.
RADII_KM <- c(3)
EXPOSURE_EPS <- 1e-10
L_PRE <- 24
H_POST <- 15
EXPANSION_LABELS <- c(
  "Central ULEZ, Apr. 2019",
  "Inner London expansion, Oct. 2021",
  "London-wide expansion, Aug. 2023"
)
EVENT_REF_BIN <- "[-6,-1]"
EVENT_FIG_RADIUS <- 3

# 2. Helpers
# ------------------------------------------------------------

clean_names_simple <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop("Could not find ", label, ". Tried:\n", paste(paths, collapse = "\n"))
  }
  hit
}

guess_col <- function(dt, candidates, label) {
  hit <- intersect(candidates, names(dt))[1]
  if (is.na(hit) || length(hit) == 0) {
    stop(
      "Could not find column for ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(names(dt), collapse = ", ")
    )
  }
  hit
}

coerce_coord <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

first_nonmissing <- function(x) {
  y <- x[!is.na(x) & trimws(as.character(x)) != ""]
  if (length(y) == 0) return(NA)
  y[1]
}

first_nonmissing_chr <- function(x) {
  y <- as.character(x)
  y <- y[!is.na(y) & trimws(y) != ""]
  if (length(y) == 0) return(NA_character_)
  y[1]
}

first_finite_numeric <- function(x) {
  z <- coerce_coord(x)
  z <- z[is.finite(z)]
  if (length(z) == 0) return(NA_real_)
  z[1]
}

to01 <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "true", "t", "yes", "y", "sim")] <- 1L
  out[z %in% c("0", "false", "f", "no", "n", "nao", "não")] <- 0L
  out
}

max01 <- function(x) {
  z <- to01(x)
  if (all(is.na(z))) return(0L)
  as.integer(max(z, na.rm = TRUE))
}

as_month_date <- function(x) {
  if (inherits(x, "Date")) return(lubridate::floor_date(x, "month"))
  
  x_chr <- as.character(x)
  non_missing <- x_chr[!is.na(x_chr) & trimws(x_chr) != ""]
  
  if (length(non_missing) > 0 && all(grepl("^\\d{4}-\\d{2}$", non_missing))) {
    return(as.Date(paste0(x_chr, "-01")))
  }
  
  if (length(non_missing) > 0 && all(grepl("^\\d{6}$", non_missing))) {
    return(as.Date(paste0(substr(x_chr, 1, 4), "-", substr(x_chr, 5, 6), "-01")))
  }
  
  out <- suppressWarnings(as.Date(x_chr))
  if (all(!is.na(out[!is.na(x_chr) & trimws(x_chr) != ""]))) {
    return(lubridate::floor_date(out, "month"))
  }
  
  stop("Could not parse date/month column.")
}

month_id <- function(date) {
  lubridate::year(date) * 12L + lubridate::month(date)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

month_date_from_id <- function(x) {
  out <- rep(as.Date(NA), length(x))
  ok <- is.finite(x)
  if (any(ok)) {
    xx <- as.integer(x[ok])
    yy <- (xx - 1L) %/% 12L
    mm <- ((xx - 1L) %% 12L) + 1L
    out[ok] <- as.Date(sprintf("%04d-%02d-01", yy, mm))
  }
  out
}

assert_no_duplicate_keys <- function(dt, keys, label) {
  dup <- dt[, .N, by = keys][N > 1]
  if (nrow(dup) > 0) {
    print(head(dup, 20))
    stop("Duplicate keys detected in ", label, ": ", paste(keys, collapse = ", "))
  }
  invisible(TRUE)
}

assert_aligned_W <- function(W) {
  if (is.null(rownames(W)) || is.null(colnames(W))) {
    stop("W must have rownames and colnames.")
  }
  if (!identical(rownames(W), colnames(W))) {
    stop("W rows and columns must be the same station order.")
  }
  invisible(TRUE)
}

pct_effect <- function(beta) {
  100 * (exp(beta) - 1)
}

coef_table <- function(model, model_name, radius_km = NA_real_) {
  ct <- as.data.table(fixest::coeftable(model), keep.rownames = "term")
  setnames(
    ct,
    old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "z value", "Pr(>|z|)"),
    new = c("estimate", "std_error", "statistic", "p_value", "statistic", "p_value"),
    skip_absent = TRUE
  )
  ct[, model := model_name]
  ct[, radius_km := radius_km]
  ct[, effect_pct := pct_effect(estimate)]
  ct[, ci_low_pct := pct_effect(estimate - 1.96 * std_error)]
  ct[, ci_high_pct := pct_effect(estimate + 1.96 * std_error)]
  ct[]
}

build_W_cutoff <- function(cutoff_km, dist_mat, station_order) {
  cutoff_m <- cutoff_km * 1000
  W_bin <- (dist_mat > 0) & (dist_mat <= cutoff_m)
  
  W <- matrix(
    0,
    nrow = nrow(W_bin),
    ncol = ncol(W_bin),
    dimnames = list(station_order, station_order)
  )
  
  row_sums <- rowSums(W_bin)
  W[row_sums > 0, ] <- sweep(
    W_bin[row_sums > 0, , drop = FALSE] * 1,
    1,
    row_sums[row_sums > 0],
    "/"
  )
  
  assert_aligned_W(W)
  W
}

compute_source_exposure <- function(W, source_stations, exposure_eps = 1e-10,
                                    prefix = "source") {
  assert_aligned_W(W)
  
  z_source <- as.numeric(colnames(W) %in% source_stations)
  S_source <- as.numeric(W %*% z_source)
  
  out <- data.table(
    station_id = rownames(W),
    S_source = S_source,
    G_source = as.integer(S_source > exposure_eps)
  )
  
  setnames(
    out,
    old = c("S_source", "G_source"),
    new = c(paste0("S_", prefix), paste0("G_", prefix))
  )
  
  out[]
}

add_main_bins <- function(dt) {
  dt[, rel_bin_main := data.table::fcase(
    rel_k <= -13, "[-24,-13]",
    rel_k >= -12 & rel_k <= -7, "[-12,-7]",
    rel_k >= -6 & rel_k <= -1, "[-6,-1]",
    rel_k >= 0 & rel_k <= 5, "[0,5]",
    rel_k >= 6 & rel_k <= 15, "[6,15]"
  )]
  
  dt[, rel_bin_main := factor(
    rel_bin_main,
    levels = c("[-24,-13]", "[-12,-7]", "[-6,-1]", "[0,5]", "[6,15]")
  )]
  
  dt[]
}

rel_bin_main_mid <- data.table(
  rel_bin_main = c("[-24,-13]", "[-12,-7]", "[-6,-1]", "[0,5]", "[6,15]"),
  rel_bin_mid = c(-18.5, -9.5, -3.5, 2.5, 10.5)
)

# ------------------------------------------------------------
# 3. Load and harmonize data
# ------------------------------------------------------------

PANEL_FILE  <- first_existing(PANEL_FILE_CANDIDATES, "monthly ULEZ panel")
ASSIGN_FILE <- first_existing(ASSIGN_FILE_CANDIDATES, "station assignment file")

panel  <- fread(PANEL_FILE)
assign <- fread(ASSIGN_FILE)

setnames(panel, clean_names_simple(names(panel)))
setnames(assign, clean_names_simple(names(assign)))

# In the current ULEZ files, the identifier is usually called code.
panel_id_col <- guess_col(
  panel,
  c("station_id", "code", "site_code", "monitoring_station", "station"),
  "station id in panel"
)
assign_id_col <- guess_col(
  assign,
  c("station_id", "code", "site_code", "monitoring_station", "station"),
  "station id in assignment file"
)

panel[, station_id := toupper(trimws(as.character(get(panel_id_col))))]
assign[, station_id := toupper(trimws(as.character(get(assign_id_col))))]

# ------------------------------------------------------------
# 3.1 Clean assignment file: coordinates + ULEZ membership
# ------------------------------------------------------------

required_assign_vars <- c(
  "station_id",
  "latitude",
  "longitude",
  "inside_ulez_2019",
  "inside_ulez_2021",
  "inside_ulez_2023"
)

missing_assign_vars <- setdiff(required_assign_vars, names(assign))
if (length(missing_assign_vars) > 0) {
  stop(
    "assign is missing required variables: ",
    paste(missing_assign_vars, collapse = ", "),
    "\nAvailable columns:\n",
    paste(names(assign), collapse = ", ")
  )
}

assign_clean <- assign[
  ,
  .(
    source = if ("source" %in% names(assign)) first_nonmissing_chr(source) else NA_character_,
    site = if ("site" %in% names(assign)) first_nonmissing_chr(site) else NA_character_,
    address = if ("address" %in% names(assign)) first_nonmissing_chr(address) else NA_character_,
    la_id = if ("la_id" %in% names(assign)) first_nonmissing_chr(la_id) else NA_character_,
    authority = if ("authority" %in% names(assign)) first_nonmissing_chr(authority) else NA_character_,
    site_type = if ("site_type" %in% names(assign)) first_nonmissing_chr(site_type) else NA_character_,
    os_grid_x = if ("os_grid_x" %in% names(assign)) first_finite_numeric(os_grid_x) else NA_real_,
    os_grid_y = if ("os_grid_y" %in% names(assign)) first_finite_numeric(os_grid_y) else NA_real_,
    latitude = first_finite_numeric(latitude),
    longitude = first_finite_numeric(longitude),
    opening_date = if ("opening_date" %in% names(assign)) first_nonmissing_chr(opening_date) else NA_character_,
    closing_date = if ("closing_date" %in% names(assign)) first_nonmissing_chr(closing_date) else NA_character_,
    inside_ulez_2019 = max01(inside_ulez_2019),
    inside_ulez_2021 = max01(inside_ulez_2021),
    inside_ulez_2023 = max01(inside_ulez_2023)
  ),
  by = station_id
]

assert_no_duplicate_keys(assign_clean, "station_id", "assign_clean")

bad_assign_coords <- assign_clean[
  !is.finite(longitude) | !is.finite(latitude),
  .(station_id, longitude, latitude)
]

if (nrow(bad_assign_coords) > 0) {
  fwrite(
    bad_assign_coords,
    file.path(TAB_DIR, "table_12_00a_assign_missing_coordinates.csv")
  )
  stop(
    "assign_clean has stations with missing/non-finite coordinates. ",
    "See table_12_00a_assign_missing_coordinates.csv."
  )
}

# Define first-treatment cohort from nested ULEZ boundaries.
assign_clean[, cohort_date := as.Date(NA)]
assign_clean[inside_ulez_2019 == 1L, cohort_date := as.Date("2019-04-01")]
assign_clean[inside_ulez_2019 != 1L & inside_ulez_2021 == 1L, cohort_date := as.Date("2021-10-01")]
assign_clean[inside_ulez_2021 != 1L & inside_ulez_2023 == 1L, cohort_date := as.Date("2023-08-01")]

assign_clean[, cohort_id := month_id(cohort_date)]
assign_clean[is.na(cohort_id), cohort_id := Inf]

assign_clean[, cohort_label := data.table::fifelse(
  cohort_id == month_id(as.Date("2019-04-01")),
  "Central ULEZ 2019",
  data.table::fifelse(
    cohort_id == month_id(as.Date("2021-10-01")),
    "Inner London 2021",
    data.table::fifelse(
      cohort_id == month_id(as.Date("2023-08-01")),
      "London-wide 2023",
      "Never-treated"
    )
  )
)]

assignment_universe <- assign_clean[, .(
  n_stations = .N,
  n_2019 = sum(cohort_label == "Central ULEZ 2019"),
  n_2021 = sum(cohort_label == "Inner London 2021"),
  n_2023 = sum(cohort_label == "London-wide 2023"),
  n_never = sum(cohort_label == "Never-treated")
)]

fwrite(
  assignment_universe,
  file.path(TAB_DIR, "table_12_00_assignment_station_universe.csv")
)

# ------------------------------------------------------------
# 3.2 Clean outcome panel and merge assignment information
# ------------------------------------------------------------

date_col <- guess_col(
  panel,
  c("month", "date", "year_month", "month_date", "time"),
  "month/date"
)

panel[, month_date := as_month_date(get(date_col))]
panel[, m_id := month_id(month_date)]

if ("log_no2" %in% names(panel)) {
  panel[, y := as.numeric(log_no2)]
} else if ("ln_no2" %in% names(panel)) {
  panel[, y := as.numeric(ln_no2)]
} else if ("log_monthly_no2" %in% names(panel)) {
  panel[, y := as.numeric(log_monthly_no2)]
} else {
  no2_col <- guess_col(
    panel,
    c("no2", "mean_no2", "monthly_no2", "no2_mean", "avg_no2", "no2_monthly_avg"),
    "NO2 outcome"
  )
  panel[, y := log(as.numeric(get(no2_col)))]
}

panel <- panel[is.finite(y)]

panel <- merge(
  panel,
  assign_clean,
  by = "station_id",
  all.x = TRUE,
  sort = FALSE
)

unmatched_panel_stations <- panel[
  is.na(latitude) | is.na(longitude) | is.na(cohort_label),
  unique(station_id)
]

if (length(unmatched_panel_stations) > 0) {
  fwrite(
    data.table(station_id = unmatched_panel_stations),
    file.path(TAB_DIR, "table_12_00b_panel_stations_not_matched_to_assignment.csv")
  )
  
  stop(
    "Some panel stations were not matched to assign: ",
    paste(unmatched_panel_stations, collapse = ", "),
    ". See table_12_00b_panel_stations_not_matched_to_assignment.csv."
  )
}

# Direct treatment status in calendar time.
panel[, D := as.integer(is.finite(cohort_id) & m_id >= cohort_id)]

assert_no_duplicate_keys(panel, c("station_id", "m_id"), "panel")

# ------------------------------------------------------------
# 4. Station universe for W
# ------------------------------------------------------------
# W is built from the assignment universe, not from station-months observed
# in the outcome panel. This preserves the full set of potential exposure
# sources even when some stations are not observed in every stack window.

station_dt <- data.table::copy(assign_clean)
station_dt <- station_dt[is.finite(longitude) & is.finite(latitude)]
station_dt <- station_dt[order(station_id)]
assert_no_duplicate_keys(station_dt, "station_id", "station_dt")

station_order <- station_dt$station_id

station_sf <- sf::st_as_sf(
  as.data.frame(station_dt),
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) |>
  sf::st_transform(27700)

stopifnot(all(station_sf$station_id == station_order))

dist_mat_units <- sf::st_distance(station_sf)
dist_mat_numeric <- as.matrix(units::drop_units(dist_mat_units))
dimnames(dist_mat_numeric) <- list(station_order, station_order)

# ------------------------------------------------------------
# 5. Build dynamic stack with source-specific exposures
# ------------------------------------------------------------

make_stack_12 <- function(panel, W, station_dt, expansion_dates, expansion_labels,
                          L_pre = 24, H_post = 15, exposure_eps = 1e-10) {
  
  assert_aligned_W(W)
  expansion_ids <- month_id(expansion_dates)
  out_candidate <- vector("list", length(expansion_ids))
  out_analysis <- vector("list", length(expansion_ids))
  
  for (s in seq_along(expansion_ids)) {
    c_id <- expansion_ids[s]
    c_date <- expansion_dates[s]
    c_label <- expansion_labels[s]
    
    # Sources are defined on the full station universe used to build W.
    entering_stations <- station_dt[cohort_id == c_id, station_id]
    prior_stations <- station_dt[is.finite(cohort_id) & cohort_id < c_id, station_id]
    
    S_enter_dt <- compute_source_exposure(
      W = W,
      source_stations = entering_stations,
      exposure_eps = exposure_eps,
      prefix = "enter"
    )
    
    S_prior_dt <- compute_source_exposure(
      W = W,
      source_stations = prior_stations,
      exposure_eps = exposure_eps,
      prefix = "prior"
    )
    
    # Candidate stack window. The calendar-time risk-set rule below determines
    # which future-treated observations remain eligible for estimation.
    dt <- copy(panel[m_id >= c_id - L_pre & m_id <= c_id + H_post])
    
    dt[, stack_id := paste0("c_", format(c_date, "%Y_%m"))]
    dt[, stack_label := c_label]
    dt[, stack_cohort_id := c_id]
    dt[, stack_cohort_date := c_date]
    dt[, rel_k := m_id - c_id]
    dt[, post := as.integer(rel_k >= 0)]
    
    dt <- merge(dt, S_enter_dt, by = "station_id", all.x = TRUE)
    dt <- merge(dt, S_prior_dt, by = "station_id", all.x = TRUE)
    
    dt[is.na(S_enter), S_enter := 0]
    dt[is.na(G_enter), G_enter := 0]
    dt[is.na(S_prior), S_prior := 0]
    dt[is.na(G_prior), G_prior := 0]
    
    # Backward-compatible aliases.
    dt[, S_c := S_enter]
    dt[, G_c := G_enter]
    
    # ----------------------------------------------------------
    # Cohort identity relative to entering cohort c.
    # These indicators do NOT by themselves define current treatment status.
    # ----------------------------------------------------------
    dt[, new_cohort := as.integer(cohort_id == c_id)]
    dt[, future_cohort := as.integer(is.finite(cohort_id) & cohort_id > c_id)]
    dt[, prior_cohort := as.integer(is.finite(cohort_id) & cohort_id < c_id)]
    dt[, never_cohort := as.integer(!is.finite(cohort_id))]
    
    # ----------------------------------------------------------
    # Calendar-time risk set.
    # A future-cohort station is currently untreated only while m_id < C_i.
    # Once its own treatment begins, that observation is censored from stack c.
    # ----------------------------------------------------------
    dt[, future_at_risk := as.integer(
      future_cohort == 1 & is.finite(cohort_id) & m_id < cohort_id
    )]
    dt[, future_after_own_treatment := as.integer(
      future_cohort == 1 & is.finite(cohort_id) & m_id >= cohort_id
    )]
    dt[, risk_set_eligible := as.integer(future_after_own_treatment == 0)]
    
    dt[, currently_untreated_receiver := as.integer(
      future_at_risk == 1 | never_cohort == 1
    )]
    
    # Preserve the old cohort-fixed indicator only for the audit/comparison.
    dt[, legacy_untreated_receiver := as.integer(
      future_cohort == 1 | never_cohort == 1
    )]
    
    # Backward-compatible alias used by the remainder of the script.
    # In v4 this alias is explicitly calendar-time correct.
    dt[, untreated_receiver := currently_untreated_receiver]
    dt[, already_treated_receiver := as.integer(prior_cohort == 1)]
    dt[, treated_entrant := as.integer(new_cohort == 1)]
    
    dt[, receiver_status := data.table::fifelse(
      new_cohort == 1, "newly-treated cohort",
      data.table::fifelse(
        prior_cohort == 1, "already-treated",
        data.table::fifelse(
          never_cohort == 1, "never-treated",
          data.table::fifelse(
            future_at_risk == 1, "not-yet-treated",
            "future cohort after own treatment (censored)"
          )
        )
      )
    )]
    
    # Entrants split by prior exposure to previously treated cohorts.
    dt[, entrant_prior_exposed := treated_entrant * G_prior]
    dt[, entrant_prior_unexposed := treated_entrant * (1 - G_prior)]
    
    # Average post-treatment terms.
    dt[, post_new_unexposed_prior := post * entrant_prior_unexposed]
    dt[, post_new_exposed_prior := post * entrant_prior_exposed]
    
    # Incremental exposure to the entering cohort, by receiver status.
    dt[, post_new_S_enter := post * treated_entrant * S_enter]
    dt[, post_already_S_enter := post * already_treated_receiver * S_enter]
    dt[, post_untreated_S_enter := post * currently_untreated_receiver * S_enter]
    
    # Exact legacy version of the fifth component, used only to verify whether
    # the old cohort-fixed implementation affected the reported application.
    dt[, legacy_post_untreated_S_enter := post * legacy_untreated_receiver * S_enter]
    
    # Event-study variables.
    dt[, S_enter_new := treated_entrant * S_enter]
    dt[, S_enter_already := already_treated_receiver * S_enter]
    dt[, S_enter_untreated := currently_untreated_receiver * S_enter]
    dt[, legacy_S_enter_untreated := legacy_untreated_receiver * S_enter]
    
    dt[, station_stack := interaction(stack_id, station_id, drop = TRUE)]
    dt[, month_stack := interaction(stack_id, m_id, drop = TRUE)]
    
    # Keep the full candidate stack for the risk-set audit and the exact
    # legacy-vs-v4 comparison. Estimation uses only eligible observations.
    out_candidate[[s]] <- dt
    out_analysis[[s]] <- dt[risk_set_eligible == 1]
  }
  
  candidate_stack <- rbindlist(out_candidate, fill = TRUE)
  analysis_stack <- rbindlist(out_analysis, fill = TRUE)
  
  candidate_stack <- add_main_bins(candidate_stack)
  analysis_stack <- add_main_bins(analysis_stack)
  
  list(
    stack = analysis_stack[],
    candidate_stack = candidate_stack[]
  )
}

# Risk-set audit requested by Reviewer 2. The max/min pair is deliberately
# reported together with the unit-level violation count: max(m_id) < min(C_i)
# is a transparent sufficient diagnostic in the present application, while the
# violation count remains valid in more general staggered settings.
audit_future_treated_risk_set_12 <- function(candidate_stack, radius_km) {
  out <- candidate_stack[, {
    future_dt <- .SD[future_cohort == 1]
    future_exp_dt <- .SD[future_cohort == 1 & S_enter > EXPOSURE_EPS]
    
    max_m_future <- safe_max(future_dt$m_id)
    min_c_future <- safe_min(future_dt$cohort_id)
    max_m_future_exp <- safe_max(future_exp_dt$m_id)
    min_c_future_exp <- safe_min(future_exp_dt$cohort_id)
    
    list(
      radius_km = radius_km,
      stack_start_m_id = safe_min(m_id),
      stack_end_m_id = safe_max(m_id),
      n_future_stations = uniqueN(future_dt$station_id),
      n_future_station_months = nrow(future_dt),
      max_m_id_future = max_m_future,
      min_cohort_id_future = min_c_future,
      gap_months_future = if (is.finite(max_m_future) && is.finite(min_c_future)) {
        min_c_future - max_m_future
      } else NA_real_,
      n_future_obs_at_or_after_own_treatment = sum(
        future_dt$m_id >= future_dt$cohort_id,
        na.rm = TRUE
      ),
      n_future_exposed_stations = uniqueN(future_exp_dt$station_id),
      max_m_id_future_exposed = max_m_future_exp,
      min_cohort_id_future_exposed = min_c_future_exp,
      gap_months_future_exposed = if (
        is.finite(max_m_future_exp) && is.finite(min_c_future_exp)
      ) {
        min_c_future_exp - max_m_future_exp
      } else NA_real_,
      n_future_exposed_obs_at_or_after_own_treatment = sum(
        future_exp_dt$m_id >= future_exp_dt$cohort_id,
        na.rm = TRUE
      ),
      n_candidate_obs = .N,
      n_censored_obs = sum(risk_set_eligible == 0)
    )
  }, by = .(stack_id, stack_label, stack_cohort_id)]
  
  date_cols <- c(
    "stack_start_m_id", "stack_end_m_id",
    "max_m_id_future", "min_cohort_id_future",
    "max_m_id_future_exposed", "min_cohort_id_future_exposed"
  )
  for (v in date_cols) {
    # Avoid NSE/get() here: use the column vector explicitly so this audit
    # remains robust to evaluation environments and package masking.
    data.table::set(
      out,
      j = paste0(v, "_date"),
      value = month_date_from_id(out[[v]])
    )
  }
  
  out[]
}

validate_stack_sources <- function(stack_dt, W, station_dt, radius_km, exposure_eps = 1e-10) {
  assert_aligned_W(W)
  
  stack_keys <- unique(stack_dt[, .(stack_id, stack_label, stack_cohort_id)])
  
  out <- rbindlist(lapply(seq_len(nrow(stack_keys)), function(k) {
    sid <- stack_keys$stack_id[k]
    slab <- stack_keys$stack_label[k]
    c_id <- stack_keys$stack_cohort_id[k]
    
    entering_stations <- station_dt[cohort_id == c_id, station_id]
    prior_stations <- station_dt[is.finite(cohort_id) & cohort_id < c_id, station_id]
    
    S_enter_check <- compute_source_exposure(W, entering_stations, exposure_eps, "enter_check")
    S_prior_check <- compute_source_exposure(W, prior_stations, exposure_eps, "prior_check")
    
    tmp <- unique(stack_dt[stack_id == sid, .(
      station_id,
      S_enter,
      S_prior
    )])
    
    tmp <- merge(
      tmp,
      S_enter_check[, .(station_id, S_enter_check)],
      by = "station_id",
      all.x = TRUE
    )
    tmp <- merge(
      tmp,
      S_prior_check[, .(station_id, S_prior_check)],
      by = "station_id",
      all.x = TRUE
    )
    
    tmp[, diff_enter := S_enter - S_enter_check]
    tmp[, diff_prior := S_prior - S_prior_check]
    
    tmp[, .(
      radius_km = radius_km,
      stack_id = sid,
      stack_label = slab,
      n_observed_stations_in_stack = .N,
      entering_sources_full = length(entering_stations),
      prior_sources_full = length(prior_stations),
      mean_S_enter = safe_mean(S_enter),
      mean_S_prior = safe_mean(S_prior),
      max_abs_diff_enter = max(abs(diff_enter), na.rm = TRUE),
      max_abs_diff_prior = max(abs(diff_prior), na.rm = TRUE)
    )]
  }), fill = TRUE)
  
  if (any(out$max_abs_diff_enter > 1e-10) || any(out$max_abs_diff_prior > 1e-10)) {
    print(out)
    stop("Source exposure validation failed.")
  }
  
  out[]
}

# ------------------------------------------------------------
# 6. Model extraction helpers
# ------------------------------------------------------------

extract_implied_effects_12 <- function(model, stack_dt, radius_km) {
  
  ct <- as.data.table(fixest::coeftable(model), keep.rownames = "term")
  setnames(
    ct,
    old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "z value", "Pr(>|z|)"),
    new = c("estimate", "std_error", "statistic", "p_value", "statistic", "p_value"),
    skip_absent = TRUE
  )
  
  term_labels <- data.table(
    term = c(
      "post_new_unexposed_prior",
      "post_new_exposed_prior",
      "post_new_S_enter",
      "post_already_S_enter",
      "post_untreated_S_enter"
    ),
    component = c(
      "Direct effect on entrants not previously exposed",
      "Direct effect on entrants previously exposed",
      "Incremental exposure on entrants",
      "Incremental exposure on already-treated receivers",
      "Incremental exposure on currently untreated receivers"
    ),
    scale_type = c(
      "binary",
      "binary",
      "mean positive S_enter among entrants",
      "mean positive S_enter among already-treated receivers",
      "mean positive S_enter among untreated receivers"
    )
  )
  
  scale_dt <- data.table(
    term = term_labels$term,
    exposure_scale = c(
      1,
      1,
      safe_mean(stack_dt[treated_entrant == 1 & S_enter > EXPOSURE_EPS, S_enter]),
      safe_mean(stack_dt[already_treated_receiver == 1 & S_enter > EXPOSURE_EPS, S_enter]),
      safe_mean(stack_dt[untreated_receiver == 1 & S_enter > EXPOSURE_EPS, S_enter])
    )
  )
  
  out <- merge(term_labels, ct, by = "term", all.x = TRUE)
  out <- merge(out, scale_dt, by = "term", all.x = TRUE)
  
  out[, radius_km := radius_km]
  out[, implied_percent := 100 * (exp(estimate * exposure_scale) - 1)]
  out[, ci_low_percent := 100 * (exp((estimate - 1.96 * std_error) * exposure_scale) - 1)]
  out[, ci_high_percent := 100 * (exp((estimate + 1.96 * std_error) * exposure_scale) - 1)]
  
  setcolorder(
    out,
    c(
      "radius_km", "term", "component", "scale_type", "exposure_scale",
      "estimate", "std_error", "statistic", "p_value",
      "implied_percent", "ci_low_percent", "ci_high_percent"
    )
  )
  
  out[]
}

support_dynamic_decomp <- function(stack_dt, radius_km) {
  
  stack_dt[, .(
    radius_km = radius_km,
    n_obs = .N,
    n_stations = uniqueN(station_id),
    
    n_new_stations = uniqueN(station_id[treated_entrant == 1]),
    n_new_prior_unexposed_stations = uniqueN(station_id[treated_entrant == 1 & G_prior == 0]),
    n_new_prior_exposed_stations = uniqueN(station_id[treated_entrant == 1 & G_prior == 1]),
    
    n_new_S_enter_pos_stations = uniqueN(station_id[treated_entrant == 1 & S_enter > EXPOSURE_EPS]),
    n_already_S_enter_pos_stations = uniqueN(station_id[already_treated_receiver == 1 & S_enter > EXPOSURE_EPS]),
    n_untreated_S_enter_pos_stations = uniqueN(station_id[untreated_receiver == 1 & S_enter > EXPOSURE_EPS]),
    n_future_S_enter_pos_stations = uniqueN(station_id[future_cohort == 1 & S_enter > EXPOSURE_EPS]),
    n_never_S_enter_pos_stations = uniqueN(station_id[never_cohort == 1 & S_enter > EXPOSURE_EPS]),
    
    mean_S_prior_new = safe_mean(S_prior[treated_entrant == 1]),
    mean_S_prior_new_positive = safe_mean(S_prior[treated_entrant == 1 & S_prior > EXPOSURE_EPS]),
    
    mean_S_enter_new = safe_mean(S_enter[treated_entrant == 1]),
    mean_S_enter_new_positive = safe_mean(S_enter[treated_entrant == 1 & S_enter > EXPOSURE_EPS]),
    
    mean_S_enter_already = safe_mean(S_enter[already_treated_receiver == 1]),
    mean_S_enter_already_positive = safe_mean(S_enter[already_treated_receiver == 1 & S_enter > EXPOSURE_EPS]),
    
    mean_S_enter_untreated = safe_mean(S_enter[untreated_receiver == 1]),
    mean_S_enter_untreated_positive = safe_mean(S_enter[untreated_receiver == 1 & S_enter > EXPOSURE_EPS])
  ), by = .(stack_id, stack_label)][]
}

# ------------------------------------------------------------
# 6b. Cohort-specific coefficients and explicit aggregation
# ------------------------------------------------------------

component_spec_12 <- data.table(
  term = c(
    "post_new_unexposed_prior",
    "post_new_exposed_prior",
    "post_new_S_enter",
    "post_already_S_enter",
    "post_untreated_S_enter"
  ),
  component = c(
    "Direct effect on entrants not previously exposed",
    "Direct effect on entrants previously exposed",
    "Incremental exposure on entrants",
    "Incremental exposure on already-treated receivers",
    "Incremental exposure on currently untreated receivers"
  ),
  support_basis = c(
    "newly treated, not previously exposed stations",
    "newly treated, previously exposed stations",
    "newly treated stations with positive entering-cohort exposure",
    "already-treated stations with positive entering-cohort exposure",
    "currently untreated stations with positive entering-cohort exposure"
  )
)

cohort_component_support_12 <- function(stack_dt, radius_km) {
  keys <- unique(stack_dt[, .(stack_id, stack_label, stack_cohort_id)])
  
  rbindlist(lapply(seq_len(nrow(keys)), function(k) {
    sid <- keys$stack_id[k]
    slab <- keys$stack_label[k]
    cid <- keys$stack_cohort_id[k]
    dt <- stack_dt[stack_id == sid]
    
    data.table(
      radius_km = radius_km,
      stack_id = sid,
      stack_label = slab,
      stack_cohort_id = cid,
      term = component_spec_12$term,
      support_n = c(
        uniqueN(dt[treated_entrant == 1 & G_prior == 0, station_id]),
        uniqueN(dt[treated_entrant == 1 & G_prior == 1, station_id]),
        uniqueN(dt[treated_entrant == 1 & S_enter > EXPOSURE_EPS, station_id]),
        uniqueN(dt[already_treated_receiver == 1 & S_enter > EXPOSURE_EPS, station_id]),
        uniqueN(dt[untreated_receiver == 1 & S_enter > EXPOSURE_EPS, station_id])
      ),
      exposure_scale = c(
        1,
        1,
        safe_mean(dt[treated_entrant == 1 & S_enter > EXPOSURE_EPS, S_enter]),
        safe_mean(dt[already_treated_receiver == 1 & S_enter > EXPOSURE_EPS, S_enter]),
        safe_mean(dt[untreated_receiver == 1 & S_enter > EXPOSURE_EPS, S_enter])
      )
    )
  }), fill = TRUE)
}

fit_cohort_specific_joint_12 <- function(stack_dt, radius_km) {
  support_dt <- cohort_component_support_12(stack_dt, radius_km)
  mapping <- merge(
    support_dt,
    component_spec_12,
    by = "term",
    all.x = TRUE,
    sort = FALSE
  )
  
  work <- copy(stack_dt)
  mapping[, regressor := paste0(term, "__", stack_id)]
  
  for (j in seq_len(nrow(mapping))) {
    vv <- mapping$regressor[j]
    tt <- mapping$term[j]
    ss <- mapping$stack_id[j]
    work[, (vv) := get(tt) * as.integer(stack_id == ss)]
  }
  
  rhs <- paste(mapping$regressor, collapse = " + ")
  fml <- as.formula(paste0(
    "y ~ ", rhs,
    " | station_stack + month_stack"
  ))
  
  model <- fixest::feols(
    fml,
    cluster = ~ station_id,
    data = work
  )
  
  ct <- as.data.table(fixest::coeftable(model), keep.rownames = "regressor")
  setnames(
    ct,
    old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "z value", "Pr(>|z|)"),
    new = c("estimate", "std_error", "statistic", "p_value", "statistic", "p_value"),
    skip_absent = TRUE
  )
  
  estimates <- merge(mapping, ct, by = "regressor", all.x = TRUE, sort = FALSE)
  estimates[, estimable := as.integer(is.finite(estimate))]
  estimates[, implied_percent := 100 * (exp(estimate * exposure_scale) - 1)]
  estimates[, ci_low_percent := 100 * (
    exp((estimate - 1.96 * std_error) * exposure_scale) - 1
  )]
  estimates[, ci_high_percent := 100 * (
    exp((estimate + 1.96 * std_error) * exposure_scale) - 1
  )]
  
  list(
    model = model,
    estimates = estimates[],
    support = support_dt[]
  )
}

aggregate_cohort_specific_12 <- function(cohort_fit, radius_km) {
  model <- cohort_fit$model
  est_dt <- copy(cohort_fit$estimates)
  b <- stats::coef(model)
  V <- stats::vcov(model)
  
  aggregate_one <- function(term_name, weight_type) {
    tmp <- est_dt[
      term == term_name &
        estimable == 1 &
        regressor %in% names(b)
    ]
    
    if (nrow(tmp) == 0) return(NULL)
    
    if (weight_type == "simple") {
      w <- rep(1 / nrow(tmp), nrow(tmp))
    } else if (weight_type == "support_weighted") {
      if (sum(tmp$support_n, na.rm = TRUE) <= 0) return(NULL)
      w <- tmp$support_n / sum(tmp$support_n, na.rm = TRUE)
    } else {
      stop("Unknown weight_type: ", weight_type)
    }
    
    names(w) <- tmp$regressor
    b_sub <- b[tmp$regressor]
    V_sub <- V[tmp$regressor, tmp$regressor, drop = FALSE]
    estimate <- sum(w * b_sub)
    std_error <- sqrt(drop(t(w) %*% V_sub %*% w))
    
    data.table(
      radius_km = radius_km,
      term = term_name,
      component = tmp$component[1],
      aggregation = weight_type,
      n_cohorts_estimable = nrow(tmp),
      total_support_n = sum(tmp$support_n, na.rm = TRUE),
      estimate = estimate,
      std_error = std_error,
      statistic = estimate / std_error,
      p_value = 2 * stats::pnorm(-abs(estimate / std_error)),
      ci_low = estimate - 1.96 * std_error,
      ci_high = estimate + 1.96 * std_error,
      mean_exposure_scale = sum(w * tmp$exposure_scale, na.rm = TRUE),
      mean_cohort_implied_percent = sum(w * tmp$implied_percent, na.rm = TRUE),
      cohort_weights = paste0(
        tmp$stack_id, "=", sprintf("%.6f", w),
        collapse = ";"
      )
    )
  }
  
  out <- rbindlist(lapply(component_spec_12$term, function(tt) {
    rbindlist(list(
      aggregate_one(tt, "simple"),
      aggregate_one(tt, "support_weighted")
    ), fill = TRUE)
  }), fill = TRUE)
  
  out[]
}

compare_pooled_cohort_aggregation_12 <- function(pooled_coef, agg_dt, radius_km) {
  pool <- copy(pooled_coef)
  pool <- pool[term %in% component_spec_12$term, .(
    term,
    pooled_estimate = estimate,
    pooled_std_error = std_error
  )]
  
  out <- merge(
    agg_dt,
    pool,
    by = "term",
    all.x = TRUE,
    sort = FALSE
  )
  out[, difference_from_pooled := estimate - pooled_estimate]
  out[, radius_km := radius_km]
  out[]
}


# ------------------------------------------------------------
# 6c. Audits motivated by Reviewer 2 and by the cohort-specific decomposition
# ------------------------------------------------------------

fit_main_decomp_12 <- function(dt) {
  fixest::feols(
    y ~ post_new_unexposed_prior +
      post_new_exposed_prior +
      post_new_S_enter +
      post_already_S_enter +
      post_untreated_S_enter |
      station_stack + month_stack,
    cluster = ~ station_id,
    data = dt
  )
}

extract_one_coef_12 <- function(model, term_name) {
  if (is.null(model)) {
    return(data.table(estimate = NA_real_, std_error = NA_real_))
  }
  b <- stats::coef(model)
  V <- stats::vcov(model)
  if (!term_name %in% names(b)) {
    return(data.table(estimate = NA_real_, std_error = NA_real_))
  }
  data.table(
    estimate = unname(b[term_name]),
    std_error = sqrt(V[term_name, term_name])
  )
}

# Why this audit exists:
# Some cohort x component cells have structural zero support. For example, the
# 2019 ULEZ stack cannot contain "already treated" receivers or entrants that
# were exposed to an earlier ULEZ cohort. Those cohort-specific coefficients are
# therefore undefined (NA), not failed estimates.
#
# In the pooled common-slope model, a zero-support stack gives NO direct
# identifying variation for the target component. It can nevertheless affect
# that pooled coefficient indirectly because the other component slopes are
# constrained to be common across stacks. This audit separates those channels.
# For every component with at least one zero-support stack we compare:
#   (A) the reported full pooled common-slope estimate;
#   (B) the same pooled specification using only stacks with target support;
#   (C) a model using the full sample in which the TARGET slope is common across
#       supported stacks but all NUISANCE component slopes are stack-specific.
# If (A) differs from (B)/(C), the zero-support stack is not identifying the
# target coefficient directly; the difference comes from cross-component
# pooling restrictions.
audit_structural_zero_support_pooled_12 <- function(stack_dt, pooled_coef, radius_km) {
  support_dt <- cohort_component_support_12(stack_dt, radius_km)
  all_terms <- component_spec_12$term
  target_terms <- unique(support_dt[support_n == 0, term])

  if (length(target_terms) == 0) {
    return(data.table())
  }

  rbindlist(lapply(target_terms, function(tt) {
    supported_ids <- support_dt[term == tt & support_n > 0, stack_id]
    unsupported_ids <- support_dt[term == tt & support_n == 0, stack_id]

    full_row <- pooled_coef[term == tt]
    full_est <- if (nrow(full_row)) full_row$estimate[1] else NA_real_
    full_se  <- if (nrow(full_row)) full_row$std_error[1] else NA_real_

    # (B) Same common-slope decomposition, but only among stacks that have
    # positive support for the target component.
    mod_supported <- tryCatch(
      fit_main_decomp_12(stack_dt[stack_id %in% supported_ids]),
      error = function(e) NULL
    )
    supported_coef <- extract_one_coef_12(mod_supported, tt)

    # (C) Keep the full sample, but prevent unsupported stacks from influencing
    # the target coefficient through common nuisance slopes. The target term is
    # common; every other decomposition term is allowed to vary by stack.
    work <- copy(stack_dt)
    nuisance_terms <- setdiff(all_terms, tt)
    nuisance_regs <- character(0)

    for (nt in nuisance_terms) {
      for (sid in unique(work$stack_id)) {
        vv <- paste0(nt, "__NUIS__", sid)
        work[, (vv) := get(nt) * as.integer(stack_id == sid)]
        nuisance_regs <- c(nuisance_regs, vv)
      }
    }

    rhs <- paste(c(tt, nuisance_regs), collapse = " + ")
    fml <- as.formula(paste0(
      "y ~ ", rhs,
      " | station_stack + month_stack"
    ))

    mod_target_common <- tryCatch(
      fixest::feols(
        fml,
        cluster = ~ station_id,
        data = work
      ),
      error = function(e) NULL
    )
    target_common_coef <- extract_one_coef_12(mod_target_common, tt)

    data.table(
      radius_km = radius_km,
      term = tt,
      component = component_spec_12[term == tt, component][1],
      n_supported_stacks = length(supported_ids),
      n_unsupported_stacks = length(unsupported_ids),
      supported_stacks = paste(supported_ids, collapse = ";"),
      unsupported_stacks = paste(unsupported_ids, collapse = ";"),
      pooled_full_estimate = full_est,
      pooled_full_std_error = full_se,
      pooled_supported_stacks_estimate = supported_coef$estimate,
      pooled_supported_stacks_std_error = supported_coef$std_error,
      diff_supported_minus_full = supported_coef$estimate - full_est,
      target_common_nuisance_stack_specific_estimate = target_common_coef$estimate,
      target_common_nuisance_stack_specific_std_error = target_common_coef$std_error,
      diff_nuisance_heterogeneous_minus_full = target_common_coef$estimate - full_est
    )
  }), fill = TRUE)
}

# Why this audit exists:
# In the 3 km benchmark the 2023 cohort-specific coefficient for exposure on
# currently untreated receivers is estimated from only two positively exposed
# receivers. A conventional station-clustered standard error can look very small
# even when the identifying support for a particular treatment/exposure contrast
# is concentrated in very few receivers. The purpose of this audit is therefore
# NOT to search for a preferred estimate or to drop 2023 mechanically. It is to
# document how much identifying support exists, inspect the two supporting
# stations, and test sensitivity to each one individually. This is essential
# before interpreting cohort heterogeneity or comparing pooled and staggered
# estimates.
audit_2023_untreated_support_12 <- function(stack_dt, pooled_coef, radius_km,
                                             exposure_eps = EXPOSURE_EPS) {
  c2023 <- month_id(as.Date("2023-08-01"))
  dt23 <- copy(stack_dt[stack_cohort_id == c2023])

  support_stations <- unique(dt23[
    untreated_receiver == 1 & S_enter > exposure_eps,
    station_id
  ])

  station_audit <- dt23[station_id %in% support_stations,
    .(
      radius_km = radius_km,
      stack_label = unique(stack_label)[1],
      S_enter = unique(S_enter)[1],
      n_obs = .N,
      n_pre = sum(post == 0),
      n_post = sum(post == 1),
      mean_y_pre = safe_mean(y[post == 0]),
      mean_y_post = safe_mean(y[post == 1]),
      raw_post_minus_pre = safe_mean(y[post == 1]) - safe_mean(y[post == 0])
    ),
    by = station_id
  ]

  mod23 <- tryCatch(fit_main_decomp_12(dt23), error = function(e) NULL)
  baseline23 <- extract_one_coef_12(mod23, "post_untreated_S_enter")

  # Leave-one-support-station-out sensitivity. With only two support stations,
  # these estimates are intentionally diagnostic rather than inferential.
  loo <- if (length(support_stations) > 0) {
    rbindlist(lapply(support_stations, function(sid) {
      mod_loo <- tryCatch(
        fit_main_decomp_12(dt23[station_id != sid]),
        error = function(e) NULL
      )
      cc <- extract_one_coef_12(mod_loo, "post_untreated_S_enter")
      data.table(
        radius_km = radius_km,
        omitted_station_id = sid,
        n_support_stations_full = length(support_stations),
        n_support_stations_remaining = length(setdiff(support_stations, sid)),
        baseline_2023_estimate = baseline23$estimate,
        baseline_2023_std_error = baseline23$std_error,
        leave_one_out_estimate = cc$estimate,
        leave_one_out_std_error = cc$std_error,
        estimate_change = cc$estimate - baseline23$estimate
      )
    }), fill = TRUE)
  } else {
    data.table()
  }

  # Pooled sensitivity to the low-support 2023 stack. This does not choose
  # between pooled and cohort-specific estimands; it only shows whether the
  # paper's pooled currently-untreated result is materially driven by 2023.
  pooled_full_row <- pooled_coef[term == "post_untreated_S_enter"]
  pooled_pre23 <- tryCatch(
    fit_main_decomp_12(stack_dt[stack_cohort_id != c2023]),
    error = function(e) NULL
  )
  pooled_pre23_coef <- extract_one_coef_12(pooled_pre23, "post_untreated_S_enter")

  pooled_sensitivity <- data.table(
    radius_km = radius_km,
    n_2023_currently_untreated_stations = uniqueN(dt23[untreated_receiver == 1, station_id]),
    n_2023_positive_exposure_stations = length(support_stations),
    pooled_all_stacks_estimate = if (nrow(pooled_full_row)) pooled_full_row$estimate[1] else NA_real_,
    pooled_all_stacks_std_error = if (nrow(pooled_full_row)) pooled_full_row$std_error[1] else NA_real_,
    pooled_without_2023_estimate = pooled_pre23_coef$estimate,
    pooled_without_2023_std_error = pooled_pre23_coef$std_error,
    difference_without_2023_minus_all = pooled_pre23_coef$estimate -
      (if (nrow(pooled_full_row)) pooled_full_row$estimate[1] else NA_real_)
  )

  list(
    station_audit = station_audit[],
    leave_one_out = loo[],
    pooled_sensitivity = pooled_sensitivity[]
  )
}

# Why this diagnostic exists:
# Raw support counts (for example, the number of positively exposed receivers)
# are necessary but not sufficient to assess whether a cohort-specific slope is
# empirically well supported. After fixed effects and the other decomposition
# components are partialled out, the variation that identifies the target slope
# may be concentrated in only one or two station clusters. This is exactly what
# the 2023 / 3 km leave-one-out audit suggests for the currently-untreated
# exposure component.
#
# For each stack c and component a, use the Frisch-Waugh-Lovell residual of the
# target regressor x_it^{a,c} after partialling out station FE, month FE, and all
# other decomposition regressors within that stack. Let x_tilde be that residual.
# Define each station's share of identifying variation as
#
#   q_i^{a,c} = sum_t x_tilde_it^2 / sum_j sum_t x_tilde_jt^2,
#
# and the effective identifying support as
#
#   N_eff^{a,c} = 1 / sum_i (q_i^{a,c})^2.
#
# N_eff is a concentration diagnostic, NOT an automatic sample-selection rule
# and NOT a replacement for the raw support counts. A low N_eff means that the
# cohort-specific coefficient relies on a small effective number of station
# clusters after residualization. This helps compare the support cost of allowing
# cohort heterogeneity relative to the pooled common-slope specification.
effective_identifying_support_12 <- function(stack_dt, radius_km,
                                              residual_tol = 1e-12) {
  support_dt <- cohort_component_support_12(stack_dt, radius_km)
  all_terms <- component_spec_12$term
  stack_keys <- unique(stack_dt[, .(stack_id, stack_label, stack_cohort_id)])

  # IMPORTANT INTERPRETATION NOTE:
  # n_stations_with_nonzero_residualized_target is NOT a count of raw-support
  # stations. After Frisch-Waugh-Lovell residualization, stations whose raw target
  # regressor is identically zero can receive nonzero x_tilde because fixed effects
  # and nuisance regressors are partialled out. Those stations are comparison
  # observations in the partial regression. Whether they materially influence the
  # coefficient is assessed below from (i) their share of sum(x_tilde^2) and
  # (ii) their contribution to sum(x_tilde*y_tilde)/sum(x_tilde^2), not from the
  # residualized-station count itself.
  support_note <- paste(
    "n_stations_with_nonzero_residualized_target is not raw support.",
    "FWL residualization can assign nonzero target residuals to comparison stations",
    "whose raw target regressor is zero. Interpret this count together with raw_support_n,",
    "raw_nonzero_regressor_stations, n_eff, q shares, and FWL coefficient contributions."
  )

  summary_list <- list()
  station_list <- list()
  kk <- 0L

  for (ss in seq_len(nrow(stack_keys))) {
    sid <- stack_keys$stack_id[ss]
    slab <- stack_keys$stack_label[ss]
    cid <- stack_keys$stack_cohort_id[ss]
    dt <- copy(stack_dt[stack_id == sid])

    for (tt in all_terms) {
      kk <- kk + 1L
      comp <- component_spec_12[term == tt, component][1]
      support_n <- support_dt[stack_id == sid & term == tt, support_n][1]
      if (length(support_n) == 0 || is.na(support_n)) support_n <- 0L

      dt[, .target_x_support := get(tt)]
      raw_nonzero_stations <- dt[
        abs(.target_x_support) > residual_tol,
        unique(station_id)
      ]

      nuisance_terms <- setdiff(all_terms, tt)
      rhs <- paste(nuisance_terms, collapse = " + ")
      fml_x <- as.formula(paste0(
        ".target_x_support ~ ", rhs,
        " | station_id + m_id"
      ))
      fml_y <- as.formula(paste0(
        "y ~ ", rhs,
        " | station_id + m_id"
      ))

      mod_x <- tryCatch(
        fixest::feols(fml_x, data = dt, fixef.rm = "none"),
        error = function(e) NULL
      )
      mod_y <- tryCatch(
        fixest::feols(fml_y, data = dt, fixef.rm = "none"),
        error = function(e) NULL
      )

      status <- "ok"
      if (is.null(mod_x) || is.null(mod_y)) status <- "auxiliary regression failed"

      if (status != "ok") {
        summary_list[[kk]] <- data.table(
          radius_km = radius_km,
          stack_id = sid,
          stack_label = slab,
          stack_cohort_id = cid,
          term = tt,
          component = comp,
          raw_support_n = as.integer(support_n),
          raw_nonzero_regressor_stations = length(raw_nonzero_stations),
          residual_ss_total = NA_real_,
          n_stations_with_nonzero_residualized_target = 0L,
          n_eff = NA_real_,
          max_share = NA_real_,
          top2_share = NA_real_,
          top3_share = NA_real_,
          max_share_station_id = NA_character_,
          residual_ss_raw_support_station_share = NA_real_,
          residual_ss_raw_zero_station_share = NA_real_,
          beta_fwl = NA_real_,
          beta_contribution_raw_support_stations = NA_real_,
          beta_contribution_raw_zero_stations = NA_real_,
          residualization_status = status,
          support_interpretation_note = support_note
        )
        next
      }

      x_tilde <- as.numeric(stats::resid(mod_x))
      y_tilde <- as.numeric(stats::resid(mod_y))
      if (length(x_tilde) != nrow(dt) || length(y_tilde) != nrow(dt)) {
        summary_list[[kk]] <- data.table(
          radius_km = radius_km,
          stack_id = sid,
          stack_label = slab,
          stack_cohort_id = cid,
          term = tt,
          component = comp,
          raw_support_n = as.integer(support_n),
          raw_nonzero_regressor_stations = length(raw_nonzero_stations),
          residual_ss_total = NA_real_,
          n_stations_with_nonzero_residualized_target = 0L,
          n_eff = NA_real_,
          max_share = NA_real_,
          top2_share = NA_real_,
          top3_share = NA_real_,
          max_share_station_id = NA_character_,
          residual_ss_raw_support_station_share = NA_real_,
          residual_ss_raw_zero_station_share = NA_real_,
          beta_fwl = NA_real_,
          beta_contribution_raw_support_stations = NA_real_,
          beta_contribution_raw_zero_stations = NA_real_,
          residualization_status = "residual length mismatch",
          support_interpretation_note = support_note
        )
        next
      }

      dt[, `:=`(
        .x_tilde_support = x_tilde,
        .y_tilde_support = y_tilde
      )]

      station_q <- dt[, .(
        residual_ss = sum(.x_tilde_support^2, na.rm = TRUE),
        residual_xy = sum(.x_tilde_support * .y_tilde_support, na.rm = TRUE),
        raw_target_ever_nonzero = as.integer(any(abs(.target_x_support) > residual_tol))
      ), by = station_id]

      total_ss <- sum(station_q$residual_ss, na.rm = TRUE)
      total_xy <- sum(station_q$residual_xy, na.rm = TRUE)

      if (!is.finite(total_ss) || total_ss <= residual_tol) {
        summary_list[[kk]] <- data.table(
          radius_km = radius_km,
          stack_id = sid,
          stack_label = slab,
          stack_cohort_id = cid,
          term = tt,
          component = comp,
          raw_support_n = as.integer(support_n),
          raw_nonzero_regressor_stations = length(raw_nonzero_stations),
          residual_ss_total = total_ss,
          n_stations_with_nonzero_residualized_target = 0L,
          n_eff = NA_real_,
          max_share = NA_real_,
          top2_share = NA_real_,
          top3_share = NA_real_,
          max_share_station_id = NA_character_,
          residual_ss_raw_support_station_share = NA_real_,
          residual_ss_raw_zero_station_share = NA_real_,
          beta_fwl = NA_real_,
          beta_contribution_raw_support_stations = NA_real_,
          beta_contribution_raw_zero_stations = NA_real_,
          residualization_status = "no residual identifying variation",
          support_interpretation_note = support_note
        )
        next
      }

      station_q[, `:=`(
        q_share = residual_ss / total_ss,
        beta_contribution = residual_xy / total_ss
      )]
      station_q <- station_q[is.finite(q_share) & q_share > residual_tol]
      setorder(station_q, -q_share, station_id)
      station_q[, q_rank := seq_len(.N)]

      q <- station_q$q_share
      n_eff <- 1 / sum(q^2)
      top_share <- function(k) sum(head(q, k))

      ss_raw_support <- sum(station_q[raw_target_ever_nonzero == 1, residual_ss], na.rm = TRUE)
      ss_raw_zero <- sum(station_q[raw_target_ever_nonzero == 0, residual_ss], na.rm = TRUE)
      beta_raw_support <- sum(station_q[raw_target_ever_nonzero == 1, residual_xy], na.rm = TRUE) / total_ss
      beta_raw_zero <- sum(station_q[raw_target_ever_nonzero == 0, residual_xy], na.rm = TRUE) / total_ss

      summary_list[[kk]] <- data.table(
        radius_km = radius_km,
        stack_id = sid,
        stack_label = slab,
        stack_cohort_id = cid,
        term = tt,
        component = comp,
        raw_support_n = as.integer(support_n),
        raw_nonzero_regressor_stations = length(raw_nonzero_stations),
        residual_ss_total = total_ss,
        n_stations_with_nonzero_residualized_target = nrow(station_q),
        n_eff = n_eff,
        max_share = if (length(q)) q[1] else NA_real_,
        top2_share = if (length(q)) top_share(2) else NA_real_,
        top3_share = if (length(q)) top_share(3) else NA_real_,
        max_share_station_id = if (nrow(station_q)) station_q$station_id[1] else NA_character_,
        residual_ss_raw_support_station_share = ss_raw_support / total_ss,
        residual_ss_raw_zero_station_share = ss_raw_zero / total_ss,
        beta_fwl = total_xy / total_ss,
        beta_contribution_raw_support_stations = beta_raw_support,
        beta_contribution_raw_zero_stations = beta_raw_zero,
        residualization_status = "ok",
        support_interpretation_note = support_note
      )

      station_q[, `:=`(
        radius_km = radius_km,
        stack_id = sid,
        stack_label = slab,
        stack_cohort_id = cid,
        term = tt,
        component = comp,
        raw_support_n = as.integer(support_n)
      )]
      setcolorder(
        station_q,
        c(
          "radius_km", "stack_id", "stack_label", "stack_cohort_id",
          "term", "component", "raw_support_n", "station_id",
          "raw_target_ever_nonzero", "residual_ss", "residual_xy",
          "q_share", "beta_contribution", "q_rank"
        )
      )
      station_list[[kk]] <- station_q[]
    }
  }

  list(
    summary = rbindlist(summary_list, fill = TRUE),
    station_shares = rbindlist(station_list, fill = TRUE)
  )
}

# Pooled analogue of the effective-support audit.
# Here the target coefficient is the reported common-slope stacked coefficient.
# The target regressor is residualized on the four other COMMON decomposition
# terms plus station-stack and month-stack fixed effects, exactly matching the
# FWL representation of the pooled specification. We decompose identifying
# variation and the FWL numerator by station and by expansion stack.
pooled_effective_identifying_support_12 <- function(stack_dt, pooled_coef, radius_km,
                                                     residual_tol = 1e-12) {
  all_terms <- component_spec_12$term
  support_dt <- cohort_component_support_12(stack_dt, radius_km)
  support_note <- paste(
    "n_stations_with_nonzero_residualized_target is not raw support.",
    "FWL residualization can assign nonzero target residuals to comparison stations",
    "whose raw target regressor is zero. Whether this expanded residualized set drives",
    "the pooled coefficient is assessed by residual-SS shares and beta contributions."
  )

  summary_list <- list()
  station_list <- list()
  stack_list <- list()

  for (jj in seq_along(all_terms)) {
    tt <- all_terms[jj]
    comp <- component_spec_12[term == tt, component][1]
    dt <- copy(stack_dt)
    dt[, .target_x_pooled := get(tt)]

    nuisance_terms <- setdiff(all_terms, tt)
    rhs <- paste(nuisance_terms, collapse = " + ")
    fml_x <- as.formula(paste0(
      ".target_x_pooled ~ ", rhs,
      " | station_stack + month_stack"
    ))
    fml_y <- as.formula(paste0(
      "y ~ ", rhs,
      " | station_stack + month_stack"
    ))

    mod_x <- tryCatch(fixest::feols(fml_x, data = dt, fixef.rm = "none"), error = function(e) NULL)
    mod_y <- tryCatch(fixest::feols(fml_y, data = dt, fixef.rm = "none"), error = function(e) NULL)

    pooled_row <- pooled_coef[term == tt]
    pooled_est <- if (nrow(pooled_row)) pooled_row$estimate[1] else NA_real_

    if (is.null(mod_x) || is.null(mod_y)) {
      summary_list[[jj]] <- data.table(
        radius_km = radius_km, term = tt, component = comp,
        pooled_estimate = pooled_est, beta_fwl = NA_real_, fwl_minus_pooled = NA_real_,
        raw_nonzero_regressor_stations = uniqueN(dt[abs(.target_x_pooled) > residual_tol, station_id]),
        n_stations_with_nonzero_residualized_target = 0L,
        n_eff_station = NA_real_, max_station_share = NA_real_, top2_station_share = NA_real_,
        residual_ss_raw_support_station_share = NA_real_, residual_ss_raw_zero_station_share = NA_real_,
        beta_contribution_raw_support_stations = NA_real_, beta_contribution_raw_zero_stations = NA_real_,
        residualization_status = "auxiliary regression failed",
        support_interpretation_note = support_note
      )
      next
    }

    x_tilde <- as.numeric(stats::resid(mod_x))
    y_tilde <- as.numeric(stats::resid(mod_y))
    if (length(x_tilde) != nrow(dt) || length(y_tilde) != nrow(dt)) {
      stop("Pooled effective-support residual length mismatch for term: ", tt)
    }

    dt[, `:=`(
      .x_tilde_pooled = x_tilde,
      .y_tilde_pooled = y_tilde
    )]

    station_q <- dt[, .(
      residual_ss = sum(.x_tilde_pooled^2, na.rm = TRUE),
      residual_xy = sum(.x_tilde_pooled * .y_tilde_pooled, na.rm = TRUE),
      raw_target_ever_nonzero = as.integer(any(abs(.target_x_pooled) > residual_tol))
    ), by = station_id]

    total_ss <- sum(station_q$residual_ss, na.rm = TRUE)
    total_xy <- sum(station_q$residual_xy, na.rm = TRUE)
    if (!is.finite(total_ss) || total_ss <= residual_tol) next

    station_q[, `:=`(
      q_share = residual_ss / total_ss,
      beta_contribution = residual_xy / total_ss
    )]
    station_q <- station_q[is.finite(q_share) & q_share > residual_tol]
    setorder(station_q, -q_share, station_id)
    station_q[, q_rank := seq_len(.N)]

    stack_q <- dt[, .(
      residual_ss = sum(.x_tilde_pooled^2, na.rm = TRUE),
      residual_xy = sum(.x_tilde_pooled * .y_tilde_pooled, na.rm = TRUE),
      raw_nonzero_regressor_stations = uniqueN(station_id[abs(.target_x_pooled) > residual_tol])
    ), by = .(stack_id, stack_label, stack_cohort_id)]
    stack_q <- merge(
      stack_q,
      support_dt[term == tt, .(stack_id, raw_support_n = support_n)],
      by = "stack_id", all.x = TRUE, sort = FALSE
    )
    stack_q[, `:=`(
      q_share = residual_ss / total_ss,
      beta_contribution = residual_xy / total_ss,
      radius_km = radius_km,
      term = tt,
      component = comp
    )]
    setorder(stack_q, -q_share, stack_id)
    stack_q[, q_rank := seq_len(.N)]

    q <- station_q$q_share
    n_eff <- 1 / sum(q^2)
    ss_raw_support <- sum(station_q[raw_target_ever_nonzero == 1, residual_ss], na.rm = TRUE)
    ss_raw_zero <- sum(station_q[raw_target_ever_nonzero == 0, residual_ss], na.rm = TRUE)
    beta_raw_support <- sum(station_q[raw_target_ever_nonzero == 1, residual_xy], na.rm = TRUE) / total_ss
    beta_raw_zero <- sum(station_q[raw_target_ever_nonzero == 0, residual_xy], na.rm = TRUE) / total_ss
    beta_fwl <- total_xy / total_ss

    summary_list[[jj]] <- data.table(
      radius_km = radius_km,
      term = tt,
      component = comp,
      pooled_estimate = pooled_est,
      beta_fwl = beta_fwl,
      fwl_minus_pooled = beta_fwl - pooled_est,
      raw_nonzero_regressor_stations = uniqueN(dt[abs(.target_x_pooled) > residual_tol, station_id]),
      n_stations_with_nonzero_residualized_target = nrow(station_q),
      n_eff_station = n_eff,
      max_station_share = if (length(q)) q[1] else NA_real_,
      top2_station_share = if (length(q)) sum(head(q, 2)) else NA_real_,
      residual_ss_raw_support_station_share = ss_raw_support / total_ss,
      residual_ss_raw_zero_station_share = ss_raw_zero / total_ss,
      beta_contribution_raw_support_stations = beta_raw_support,
      beta_contribution_raw_zero_stations = beta_raw_zero,
      residualization_status = "ok",
      support_interpretation_note = support_note
    )

    station_q[, `:=`(
      radius_km = radius_km,
      term = tt,
      component = comp
    )]
    setcolorder(
      station_q,
      c("radius_km", "term", "component", "station_id", "raw_target_ever_nonzero",
        "residual_ss", "residual_xy", "q_share", "beta_contribution", "q_rank")
    )
    station_list[[jj]] <- station_q[]

    setcolorder(
      stack_q,
      c("radius_km", "term", "component", "stack_id", "stack_label", "stack_cohort_id",
        "raw_support_n", "raw_nonzero_regressor_stations", "residual_ss", "residual_xy",
        "q_share", "beta_contribution", "q_rank")
    )
    stack_list[[jj]] <- stack_q[]
  }

  list(
    summary = rbindlist(summary_list, fill = TRUE),
    station_shares = rbindlist(station_list, fill = TRUE),
    stack_shares = rbindlist(stack_list, fill = TRUE)
  )
}

compare_legacy_riskset_12 <- function(candidate_stack, riskset_stack, radius_km) {
  legacy_model <- fixest::feols(
    y ~ post_new_unexposed_prior +
      post_new_exposed_prior +
      post_new_S_enter +
      post_already_S_enter +
      legacy_post_untreated_S_enter |
      station_stack + month_stack,
    cluster = ~ station_id,
    data = candidate_stack
  )
  
  risk_model <- fixest::feols(
    y ~ post_new_unexposed_prior +
      post_new_exposed_prior +
      post_new_S_enter +
      post_already_S_enter +
      post_untreated_S_enter |
      station_stack + month_stack,
    cluster = ~ station_id,
    data = riskset_stack
  )
  
  legacy <- coef_table(
    legacy_model,
    model_name = "legacy_cohort_fixed_untreated",
    radius_km = radius_km
  )
  legacy[term == "legacy_post_untreated_S_enter", term := "post_untreated_S_enter"]
  legacy <- legacy[, .(
    term,
    legacy_estimate = estimate,
    legacy_std_error = std_error
  )]
  
  current <- coef_table(
    risk_model,
    model_name = "calendar_time_risk_set",
    radius_km = radius_km
  )[, .(
    term,
    riskset_estimate = estimate,
    riskset_std_error = std_error
  )]
  
  out <- merge(legacy, current, by = "term", all = TRUE)
  out[, radius_km := radius_km]
  out[, estimate_difference := riskset_estimate - legacy_estimate]
  out[, std_error_difference := riskset_std_error - legacy_std_error]
  out[, n_candidate_obs := nobs(legacy_model)]
  out[, n_riskset_obs := nobs(risk_model)]
  out[, n_obs_removed := n_candidate_obs - n_riskset_obs]
  
  list(
    table = out[],
    legacy_model = legacy_model,
    risk_model = risk_model
  )
}

# ------------------------------------------------------------
# 7. Radius-level routine
# ------------------------------------------------------------

run_pipeline_12_radius <- function(radius_km) {
  
  message("Running pipeline 12 v4 for radius: ", radius_km, " km")
  
  W_r <- build_W_cutoff(
    cutoff_km = radius_km,
    dist_mat = dist_mat_numeric,
    station_order = station_order
  )
  
  stack_bundle_r <- make_stack_12(
    panel = panel,
    W = W_r,
    station_dt = station_dt,
    expansion_dates = EXPANSION_DATES,
    expansion_labels = EXPANSION_LABELS,
    L_pre = L_PRE,
    H_post = H_POST,
    exposure_eps = EXPOSURE_EPS
  )
  
  stack_r <- stack_bundle_r$stack
  candidate_stack_r <- stack_bundle_r$candidate_stack
  
  riskset_audit_r <- audit_future_treated_risk_set_12(
    candidate_stack = candidate_stack_r,
    radius_km = radius_km
  )
  
  validation_r <- validate_stack_sources(
    stack_dt = stack_r,
    W = W_r,
    station_dt = station_dt,
    radius_km = radius_km,
    exposure_eps = EXPOSURE_EPS
  )
  
  support_r <- support_dynamic_decomp(stack_r, radius_km)
  
  # Main pooled decomposition: common coefficients across stacks.
  mod_decomp_r <- fixest::feols(
    y ~ post_new_unexposed_prior +
      post_new_exposed_prior +
      post_new_S_enter +
      post_already_S_enter +
      post_untreated_S_enter |
      station_stack + month_stack,
    cluster = ~ station_id,
    data = stack_r
  )
  
  coef_r <- coef_table(
    model = mod_decomp_r,
    model_name = "dynamic_static_aligned_decomposition_v4_riskset",
    radius_km = radius_km
  )
  
  implied_r <- extract_implied_effects_12(
    model = mod_decomp_r,
    stack_dt = stack_r,
    radius_km = radius_km
  )
  
  # Exact comparison with the previous cohort-fixed untreated definition.
  legacy_compare_r <- compare_legacy_riskset_12(
    candidate_stack = candidate_stack_r,
    riskset_stack = stack_r,
    radius_km = radius_km
  )
  
  # If no candidate observation violates the calendar-time risk set, the v4
  # point estimates should reproduce the legacy implementation numerically.
  if (sum(riskset_audit_r$n_censored_obs, na.rm = TRUE) == 0) {
    max_diff <- max(abs(legacy_compare_r$table$estimate_difference), na.rm = TRUE)
    if (is.finite(max_diff) && max_diff > 1e-10) {
      stop(
        "Risk-set audit found zero censored observations, but legacy and v4 ",
        "coefficients differ by more than tolerance. Max difference: ", max_diff
      )
    }
  }
  
  # Cohort-specific coefficients. The point estimates are stack-specific, while
  # the joint model retains cross-stack covariance under station clustering.
  cohort_fit_r <- fit_cohort_specific_joint_12(
    stack_dt = stack_r,
    radius_km = radius_km
  )
  
  cohort_agg_r <- aggregate_cohort_specific_12(
    cohort_fit = cohort_fit_r,
    radius_km = radius_km
  )
  
  pooled_vs_agg_r <- compare_pooled_cohort_aggregation_12(
    pooled_coef = coef_r,
    agg_dt = cohort_agg_r,
    radius_km = radius_km
  )
  
  # Structural-zero-support audit: verifies whether cohort/component cells that
  # are undefined by construction can influence the pooled coefficient only
  # through the model's common-slope restrictions on other components.
  zero_support_pooled_audit_r <- audit_structural_zero_support_pooled_12(
    stack_dt = stack_r,
    pooled_coef = coef_r,
    radius_km = radius_km
  )
  
  # 2023 support audit: motivated by the unusually small clustered SE for the
  # currently-untreated exposure coefficient despite very limited positive
  # exposure support. This is a support/sensitivity diagnostic, not a rule for
  # excluding the 2023 cohort.
  audit_2023_r <- audit_2023_untreated_support_12(
    stack_dt = stack_r,
    pooled_coef = coef_r,
    radius_km = radius_km,
    exposure_eps = EXPOSURE_EPS
  )
  
  # Effective identifying support: after partialling out fixed effects and all
  # other decomposition components, measure how concentrated the target
  # regressor's identifying variation is across station clusters. This is a
  # support diagnostic for cohort-specific heterogeneity, not an exclusion rule.
  effective_support_r <- effective_identifying_support_12(
    stack_dt = stack_r,
    radius_km = radius_km
  )
  
  # Pooled FWL support audit: residualize the reported common-slope target on
  # the other common decomposition terms and stack-specific fixed effects. This
  # shows whether the larger residualized station set contributes meaningful
  # identifying variance or actually pushes the pooled coefficient.
  pooled_effective_support_r <- pooled_effective_identifying_support_12(
    stack_dt = stack_r,
    pooled_coef = coef_r,
    radius_km = radius_km
  )
  
  # Event-study version: diagnostic only.
  mod_event_r <- fixest::feols(
    y ~
      i(rel_bin_main, entrant_prior_unexposed, ref = EVENT_REF_BIN) +
      i(rel_bin_main, entrant_prior_exposed, ref = EVENT_REF_BIN) +
      i(rel_bin_main, S_enter_new, ref = EVENT_REF_BIN) +
      i(rel_bin_main, S_enter_already, ref = EVENT_REF_BIN) +
      i(rel_bin_main, S_enter_untreated, ref = EVENT_REF_BIN) |
      station_stack + month_stack,
    cluster = ~ station_id,
    data = stack_r
  )
  
  event_r <- coef_table(
    model = mod_event_r,
    model_name = "event_dynamic_static_aligned_decomposition_v4_riskset",
    radius_km = radius_km
  )
  
  event_r[, rel_bin_main := stringr::str_extract(term, "\\[[^\\]]+\\]")]
  event_r <- merge(event_r, rel_bin_main_mid, by = "rel_bin_main", all.x = TRUE)
  
  event_r[, component := data.table::fcase(
    grepl("entrant_prior_unexposed", term),
    "Entrants not previously exposed",
    grepl("entrant_prior_exposed", term),
    "Entrants previously exposed",
    grepl("S_enter_new", term),
    "Incremental exposure on entrants",
    grepl("S_enter_already", term),
    "Incremental exposure on already-treated receivers",
    grepl("S_enter_untreated", term),
    "Incremental exposure on currently untreated receivers",
    default = "Other"
  )]
  
  list(
    radius_km = radius_km,
    W = W_r,
    stack = stack_r,
    candidate_stack = candidate_stack_r,
    riskset_audit = riskset_audit_r,
    validation = validation_r,
    support = support_r,
    model = mod_decomp_r,
    coef = coef_r,
    implied = implied_r,
    legacy_comparison = legacy_compare_r$table,
    legacy_model = legacy_compare_r$legacy_model,
    cohort_model = cohort_fit_r$model,
    cohort_coef = cohort_fit_r$estimates,
    cohort_aggregation = cohort_agg_r,
    pooled_vs_cohort_aggregation = pooled_vs_agg_r,
    zero_support_pooled_audit = zero_support_pooled_audit_r,
    audit_2023_station = audit_2023_r$station_audit,
    audit_2023_leave_one_out = audit_2023_r$leave_one_out,
    audit_2023_pooled_sensitivity = audit_2023_r$pooled_sensitivity,
    effective_support_summary = effective_support_r$summary,
    effective_support_station_shares = effective_support_r$station_shares,
    pooled_effective_support_summary = pooled_effective_support_r$summary,
    pooled_effective_support_station_shares = pooled_effective_support_r$station_shares,
    pooled_effective_support_stack_shares = pooled_effective_support_r$stack_shares,
    event_model = mod_event_r,
    event = event_r
  )
}

# ------------------------------------------------------------
# 8. Run
# ------------------------------------------------------------

radius_results_12 <- lapply(RADII_KM, run_pipeline_12_radius)

source_validation_12 <- rbindlist(lapply(radius_results_12, `[[`, "validation"), fill = TRUE)
support_12 <- rbindlist(lapply(radius_results_12, `[[`, "support"), fill = TRUE)
coef_12 <- rbindlist(lapply(radius_results_12, `[[`, "coef"), fill = TRUE)
implied_12 <- rbindlist(lapply(radius_results_12, `[[`, "implied"), fill = TRUE)
event_12 <- rbindlist(lapply(radius_results_12, `[[`, "event"), fill = TRUE)
riskset_audit_12 <- rbindlist(lapply(radius_results_12, `[[`, "riskset_audit"), fill = TRUE)
legacy_comparison_12 <- rbindlist(lapply(radius_results_12, `[[`, "legacy_comparison"), fill = TRUE)
cohort_coef_12 <- rbindlist(lapply(radius_results_12, `[[`, "cohort_coef"), fill = TRUE)
cohort_aggregation_12 <- rbindlist(lapply(radius_results_12, `[[`, "cohort_aggregation"), fill = TRUE)
pooled_vs_cohort_aggregation_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "pooled_vs_cohort_aggregation"),
  fill = TRUE
)
zero_support_pooled_audit_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "zero_support_pooled_audit"),
  fill = TRUE
)
audit_2023_station_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "audit_2023_station"),
  fill = TRUE
)
audit_2023_leave_one_out_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "audit_2023_leave_one_out"),
  fill = TRUE
)
audit_2023_pooled_sensitivity_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "audit_2023_pooled_sensitivity"),
  fill = TRUE
)
effective_support_summary_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "effective_support_summary"),
  fill = TRUE
)
effective_support_station_shares_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "effective_support_station_shares"),
  fill = TRUE
)
pooled_effective_support_summary_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "pooled_effective_support_summary"),
  fill = TRUE
)
pooled_effective_support_station_shares_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "pooled_effective_support_station_shares"),
  fill = TRUE
)
pooled_effective_support_stack_shares_12 <- rbindlist(
  lapply(radius_results_12, `[[`, "pooled_effective_support_stack_shares"),
  fill = TRUE
)

fwrite(
  source_validation_12,
  file.path(TAB_DIR, "table_12_01_source_exposure_validation_by_radius.csv")
)

fwrite(
  support_12,
  file.path(TAB_DIR, "table_12_02_dynamic_static_decomp_support_by_radius_stack.csv")
)

fwrite(
  coef_12,
  file.path(TAB_DIR, "table_12_03_dynamic_static_decomp_coefficients_by_radius.csv")
)

fwrite(
  implied_12,
  file.path(TAB_DIR, "table_12_04_dynamic_static_decomp_implied_by_radius.csv")
)

fwrite(
  event_12,
  file.path(TAB_DIR, "table_12_05_dynamic_static_decomp_event_study_by_radius.csv")
)

fwrite(
  riskset_audit_12,
  file.path(TAB_DIR, "table_12_06_future_treated_risk_set_audit.csv")
)

fwrite(
  legacy_comparison_12,
  file.path(TAB_DIR, "table_12_07_legacy_vs_calendar_riskset_comparison.csv")
)

fwrite(
  cohort_coef_12,
  file.path(TAB_DIR, "table_12_08_cohort_specific_coefficients.csv")
)

fwrite(
  cohort_aggregation_12,
  file.path(TAB_DIR, "table_12_09_cohort_aggregation.csv")
)

fwrite(
  pooled_vs_cohort_aggregation_12,
  file.path(TAB_DIR, "table_12_10_pooled_vs_cohort_aggregation.csv")
)

fwrite(
  zero_support_pooled_audit_12,
  file.path(TAB_DIR, "table_12_11_structural_zero_support_pooled_audit.csv")
)

fwrite(
  audit_2023_station_12,
  file.path(TAB_DIR, "table_12_12_2023_untreated_support_station_audit.csv")
)

fwrite(
  audit_2023_leave_one_out_12,
  file.path(TAB_DIR, "table_12_13_2023_untreated_leave_one_out.csv")
)

fwrite(
  audit_2023_pooled_sensitivity_12,
  file.path(TAB_DIR, "table_12_14_2023_untreated_pooled_sensitivity.csv")
)

fwrite(
  effective_support_summary_12,
  file.path(TAB_DIR, "table_12_15_effective_identifying_support_summary.csv")
)

fwrite(
  effective_support_station_shares_12,
  file.path(TAB_DIR, "table_12_16_effective_identifying_support_station_shares.csv")
)

fwrite(
  pooled_effective_support_summary_12,
  file.path(TAB_DIR, "table_12_17_pooled_effective_identifying_support_summary.csv")
)

fwrite(
  pooled_effective_support_station_shares_12,
  file.path(TAB_DIR, "table_12_18_pooled_effective_identifying_support_station_shares.csv")
)

fwrite(
  pooled_effective_support_stack_shares_12,
  file.path(TAB_DIR, "table_12_19_pooled_effective_identifying_support_stack_shares.csv")
)

# ------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------

plot_implied <- implied_12[!is.na(estimate)]
plot_implied[, component := factor(
  component,
  levels = c(
    "Direct effect on entrants not previously exposed",
    "Direct effect on entrants previously exposed",
    "Incremental exposure on entrants",
    "Incremental exposure on already-treated receivers",
    "Incremental exposure on currently untreated receivers"
  )
)]

if (nrow(plot_implied) > 0) {
  fig_implied <- ggplot(
    plot_implied,
    aes(
      x = factor(radius_km),
      y = implied_percent,
      ymin = ci_low_percent,
      ymax = ci_high_percent
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_pointrange(position = position_dodge(width = 0.45)) +
    facet_wrap(~ component, scales = "free_y") +
    labs(
      x = "Exposure radius (km)",
      y = "Implied percent effect",
      title = "Dynamic decomposition aligned with the static treatment-exposure model",
      subtitle = "Exposure components are evaluated at the corresponding mean positive S_enter"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(FIG_DIR, "fig_12_01_dynamic_static_decomp_implied_by_radius.pdf"),
    fig_implied,
    width = 220,
    height = 145,
    units = "mm",
    device = grDevices::cairo_pdf
  )
  
  ggsave(
    file.path(FIG_DIR, "fig_12_01_dynamic_static_decomp_implied_by_radius.png"),
    fig_implied,
    width = 220,
    height = 145,
    units = "mm",
    dpi = 400,
    bg = "white"
  )
}

# Event-study figure for the selected radius.
event_plot_dt <- event_12[
  radius_km == EVENT_FIG_RADIUS &
    !is.na(rel_bin_mid) &
    component != "Other"
]

if (nrow(event_plot_dt) > 0) {
  event_plot_dt[, component := factor(
    component,
    levels = c(
      "Entrants not previously exposed",
      "Entrants previously exposed",
      "Incremental exposure on entrants",
      "Incremental exposure on already-treated receivers",
      "Incremental exposure on currently untreated receivers"
    )
  )]
  
  fig_event <- ggplot(
    event_plot_dt,
    aes(
      x = rel_bin_mid,
      y = estimate,
      ymin = estimate - 1.96 * std_error,
      ymax = estimate + 1.96 * std_error
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = -0.5, linewidth = 0.25, linetype = "dashed") +
    geom_pointrange() +
    facet_wrap(~ component, scales = "free_y") +
    scale_x_continuous(
      breaks = rel_bin_main_mid$rel_bin_mid,
      labels = rel_bin_main_mid$rel_bin_main
    ) +
    labs(
      x = "Event-time bin",
      y = "Coefficient relative to [-6,-1]",
      title = paste0("Event-study diagnostic for the dynamic decomposition, ", EVENT_FIG_RADIUS, " km"),
      subtitle = "Reference bin: [-6,-1]"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(FIG_DIR, paste0("fig_12_02_dynamic_static_decomp_event_study_", EVENT_FIG_RADIUS, "km.pdf")),
    fig_event,
    width = 225,
    height = 145,
    units = "mm",
    device = grDevices::cairo_pdf
  )
  
  ggsave(
    file.path(FIG_DIR, paste0("fig_12_02_dynamic_static_decomp_event_study_", EVENT_FIG_RADIUS, "km.png")),
    fig_event,
    width = 225,
    height = 145,
    units = "mm",
    dpi = 400,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

message("Pipeline 12 v4 completed.")
message("Tables saved to: ", TAB_DIR)
message("Figures saved to: ", FIG_DIR)

message("\nRisk-set audit (benchmark logic; zero violations expected in current ULEZ windows):")
print(riskset_audit_12[, .(
  radius_km, stack_label,
  max_m_id_future_date, min_cohort_id_future_date, gap_months_future,
  n_future_obs_at_or_after_own_treatment, n_censored_obs
)])

message("\nLegacy vs calendar-time risk-set coefficient comparison:")
print(legacy_comparison_12[, .(
  radius_km, term, legacy_estimate, riskset_estimate, estimate_difference,
  n_obs_removed
)])

message("\nPooled vs explicit cohort aggregation:")
print(pooled_vs_cohort_aggregation_12[, .(
  radius_km, term, aggregation, pooled_estimate, estimate,
  difference_from_pooled, std_error, n_cohorts_estimable, total_support_n
)])

message("\nStructural-zero-support audit for pooled coefficients:")
print(zero_support_pooled_audit_12[, .(
  radius_km, term, unsupported_stacks, pooled_full_estimate,
  pooled_supported_stacks_estimate, diff_supported_minus_full,
  target_common_nuisance_stack_specific_estimate,
  diff_nuisance_heterogeneous_minus_full
)])

message("\n2023 currently-untreated support audit: pooled sensitivity")
print(audit_2023_pooled_sensitivity_12)

message("\n2023 currently-untreated support audit: leave-one-support-station-out")
print(audit_2023_leave_one_out_12)

message("\nEffective identifying support: currently-untreated exposure by stack")
print(effective_support_summary_12[
  term == "post_untreated_S_enter",
  .(
    radius_km, stack_label, raw_support_n, raw_nonzero_regressor_stations,
    n_stations_with_nonzero_residualized_target, n_eff, max_share, top2_share, top3_share,
    max_share_station_id, residualization_status
  )
])

message("\nBenchmark 3 km / 2023 station shares for currently-untreated exposure")
print(effective_support_station_shares_12[
  radius_km == 3 &
    stack_cohort_id == month_id(as.Date("2023-08-01")) &
    term == "post_untreated_S_enter",
  .(
    station_id, raw_target_ever_nonzero, residual_ss, residual_xy,
    q_share, beta_contribution, q_rank
  )
][order(q_rank)], nrows = Inf)

message("\nPooled effective identifying support: currently-untreated exposure")
print(pooled_effective_support_summary_12[
  term == "post_untreated_S_enter",
  .(
    radius_km, pooled_estimate, beta_fwl, fwl_minus_pooled,
    raw_nonzero_regressor_stations,
    n_stations_with_nonzero_residualized_target,
    n_eff_station, max_station_share, top2_station_share,
    residual_ss_raw_support_station_share,
    residual_ss_raw_zero_station_share,
    beta_contribution_raw_support_stations,
    beta_contribution_raw_zero_stations,
    residualization_status
  )
])

message("\nPooled stack shares: currently-untreated exposure")
print(pooled_effective_support_stack_shares_12[
  term == "post_untreated_S_enter",
  .(
    radius_km, stack_label, raw_support_n, raw_nonzero_regressor_stations,
    q_share, beta_contribution, q_rank
  )
][order(radius_km, q_rank)])

print(implied_12[, .(
  radius_km,
  component,
  exposure_scale,
  estimate,
  std_error,
  implied_percent,
  ci_low_percent,
  ci_high_percent
)])


# ============================================================
# Diagnostic: receiver cohorts observed in each stack
# Pipeline 12 object name: radius_results_12
# ============================================================

if (!exists("radius_results_12")) {
  stop(
    "Object radius_results_12 not found. ",
    "Run/source pipeline 12 first, or add this diagnostic inside the pipeline after radius_results_12 is created."
  )
}

rr12 <- radius_results_12

debug_receiver_cohorts <- rbindlist(lapply(rr12, function(res) {
  
  r <- res$radius_km
  
  st <- unique(res$stack[, .(
    station_id,
    stack_label,
    stack_cohort_id,
    cohort_id,
    receiver_status,
    new_cohort,
    prior_cohort,
    future_cohort,
    never_cohort,
    S_enter,
    S_prior
  )])
  
  st[, radius_km := r]
  
  st[, receiver_cohort := fifelse(
    cohort_id == month_id(as.Date("2019-04-01")),
    "2019 Central ULEZ",
    fifelse(
      cohort_id == month_id(as.Date("2021-10-01")),
      "2021 Inner London",
      fifelse(
        cohort_id == month_id(as.Date("2023-08-01")),
        "2023 London-wide",
        "Never-treated"
      )
    )
  )]
  
  st[, receiver_status_check := fifelse(
    new_cohort == 1,
    "newly-treated",
    fifelse(
      prior_cohort == 1,
      "already-treated",
      fifelse(
        future_cohort == 1,
        "not-yet-treated",
        "never-treated"
      )
    )
  )]
  
  st[
    ,
    .(
      n_stations = uniqueN(station_id),
      n_S_enter_positive = uniqueN(station_id[S_enter > EXPOSURE_EPS]),
      mean_S_enter = mean(S_enter, na.rm = TRUE),
      mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE)
    ),
    by = .(
      radius_km,
      stack_label,
      receiver_cohort,
      receiver_status_check
    )
  ]
  
}), fill = TRUE)

fwrite(
  debug_receiver_cohorts,
  file.path(TAB_DIR, "table_12_06_receiver_cohorts_observed_by_stack.csv")
)

debug_receiver_cohorts[
  radius_km == 6 &
    stack_label == "Inner London expansion, Oct. 2021"
]

# ============================================================
# 11. Benchmark outputs: 3 km
# ============================================================

BENCHMARK_RADIUS_KM <- 3

if (!exists("radius_results_12")) {
  stop("radius_results_12 not found. Run pipeline 12 first.")
}

idx_bench <- which(
  vapply(radius_results_12, function(x) x$radius_km, numeric(1)) ==
    BENCHMARK_RADIUS_KM
)

if (length(idx_bench) != 1) {
  stop("Could not find benchmark radius: ", BENCHMARK_RADIUS_KM, " km.")
}

res_bench <- radius_results_12[[idx_bench]]
stack_bench <- data.table::copy(res_bench$stack)

# ------------------------------------------------------------
# 11.1 Benchmark tables
# ------------------------------------------------------------

validation_bench <- source_validation_12[radius_km == BENCHMARK_RADIUS_KM]
support_bench <- support_12[radius_km == BENCHMARK_RADIUS_KM]
coef_bench <- coef_12[radius_km == BENCHMARK_RADIUS_KM]
implied_bench <- implied_12[radius_km == BENCHMARK_RADIUS_KM]
event_bench <- event_12[radius_km == BENCHMARK_RADIUS_KM]

data.table::fwrite(
  validation_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_01a_source_exposure_validation_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

data.table::fwrite(
  support_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_02a_dynamic_static_decomp_support_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

data.table::fwrite(
  coef_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_03a_dynamic_static_decomp_coefficients_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

data.table::fwrite(
  implied_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_04a_dynamic_static_decomp_implied_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

data.table::fwrite(
  event_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_05a_dynamic_static_decomp_event_study_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

# ------------------------------------------------------------
# 11.2 Receiver cohorts observed in each stack, benchmark radius
# ------------------------------------------------------------

receiver_cohorts_bench <- unique(stack_bench[, .(
  station_id,
  stack_label,
  stack_cohort_id,
  cohort_id,
  receiver_status,
  new_cohort,
  prior_cohort,
  future_cohort,
  never_cohort,
  S_enter,
  S_prior
)])

receiver_cohorts_bench[, radius_km := BENCHMARK_RADIUS_KM]

receiver_cohorts_bench[, receiver_cohort := fifelse(
  cohort_id == month_id(as.Date("2019-04-01")),
  "2019 Central ULEZ",
  fifelse(
    cohort_id == month_id(as.Date("2021-10-01")),
    "2021 Inner London",
    fifelse(
      cohort_id == month_id(as.Date("2023-08-01")),
      "2023 London-wide",
      "Never-treated"
    )
  )
)]

receiver_cohorts_bench[, receiver_status_check := fifelse(
  new_cohort == 1,
  "newly-treated",
  fifelse(
    prior_cohort == 1,
    "already-treated",
    fifelse(
      future_cohort == 1,
      "not-yet-treated",
      "never-treated"
    )
  )
)]

receiver_cohorts_observed_bench <- receiver_cohorts_bench[
  ,
  .(
    n_stations = uniqueN(station_id),
    n_S_enter_positive = uniqueN(station_id[S_enter > EXPOSURE_EPS]),
    mean_S_enter = mean(S_enter, na.rm = TRUE),
    mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE),
    mean_S_prior = mean(S_prior, na.rm = TRUE),
    mean_S_prior_positive = mean(S_prior[S_prior > EXPOSURE_EPS], na.rm = TRUE)
  ),
  by = .(
    radius_km,
    stack_label,
    receiver_cohort,
    receiver_status_check
  )
][order(stack_label, receiver_status_check, receiver_cohort)]

data.table::fwrite(
  receiver_cohorts_observed_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_06_receiver_cohorts_observed_by_stack_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

# ------------------------------------------------------------
# 11.3 Clean support table for benchmark radius
# ------------------------------------------------------------

support_clean_bench <- support_bench[, .(
  radius_km,
  stack_label,
  
  n_obs,
  n_stations,
  
  newly_treated_receivers_observed =
    n_new_stations,
  
  newly_treated_receivers_previously_unexposed =
    n_new_prior_unexposed_stations,
  
  newly_treated_receivers_previously_exposed =
    n_new_prior_exposed_stations,
  
  newly_treated_receivers_exposed_to_entering_cohort =
    n_new_S_enter_pos_stations,
  
  already_treated_receivers_exposed_to_entering_cohort =
    n_already_S_enter_pos_stations,
  
  currently_untreated_receivers_exposed_to_entering_cohort =
    n_untreated_S_enter_pos_stations,
  
  not_yet_treated_receivers_exposed_to_entering_cohort =
    n_future_S_enter_pos_stations,
  
  never_treated_receivers_exposed_to_entering_cohort =
    n_never_S_enter_pos_stations,
  
  mean_S_prior_new,
  mean_S_prior_new_positive,
  
  mean_S_enter_new,
  mean_S_enter_new_positive,
  
  mean_S_enter_already,
  mean_S_enter_already_positive,
  
  mean_S_enter_untreated,
  mean_S_enter_untreated_positive
)]

data.table::fwrite(
  support_clean_bench,
  file.path(
    TAB_DIR,
    paste0("table_12_02b_dynamic_static_decomp_support_clean_", BENCHMARK_RADIUS_KM, "km.csv")
  )
)

# ------------------------------------------------------------
# 11.4 Main implied-effect figure for benchmark radius
# ------------------------------------------------------------

implied_plot_bench <- data.table::copy(implied_bench)

implied_plot_bench[, component_clean := data.table::fcase(
  component == "Direct effect on entrants not previously exposed",
  "Direct effect on entrants\nnot previously exposed",
  
  component == "Direct effect on entrants previously exposed",
  "Direct effect on entrants\npreviously exposed",
  
  component == "Incremental exposure on entrants",
  "Incremental exposure\non entrants",
  
  component == "Incremental exposure on already-treated receivers",
  "Incremental exposure\non already-treated receivers",
  
  component == "Incremental exposure on currently untreated receivers",
  "Incremental exposure\non untreated receivers",
  
  default = as.character(component)
)]

implied_plot_bench[, component_clean := factor(
  component_clean,
  levels = c(
    "Direct effect on entrants\nnot previously exposed",
    "Direct effect on entrants\npreviously exposed",
    "Incremental exposure\non entrants",
    "Incremental exposure\non already-treated receivers",
    "Incremental exposure\non untreated receivers"
  )
)]

fig_implied_bench <- ggplot(
  implied_plot_bench,
  aes(
    x = component_clean,
    y = implied_percent,
    ymin = ci_low_percent,
    ymax = ci_high_percent
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_pointrange() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Implied percent effect",
    title = paste0(
      "Dynamic decomposition aligned with the static model, ",
      BENCHMARK_RADIUS_KM,
      " km benchmark"
    ),
    subtitle = "Exposure components are evaluated at the corresponding mean positive S_enter"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.5)
  )

ggsave(
  file.path(
    FIG_DIR,
    paste0("fig_12_01a_dynamic_static_decomp_implied_", BENCHMARK_RADIUS_KM, "km.pdf")
  ),
  fig_implied_bench,
  width = 175,
  height = 115,
  units = "mm",
  device = grDevices::cairo_pdf
)

ggsave(
  file.path(
    FIG_DIR,
    paste0("fig_12_01a_dynamic_static_decomp_implied_", BENCHMARK_RADIUS_KM, "km.png")
  ),
  fig_implied_bench,
  width = 175,
  height = 115,
  units = "mm",
  dpi = 400,
  bg = "white"
)

# ------------------------------------------------------------
# 11.5 Event-study figure for benchmark radius
# ------------------------------------------------------------

event_plot_bench <- event_bench[
  !is.na(rel_bin_mid) &
    component != "Other"
]

event_plot_bench[, component := factor(
  component,
  levels = c(
    "Entrants not previously exposed",
    "Entrants previously exposed",
    "Incremental exposure on entrants",
    "Incremental exposure on already-treated receivers",
    "Incremental exposure on currently untreated receivers"
  )
)]

event_plot_bench[, component_wrapped := stringr::str_wrap(
  as.character(component),
  width = 32
)]

fig_event_bench <- ggplot(
  event_plot_bench,
  aes(
    x = rel_bin_mid,
    y = estimate,
    ymin = estimate - 1.96 * std_error,
    ymax = estimate + 1.96 * std_error
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = -0.5, linewidth = 0.25, linetype = "dashed") +
  geom_pointrange() +
  facet_wrap(~ component_wrapped, scales = "free_y") +
  scale_x_continuous(
    breaks = rel_bin_main_mid$rel_bin_mid,
    labels = rel_bin_main_mid$rel_bin_main
  ) +
  labs(
    x = "Event-time bin",
    y = "Coefficient relative to [-6,-1]",
    title = paste0(
      "Event-study diagnostic for the dynamic decomposition, ",
      BENCHMARK_RADIUS_KM,
      " km benchmark"
    ),
    subtitle = "Reference bin: [-6,-1]"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = grid::unit(4, "mm")
  )

ggsave(
  file.path(
    FIG_DIR,
    paste0("fig_12_02a_dynamic_static_decomp_event_study_", BENCHMARK_RADIUS_KM, "km.pdf")
  ),
  fig_event_bench,
  width = 225,
  height = 145,
  units = "mm",
  device = grDevices::cairo_pdf
)

ggsave(
  file.path(
    FIG_DIR,
    paste0("fig_12_02a_dynamic_static_decomp_event_study_", BENCHMARK_RADIUS_KM, "km.png")
  ),
  fig_event_bench,
  width = 225,
  height = 145,
  units = "mm",
  dpi = 400,
  bg = "white"
)

# ============================================================
# 12. Focused event-study:
# Exposure among currently untreated exposed units
# Benchmark radius: 3 km
# ============================================================

BENCHMARK_RADIUS_KM <- 3

idx_bench <- which(
  vapply(radius_results_12, function(x) x$radius_km, numeric(1)) ==
    BENCHMARK_RADIUS_KM
)

if (length(idx_bench) != 1) {
  stop("Could not find benchmark radius: ", BENCHMARK_RADIUS_KM, " km.")
}

stack_bench <- data.table::copy(radius_results_12[[idx_bench]]$stack)

# Standardized language:
# untreated_receiver = currently untreated unit in the stack:
# not-yet-treated or never-treated.
stack_bench[, currently_untreated_unit := untreated_receiver]

stack_bench[, currently_untreated_exposed :=
              as.integer(currently_untreated_unit == 1 & S_enter > EXPOSURE_EPS)
]

# Focused sample: units not directly treated in the stack.
untreated_es_dt <- stack_bench[currently_untreated_unit == 1]

# Safety support check by stack and event bin.
support_untreated_es <- untreated_es_dt[
  ,
  .(
    n_obs = .N,
    n_stations = uniqueN(station_id),
    n_exposed_obs = sum(S_enter > EXPOSURE_EPS),
    n_exposed_stations = uniqueN(station_id[S_enter > EXPOSURE_EPS]),
    n_unexposed_stations = uniqueN(station_id[S_enter <= EXPOSURE_EPS]),
    mean_S_enter = mean(S_enter, na.rm = TRUE),
    mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE)
  ),
  by = .(stack_label, rel_bin_main)
][order(stack_label, rel_bin_main)]

data.table::fwrite(
  support_untreated_es,
  file.path(
    TAB_DIR,
    paste0("table_12_09_focused_untreated_exposure_event_support_",
           BENCHMARK_RADIUS_KM, "km.csv")
  )
)

untreated_es_dt_2019_2021 <- untreated_es_dt[
  stack_label != "London-wide expansion, Aug. 2023"
]

mod_es_untreated_S_2019_2021 <- fixest::feols(
  as.formula(
    paste0(
      "y ~ i(rel_bin_main, S_enter, ref = '[-6,-1]') | ",
      "station_stack + month_stack"
    )
  ),
  data = untreated_es_dt_2019_2021,
  cluster = ~ station_id
)

mod_post_untreated_2019_2021 <- fixest::feols(
  as.formula(
    paste0(
      "y ~ post_untreated_S_enter | ",
      "station_stack + month_stack"
    )
  ),
  data = stack_bench[
    untreated_receiver == 1 &
      stack_label != "London-wide expansion, Aug. 2023"
  ],
  cluster = ~ station_id
)


untreated_2019_comp <- unique(stack_bench[
  stack_label == "Central ULEZ, Apr. 2019" &
    untreated_receiver == 1,
  .(
    station_id,
    cohort_id,
    cohort_label,
    S_enter,
    exposed_to_2019 = S_enter > EXPOSURE_EPS
  )
])

untreated_2019_comp[
  ,
  .(
    n_stations = uniqueN(station_id),
    n_exposed = uniqueN(station_id[exposed_to_2019]),
    mean_S_enter = mean(S_enter, na.rm = TRUE),
    mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE)
  ),
  by = cohort_label
][order(cohort_label)]


# Stack 2019, only currently untreated units
dt_2019_untreated <- data.table::copy(stack_bench[
  stack_cohort_id == month_id(as.Date("2019-04-01")) &
    untreated_receiver == 1
])

dt_2019_untreated[, future_group := data.table::fcase(
  cohort_id == month_id(as.Date("2021-10-01")),
  "Future 2021",
  
  cohort_id == month_id(as.Date("2023-08-01")),
  "Future 2023",
  
  !is.finite(cohort_id),
  "Never treated",
  
  default = "Other"
)]

dt_2019_untreated[, future_group := factor(
  future_group,
  levels = c("Future 2021", "Future 2023", "Never treated", "Other")
)]

dt_2019_untreated[
  ,
  .(
    n_stations = uniqueN(station_id),
    n_exposed_stations = uniqueN(station_id[S_enter > EXPOSURE_EPS]),
    n_unexposed_stations = uniqueN(station_id[S_enter <= EXPOSURE_EPS]),
    mean_S_enter = mean(S_enter, na.rm = TRUE),
    mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE)
  ),
  by = future_group
][order(future_group)]

dt_2019_untreated[, S_enter_future_2021 :=
                    fifelse(future_group == "Future 2021", S_enter, 0)
]

dt_2019_untreated[, S_enter_future_2023 :=
                    fifelse(future_group == "Future 2023", S_enter, 0)
]

dt_2019_untreated[, S_enter_never :=
                    fifelse(future_group == "Never treated", S_enter, 0)
]

mod_2019_future_group <- fixest::feols(
  y ~
    i(rel_bin_main, S_enter_future_2021, ref = "[-6,-1]") +
    i(rel_bin_main, S_enter_future_2023, ref = "[-6,-1]") +
    i(rel_bin_main, S_enter_never, ref = "[-6,-1]") |
    station_id + m_id,
  data = dt_2019_untreated,
  cluster = ~ station_id
)

summary(mod_2019_future_group)


dt_2019_untreated[, announcement_period := data.table::fcase(
  month_date < as.Date("2018-06-01"),
  "before_2021_expansion_confirmation",
  
  month_date >= as.Date("2018-06-01") &
    month_date < as.Date("2019-04-01"),
  "after_2021_confirmation_before_2019_implementation",
  
  month_date >= as.Date("2019-04-01") &
    month_date <= as.Date("2020-06-30"),
  "after_2019_implementation",
  
  default = NA_character_
)]

dt_2019_untreated <- dt_2019_untreated[!is.na(announcement_period)]

dt_2019_untreated[, announcement_period := factor(
  announcement_period,
  levels = c(
    "before_2021_expansion_confirmation",
    "after_2021_confirmation_before_2019_implementation",
    "after_2019_implementation"
  )
)]

mod_announcement_2019 <- fixest::feols(
  y ~ i(
    announcement_period,
    S_enter,
    ref = "before_2021_expansion_confirmation"
  ) |
    station_id + m_id,
  data = dt_2019_untreated,
  cluster = ~ station_id
)

summary(mod_announcement_2019)


b <- coef(mod_announcement_2019)
V <- vcov(mod_announcement_2019)

nm_ant <- grep(
  "after_2021_confirmation_before_2019_implementation",
  names(b),
  value = TRUE
)

nm_post <- grep(
  "after_2019_implementation",
  names(b),
  value = TRUE
)

diff_est <- b[nm_post] - b[nm_ant]

diff_se <- sqrt(
  V[nm_post, nm_post] +
    V[nm_ant, nm_ant] -
    2 * V[nm_post, nm_ant]
)

diff_t <- diff_est / diff_se
diff_p <- 2 * pt(
  -abs(diff_t),
  df = data.table::uniqueN(dt_2019_untreated$station_id) - 1
)

data.table::data.table(
  contrast = "after_2019_implementation minus anticipation_window",
  estimate = diff_est,
  std_error = diff_se,
  t_value = diff_t,
  p_value = diff_p
)

untreated_2021_comp <- unique(stack_bench[
  stack_cohort_id == month_id(as.Date("2021-10-01")) &
    untreated_receiver == 1,
  .(
    station_id,
    cohort_id,
    cohort_label,
    S_enter,
    exposed_to_2021 = S_enter > EXPOSURE_EPS
  )
])

untreated_2021_comp[
  ,
  .(
    n_stations = uniqueN(station_id),
    n_exposed = uniqueN(station_id[exposed_to_2021]),
    mean_S_enter = mean(S_enter, na.rm = TRUE),
    mean_S_enter_positive = mean(S_enter[S_enter > EXPOSURE_EPS], na.rm = TRUE)
  ),
  by = cohort_label
][order(cohort_label)]


dt_2021_untreated <- data.table::copy(stack_bench[
  stack_cohort_id == month_id(as.Date("2021-10-01")) &
    untreated_receiver == 1
])

dt_2021_untreated[, period_2023_ant := data.table::fcase(
  rel_k >= -24 & rel_k <= -13,
  "early_pre_2021",
  
  rel_k >= -12 & rel_k <= -7,
  "mid_pre_2021",
  
  rel_k >= -6 & rel_k <= -1,
  "immediate_pre_2021",
  
  month_date >= as.Date("2021-10-01") &
    month_date < as.Date("2022-05-01"),
  "after_2021_before_2023_consultation",
  
  month_date >= as.Date("2022-05-01") &
    month_date < as.Date("2022-11-01"),
  "during_after_2023_consultation",
  
  month_date >= as.Date("2022-11-01") &
    month_date < as.Date("2023-08-01"),
  "after_2023_confirmation_before_implementation",
  
  default = NA_character_
)]

dt_2021_untreated <- dt_2021_untreated[!is.na(period_2023_ant)]

dt_2021_untreated[, period_2023_ant := factor(
  period_2023_ant,
  levels = c(
    "early_pre_2021",
    "mid_pre_2021",
    "immediate_pre_2021",
    "after_2021_before_2023_consultation",
    "during_after_2023_consultation",
    "after_2023_confirmation_before_implementation"
  )
)]

mod_2021_periods <- fixest::feols(
  y ~ i(period_2023_ant, S_enter, ref = "mid_pre_2021") |
    station_id + m_id,
  data = dt_2021_untreated,
  cluster = ~ station_id
)

summary(mod_2021_periods)

b <- coef(mod_2021_periods)
V <- vcov(mod_2021_periods)

make_contrast <- function(name_a, name_b, label) {
  est <- b[name_b] - b[name_a]
  se <- sqrt(
    V[name_b, name_b] +
      V[name_a, name_a] -
      2 * V[name_b, name_a]
  )
  tval <- est / se
  pval <- 2 * pt(
    -abs(tval),
    df = data.table::uniqueN(dt_2021_untreated$station_id) - 1
  )
  
  data.table::data.table(
    contrast = label,
    estimate = est,
    std_error = se,
    t_value = tval,
    p_value = pval
  )
}

nm_early <- grep("early_pre_2021", names(b), value = TRUE)
nm_immediate <- grep("immediate_pre_2021", names(b), value = TRUE)
nm_post_pre2023 <- grep("after_2021_before_2023_consultation", names(b), value = TRUE)
nm_consult <- grep("during_after_2023_consultation", names(b), value = TRUE)
nm_confirm <- grep("after_2023_confirmation_before_implementation", names(b), value = TRUE)

rbind(
  make_contrast(
    nm_early,
    nm_immediate,
    "immediate_pre_2021 minus early_pre_2021"
  ),
  make_contrast(
    nm_immediate,
    nm_post_pre2023,
    "after_2021_before_2023_consultation minus immediate_pre_2021"
  ),
  make_contrast(
    nm_post_pre2023,
    nm_consult,
    "during_2023_consultation minus after_2021_before_2023_consultation"
  ),
  make_contrast(
    nm_post_pre2023,
    nm_confirm,
    "after_2023_confirmation minus after_2021_before_2023_consultation"
  )
)


dt_2019_untreated <- data.table::copy(stack_bench[
  stack_cohort_id == month_id(as.Date("2019-04-01")) &
    untreated_receiver == 1
])

dt_2019_untreated[, F2021 := as.integer(
  cohort_id == month_id(as.Date("2021-10-01"))
)]

dt_2019_untreated[, G_enter := as.integer(S_enter > EXPOSURE_EPS)]

dt_2019_untreated[, period_2019 := data.table::fcase(
  month_date < as.Date("2018-06-01"),
  "before_2021_confirmation",
  
  month_date >= as.Date("2018-06-01") &
    month_date < as.Date("2019-04-01"),
  "after_2021_confirmation_before_2019",
  
  month_date >= as.Date("2019-04-01"),
  "after_2019_implementation",
  
  default = NA_character_
)]

dt_2019_untreated <- dt_2019_untreated[!is.na(period_2019)]

dt_2019_untreated[, period_2019 := factor(
  period_2019,
  levels = c(
    "before_2021_confirmation",
    "after_2021_confirmation_before_2019",
    "after_2019_implementation"
  )
)]

mod_2019_anticipation_vs_exposure <- fixest::feols(
  y ~
    i(period_2019, F2021, ref = "before_2021_confirmation") +
    i(period_2019, S_enter, ref = "before_2021_confirmation") |
    station_id + m_id,
  data = dt_2019_untreated,
  cluster = ~ station_id
)

summary(mod_2019_anticipation_vs_exposure)

b <- coef(mod_2019_anticipation_vs_exposure)
V <- vcov(mod_2019_anticipation_vs_exposure)

nm_ant_S <- grep(
  "after_2021_confirmation_before_2019.*S_enter",
  names(b),
  value = TRUE
)

nm_post_S <- grep(
  "after_2019_implementation.*S_enter",
  names(b),
  value = TRUE
)

diff_est <- b[nm_post_S] - b[nm_ant_S]

diff_se <- sqrt(
  V[nm_post_S, nm_post_S] +
    V[nm_ant_S, nm_ant_S] -
    2 * V[nm_post_S, nm_ant_S]
)

diff_t <- diff_est / diff_se
diff_p <- 2 * pt(
  -abs(diff_t),
  df = data.table::uniqueN(dt_2019_untreated$station_id) - 1
)

data.table::data.table(
  contrast = "additional spatial-exposure change after 2019 implementation",
  estimate = diff_est,
  std_error = diff_se,
  t_value = diff_t,
  p_value = diff_p
)

dt_2019_f2021 <- dt_2019_untreated[F2021 == 1]

mod_2019_f2021_exposed_vs_unexposed <- fixest::feols(
  y ~ i(period_2019, G_enter, ref = "before_2021_confirmation") |
    station_id + m_id,
  data = dt_2019_f2021,
  cluster = ~ station_id
)

summary(mod_2019_f2021_exposed_vs_unexposed)


mod_2019_f2021_S <- fixest::feols(
  y ~ i(period_2019, S_enter, ref = "before_2021_confirmation") |
    station_id + m_id,
  data = dt_2019_f2021,
  cluster = ~ station_id
)

summary(mod_2019_f2021_S)

# ============================================================
# Exportações canônicas do paper
# ============================================================

# Tabela de suporte dinâmico do benchmark de 3 km.
stack_order <- c(
  "Central ULEZ, Apr. 2019",
  "Inner London expansion, Oct. 2021",
  "London-wide expansion, Aug. 2023"
)
support_paper <- data.table::copy(support_clean_bench)
support_paper[, stack_label := factor(stack_label, levels = stack_order)]
data.table::setorder(support_paper, stack_label)

support_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Support for the dynamic treatment-exposure decomposition, 0--3 km benchmark}",
  "\\label{tab:ulez_dynamic_support_3km}",
  "\\small",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Stack & New, prior unexposed & New, prior exposed & New exposed & Already treated exposed & Currently untreated exposed \\\\",
  "\\midrule",
  paste0(
    as.character(support_paper$stack_label), " & ",
    support_paper$newly_treated_receivers_previously_unexposed, " & ",
    support_paper$newly_treated_receivers_previously_exposed, " & ",
    support_paper$newly_treated_receivers_exposed_to_entering_cohort, " & ",
    support_paper$already_treated_receivers_exposed_to_entering_cohort, " & ",
    support_paper$currently_untreated_receivers_exposed_to_entering_cohort, " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  "\\item Notes: Entries are numbers of monitoring stations supporting each component in the corresponding stack. Exposure refers to positive exposure to the entering cohort under the 0--3 km row-standardized mapping. Future-treated stations enter the currently untreated group only while they remain untreated in calendar time.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(support_lines, file.path(TABLE_DIR, "table_ulez_dynamic_support_3km.tex"))

# Effective identifying support: only the component reported in the paper.
effective_support_paper <- effective_support_summary_12[
  radius_km == 3 & term == "post_untreated_S_enter",
  .(
    stack_label,
    nominal_exposed_receivers = raw_support_n,
    raw_target_stations = raw_nonzero_regressor_stations,
    nonzero_residualized_target_stations = n_stations_with_nonzero_residualized_target,
    n_eff,
    max_share,
    top2_share,
    max_share_station_id
  )
]
pooled_effective_support_paper <- pooled_effective_support_summary_12[
  radius_km == 3 & term == "post_untreated_S_enter",
  .(
    pooled_estimate,
    beta_fwl,
    raw_nonzero_regressor_stations,
    nonzero_residualized_target_stations = n_stations_with_nonzero_residualized_target,
    n_eff = n_eff_station,
    max_share = max_station_share,
    top2_share = top2_station_share
  )
]

# Announcement-calendar diagnostic used in the appendix.
# Re-estimate in a compact block so the reported table is not tied to console output.
dt19 <- data.table::copy(stack_bench[
  stack_cohort_id == month_id(as.Date("2019-04-01")) & untreated_receiver == 1
])
dt19[, F2021 := as.integer(cohort_id == month_id(as.Date("2021-10-01")))]
dt19[, period_2019 := data.table::fcase(
  month_date < as.Date("2018-06-01"), "before_2021_confirmation",
  month_date >= as.Date("2018-06-01") & month_date < as.Date("2019-04-01"),
    "after_2021_confirmation_before_2019",
  month_date >= as.Date("2019-04-01"), "after_2019_implementation",
  default = NA_character_
)]
dt19 <- dt19[!is.na(period_2019)]
dt19[, period_2019 := factor(period_2019, levels = c(
  "before_2021_confirmation",
  "after_2021_confirmation_before_2019",
  "after_2019_implementation"
))]
mod19 <- fixest::feols(
  y ~ i(period_2019, F2021, ref = "before_2021_confirmation") +
    i(period_2019, S_enter, ref = "before_2021_confirmation") |
    station_id + m_id,
  data = dt19, cluster = ~ station_id
)
ct19 <- data.table::as.data.table(fixest::coeftable(mod19), keep.rownames = "term")
data.table::setnames(
  ct19,
  old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "z value", "Pr(>|z|)"),
  new = c("estimate", "std_error", "statistic", "p_value", "statistic", "p_value"),
  skip_absent = TRUE
)
get19 <- function(pattern, label) {
  z <- ct19[grepl(pattern, term)]
  if (nrow(z) != 1L) stop("Termo inesperado no diagnóstico 2019: ", pattern)
  z[, .(component = label, estimate, std_error, p_value)]
}
ann19 <- data.table::rbindlist(list(
  get19("after_2021_confirmation_before_2019.*F2021", "$A_tF_i^{2021}$"),
  get19("after_2019_implementation.*F2021", "$P_tF_i^{2021}$"),
  get19("after_2021_confirmation_before_2019.*S_enter", "$A_tS_i^{\\mathrm{enter},2019}$"),
  get19("after_2019_implementation.*S_enter", "$P_tS_i^{\\mathrm{enter},2019}$")
))
b19 <- stats::coef(mod19)
V19 <- stats::vcov(mod19)
nmA <- grep("after_2021_confirmation_before_2019.*S_enter", names(b19), value = TRUE)
nmP <- grep("after_2019_implementation.*S_enter", names(b19), value = TRUE)
diff19 <- unname(b19[nmP] - b19[nmA])
se19 <- sqrt(unname(V19[nmP, nmP] + V19[nmA, nmA] - 2 * V19[nmP, nmA]))
p19 <- 2 * stats::pt(-abs(diff19 / se19), df = data.table::uniqueN(dt19$station_id) - 1)
ann19 <- data.table::rbindlist(list(
  ann19,
  data.table::data.table(component = "$\\kappa_P-\\kappa_A$", estimate = diff19, std_error = se19, p_value = p19)
))

# Exposure-timing contrasts used for the 2021 stack.
dt21 <- data.table::copy(stack_bench[
  stack_cohort_id == month_id(as.Date("2021-10-01")) & untreated_receiver == 1
])
dt21[, period_2023_ant := data.table::fcase(
  rel_k >= -24 & rel_k <= -13, "early_pre_2021",
  rel_k >= -12 & rel_k <= -7, "mid_pre_2021",
  rel_k >= -6 & rel_k <= -1, "immediate_pre_2021",
  month_date >= as.Date("2021-10-01") & month_date < as.Date("2022-05-01"),
    "after_2021_before_2023_consultation",
  month_date >= as.Date("2022-05-01") & month_date < as.Date("2022-11-01"),
    "during_after_2023_consultation",
  month_date >= as.Date("2022-11-01") & month_date < as.Date("2023-08-01"),
    "after_2023_confirmation_before_implementation",
  default = NA_character_
)]
dt21 <- dt21[!is.na(period_2023_ant)]
dt21[, period_2023_ant := factor(period_2023_ant, levels = c(
  "early_pre_2021", "mid_pre_2021", "immediate_pre_2021",
  "after_2021_before_2023_consultation", "during_after_2023_consultation",
  "after_2023_confirmation_before_implementation"
))]
mod21 <- fixest::feols(
  y ~ i(period_2023_ant, S_enter, ref = "mid_pre_2021") | station_id + m_id,
  data = dt21, cluster = ~ station_id
)
b21 <- stats::coef(mod21)
V21 <- stats::vcov(mod21)
find21 <- function(txt) grep(txt, names(b21), value = TRUE)
contrast21 <- function(a, b, label) {
  est <- unname(b21[b] - b21[a])
  se <- sqrt(unname(V21[b, b] + V21[a, a] - 2 * V21[b, a]))
  data.table::data.table(
    contrast = label, estimate = est, std_error = se,
    p_value = 2 * stats::pt(-abs(est / se), df = data.table::uniqueN(dt21$station_id) - 1)
  )
}
nm_early <- find21("early_pre_2021")
nm_immediate <- find21("immediate_pre_2021")
nm_post <- find21("after_2021_before_2023_consultation")
nm_consult <- find21("during_after_2023_consultation")
nm_confirm <- find21("after_2023_confirmation_before_implementation")
timing21 <- data.table::rbindlist(list(
  contrast21(nm_early, nm_immediate, "Immediate pre-2021 minus early pre-2021"),
  contrast21(nm_immediate, nm_post, "After 2021, before 2023 consultation minus immediate pre-2021"),
  contrast21(nm_post, nm_consult, "During 2023 consultation minus after 2021, before 2023 consultation"),
  contrast21(nm_post, nm_confirm, "After 2023 confirmation minus after 2021, before 2023 consultation")
))

# Tabelas LaTeX dos dois calendários.
fmt3 <- function(x) formatC(x, digits = 3, format = "f")
writeLines(c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Announcement-calendar diagnostic for the 2019 stack}",
  "\\label{tab:ulez_announcement_2019}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}", "\\toprule", "Component & Estimate & Std. error & $p$-value \\\\", "\\midrule",
  paste0(ann19$component[1:4], " & ", fmt3(ann19$estimate[1:4]), " & ", fmt3(ann19$std_error[1:4]), " & ", fmt3(ann19$p_value[1:4]), " \\\\"),
  "\\midrule",
  paste0(ann19$component[5], " & ", fmt3(ann19$estimate[5]), " & ", fmt3(ann19$std_error[5]), " & ", fmt3(ann19$p_value[5]), " \\\\"),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
  "\\item \\emph{Notes:} The table reports the announcement-calendar diagnostic for currently untreated observations in the 2019 stack. $F_i^{2021}$ indicates stations belonging to the future 2021 expansion cohort. $A_t$ denotes the period after the 2021 expansion had been confirmed but before the April 2019 implementation. $P_t$ denotes the period after the April 2019 implementation. Standard errors are clustered by monitoring station. The contrast $\\kappa_P-\\kappa_A$ measures the additional exposure-correlated change after implementation relative to the pre-implementation announcement period.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_announcement_2019.tex"))

writeLines(c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Exposure-timing contrasts for the 2021 stack}",
  "\\label{tab:ulez_exposure_timing_2021}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabularx}{\\textwidth}{Xccc}", "\\toprule", "Contrast & Estimate & Std. error & $p$-value \\\\", "\\midrule",
  paste0(timing21$contrast, " & ", fmt3(timing21$estimate), " & ", fmt3(timing21$std_error), " & ", fmt3(timing21$p_value), " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}", "\\footnotesize",
  "\\item \\emph{Notes:} The table reports calendar-based contrasts for the exposure gradient among currently untreated observations in the 2021 stack. The currently untreated exposed stations in this stack belong to the future 2023 expansion cohort. Standard errors are clustered by monitoring station. The contrasts are designed to distinguish changes associated with the 2021 implementation from changes associated with the formal consultation and confirmation of the 2023 London-wide expansion.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_exposure_timing_2021.tex"))

# Copia apenas as duas figuras dinâmicas que aparecem no manuscrito.
file.copy(
  file.path(FIG_DIR, "fig_12_01a_dynamic_static_decomp_implied_3km.png"),
  file.path(FIGURE_DIR, "fig_12_01a_dynamic_static_decomp_implied_3km.png"),
  overwrite = TRUE
)
file.copy(
  file.path(FIG_DIR, "fig_12_02a_dynamic_static_decomp_event_study_3km.png"),
  file.path(FIGURE_DIR, "fig_12_02a_dynamic_static_decomp_event_study_3km.png"),
  overwrite = TRUE
)

saveRDS(
  list(
    support = support_paper,
    coefficients = coef_bench,
    implied = implied_bench,
    event_study = event_bench,
    effective_support = effective_support_paper,
    pooled_effective_support = pooled_effective_support_paper,
    announcement_2019 = ann19,
    timing_2021 = timing21,
    stack = stack_bench
  ),
  DYNAMIC_RESULTS_RDS
)
message("[03] Resultados dinâmicos concluídos.")
