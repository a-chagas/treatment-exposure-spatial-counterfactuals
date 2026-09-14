# ============================================================
# 04_boundary_network_results.R
# Policy-boundary exposure e predetermined-monitor network
# ============================================================
#
# Este módulo reproduz os diagnósticos incorporados após o segundo parecer:
# (i) representação de exposição baseada diretamente na área regulada e
# (ii) restrição da rede de monitoramento a estações abertas antes de
# 13 de fevereiro de 2013.
# ============================================================

message("[04] Estimando diagnósticos de policy geography e monitoring network...")

OUT_DIR <- file.path(DERIVED_DIR, "boundary_network_core")
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "rds")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

RADII_KM <- c(3)
EXPOSURE_EPS <- 1e-10
TARGET_CRS <- 27700L
EXPANSION_YEARS <- c(2019L, 2021L, 2023L)
L_PRE <- 24L
H_POST <- 15L

# 2. Helpers
# ------------------------------------------------------------

clean_names_simple <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

month_id <- function(date) {
  lubridate::year(date) * 12L + lubridate::month(date)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) return(NA_real_)
  stats::cor(x[ok], y[ok])
}

pct_effect <- function(beta) {
  100 * (exp(beta) - 1)
}

to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "true", "t", "yes", "y", "sim")] <- 1L
  out[z %in% c("0", "false", "f", "no", "n", "nao", "não")] <- 0L
  out
}

parse_monitor_date_13 <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_

  parsed <- suppressWarnings(lubridate::parse_date_time(
    z,
    orders = c(
      "Ymd", "Y-m-d", "Y/m/d",
      "dmy", "d/m/Y", "d-m-Y",
      "mdy", "m/d/Y", "m-d-Y",
      "Y", "Ym", "Y-m", "Y/m"
    ),
    exact = FALSE,
    quiet = TRUE
  ))
  as.Date(parsed)
}

coef_dt <- function(model) {
  ct <- data.table::as.data.table(fixest::coeftable(model), keep.rownames = "term")
  data.table::setnames(
    ct,
    old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "z value", "Pr(>|z|)"),
    new = c("estimate", "std_error", "statistic", "p_value", "statistic", "p_value"),
    skip_absent = TRUE
  )
  ct[]
}

read_boundary_polygon <- function(file, target_crs = TARGET_CRS) {
  x <- suppressWarnings(sf::st_read(file, quiet = TRUE))
  x <- sf::st_make_valid(x)
  if (is.na(sf::st_crs(x))) {
    stop("Boundary file has no CRS: ", file)
  }
  x <- sf::st_transform(x, target_crs)
  geom <- sf::st_union(sf::st_geometry(x))
  geom <- sf::st_make_valid(geom)
  sf::st_sf(boundary_id = basename(file), geometry = geom)
}

empty_polygon_sf <- function(crs = TARGET_CRS, boundary_id = "empty") {
  # Keep a one-row POLYGON EMPTY rather than a zero-row sf object. This makes
  # downstream area and intersection diagnostics return 0 cleanly.
  empty_geom <- sf::st_sfc(sf::st_polygon(list()), crs = crs)
  sf::st_sf(boundary_id = boundary_id, geometry = empty_geom)
}

polygon_difference <- function(x, y, boundary_id = "difference") {
  # x and y are expected to be single-feature cumulative policy polygons in the
  # same projected CRS. A nested expansion can imply an EMPTY difference (for
  # example, 2019 ULEZ minus the cumulative 2021 ULEZ). Some sf/GEOS versions
  # represent that result as an sfc of length zero, which cannot be combined
  # with a length-one attribute column in st_sf(). Handle it explicitly.
  if (nrow(x) != 1L || nrow(y) != 1L) {
    stop("polygon_difference() expects one feature in x and one feature in y.")
  }

  if (!identical(sf::st_crs(x), sf::st_crs(y))) {
    stop("polygon_difference() requires x and y to have the same CRS.")
  }

  g <- suppressWarnings(
    sf::st_difference(sf::st_geometry(x), sf::st_geometry(y))
  )

  if (length(g) == 0L) {
    return(empty_polygon_sf(sf::st_crs(x), boundary_id = boundary_id))
  }

  g <- sf::st_make_valid(g)

  if (length(g) == 0L || all(sf::st_is_empty(g))) {
    return(empty_polygon_sf(sf::st_crs(x), boundary_id = boundary_id))
  }

  # Collapse multiple pieces to one feature so all subsequent diagnostics have
  # a stable one-row polygon object.
  g <- suppressWarnings(sf::st_union(g))

  if (length(g) == 0L || all(sf::st_is_empty(g))) {
    return(empty_polygon_sf(sf::st_crs(x), boundary_id = boundary_id))
  }

  sf::st_sf(boundary_id = boundary_id, geometry = g)
}

polygon_area_km2 <- function(x) {
  if (nrow(x) == 0L || length(sf::st_geometry(x)) == 0L ||
      all(sf::st_is_empty(x))) {
    return(0)
  }
  sum(as.numeric(sf::st_area(x)), na.rm = TRUE) / 1e6
}

compute_buffer_share <- function(receiver_sf, polygon_sf, radius_km, share_name) {
  radius_m <- radius_km * 1000

  buffers <- sf::st_buffer(receiver_sf[, "station_id", drop = FALSE], dist = radius_m)
  buffers$buffer_area_m2 <- as.numeric(sf::st_area(buffers))

  if (nrow(polygon_sf) == 0 || all(sf::st_is_empty(polygon_sf))) {
    out <- data.table(
      station_id = receiver_sf$station_id,
      share = 0
    )
    data.table::setnames(out, "share", share_name)
    return(out[])
  }

  inter <- suppressWarnings(
    sf::st_intersection(
      buffers[, "station_id", drop = FALSE],
      polygon_sf[, "boundary_id", drop = FALSE]
    )
  )

  if (nrow(inter) == 0) {
    out <- data.table(
      station_id = receiver_sf$station_id,
      share = 0
    )
    data.table::setnames(out, "share", share_name)
    return(out[])
  }

  inter$intersection_area_m2 <- as.numeric(sf::st_area(inter))
  int_dt <- data.table::as.data.table(sf::st_drop_geometry(inter))[
    , .(intersection_area_m2 = sum(intersection_area_m2, na.rm = TRUE)),
    by = station_id
  ]

  den_dt <- data.table(
    station_id = buffers$station_id,
    buffer_area_m2 = buffers$buffer_area_m2
  )

  out <- merge(den_dt, int_dt, by = "station_id", all.x = TRUE)
  out[is.na(intersection_area_m2), intersection_area_m2 := 0]
  out[, share := pmax(0, pmin(1, intersection_area_m2 / buffer_area_m2))]
  out[, c("buffer_area_m2", "intersection_area_m2") := NULL]
  data.table::setnames(out, "share", share_name)
  out[]
}

build_W_cutoff <- function(stations_sf, radius_km) {
  station_order <- stations_sf$station_id
  d <- units::drop_units(sf::st_distance(stations_sf))
  W_bin <- (d > 0) & (d <= radius_km * 1000)

  W <- matrix(
    0,
    nrow = nrow(W_bin),
    ncol = ncol(W_bin),
    dimnames = list(station_order, station_order)
  )

  rs <- rowSums(W_bin)
  if (any(rs > 0)) {
    W[rs > 0, ] <- sweep(
      W_bin[rs > 0, , drop = FALSE] * 1,
      1,
      rs[rs > 0],
      "/"
    )
  }
  W
}

component_labels <- data.table(
  term = c(
    "post_new_unexposed_prior",
    "post_new_exposed_prior",
    "post_new_S_enter",
    "post_already_S_enter",
    "post_untreated_S_enter"
  ),
  component = c(
    "Direct treatment: entrants not previously exposed",
    "Direct treatment: entrants previously exposed",
    "Entering-cohort exposure: newly treated receivers",
    "Entering-cohort exposure: already-treated receivers",
    "Entering-cohort exposure: currently untreated receivers"
  )
)


dynamic_terms_13 <- component_labels$term

semantic_support_n_13 <- function(dt, term_name) {
  switch(
    term_name,
    post_new_unexposed_prior = data.table::uniqueN(
      dt[treated_entrant == 1L & G_prior == 0L, station_id]
    ),
    post_new_exposed_prior = data.table::uniqueN(
      dt[treated_entrant == 1L & G_prior == 1L, station_id]
    ),
    post_new_S_enter = data.table::uniqueN(
      dt[treated_entrant == 1L & S_enter > EXPOSURE_EPS, station_id]
    ),
    post_already_S_enter = data.table::uniqueN(
      dt[already_treated_receiver == 1L & S_enter > EXPOSURE_EPS, station_id]
    ),
    post_untreated_S_enter = data.table::uniqueN(
      dt[untreated_receiver == 1L & S_enter > EXPOSURE_EPS, station_id]
    ),
    NA_integer_
  )
}

effective_support_one_design_13 <- function(dt, term_name, fe_rhs,
                                            mapping_name, radius_km,
                                            stack_year = NA_integer_,
                                            stack_label = NA_character_,
                                            pooled_model = NULL,
                                            residual_tol = 1e-12) {
  work <- data.table::copy(dt)
  work[, .target_x_13 := get(term_name)]

  raw_nonzero <- work[
    abs(.target_x_13) > residual_tol,
    data.table::uniqueN(station_id)
  ]
  raw_support <- semantic_support_n_13(work, term_name)

  nuisance <- setdiff(dynamic_terms_13, term_name)
  rhs <- paste(nuisance, collapse = " + ")

  fml_x <- stats::as.formula(
    paste0(".target_x_13 ~ ", rhs, " | ", fe_rhs)
  )
  fml_y <- stats::as.formula(
    paste0("y ~ ", rhs, " | ", fe_rhs)
  )

  mod_x <- tryCatch(
    fixest::feols(fml_x, data = work, fixef.rm = "none"),
    error = function(e) NULL
  )
  mod_y <- tryCatch(
    fixest::feols(fml_y, data = work, fixef.rm = "none"),
    error = function(e) NULL
  )

  component_name <- component_labels[term == term_name, component][1]
  note <- paste(
    "n_stations_with_nonzero_residualized_target is not raw support.",
    "FWL residualization can assign nonzero target residuals to comparison stations.",
    "Interpret n_eff and concentration shares jointly with raw support."
  )

  empty_row <- function(status) {
    data.table::data.table(
      mapping = mapping_name,
      radius_km = radius_km,
      stack_year = stack_year,
      stack_label = stack_label,
      term = term_name,
      component = component_name,
      raw_support_n = raw_support,
      raw_nonzero_regressor_stations = raw_nonzero,
      n_stations_with_nonzero_residualized_target = 0L,
      residual_ss_total = NA_real_,
      n_eff = NA_real_,
      max_share = NA_real_,
      top2_share = NA_real_,
      top3_share = NA_real_,
      max_share_station_id = NA_character_,
      beta_fwl = NA_real_,
      pooled_estimate = NA_real_,
      fwl_minus_pooled = NA_real_,
      residualization_status = status,
      support_interpretation_note = note
    )
  }

  if (is.null(mod_x) || is.null(mod_y)) {
    return(empty_row("auxiliary regression failed"))
  }

  xt <- as.numeric(stats::resid(mod_x))
  yt <- as.numeric(stats::resid(mod_y))
  if (length(xt) != nrow(work) || length(yt) != nrow(work)) {
    return(empty_row("residual length mismatch"))
  }

  work[, `:=`(.x_tilde_13 = xt, .y_tilde_13 = yt)]
  station_q <- work[, .(
    residual_ss = sum(.x_tilde_13^2, na.rm = TRUE),
    residual_xy = sum(.x_tilde_13 * .y_tilde_13, na.rm = TRUE)
  ), by = station_id]

  total_ss <- sum(station_q$residual_ss, na.rm = TRUE)
  total_xy <- sum(station_q$residual_xy, na.rm = TRUE)

  if (!is.finite(total_ss) || total_ss <= residual_tol) {
    return(empty_row("no residual identifying variation"))
  }

  station_q[, q_share := residual_ss / total_ss]
  station_q <- station_q[is.finite(q_share) & q_share > residual_tol]
  data.table::setorder(station_q, -q_share, station_id)

  q <- station_q$q_share
  n_eff <- 1 / sum(q^2)
  beta_fwl <- total_xy / total_ss

  pooled_est <- NA_real_
  if (!is.null(pooled_model) && term_name %in% names(stats::coef(pooled_model))) {
    pooled_est <- unname(stats::coef(pooled_model)[term_name])
  }

  data.table::data.table(
    mapping = mapping_name,
    radius_km = radius_km,
    stack_year = stack_year,
    stack_label = stack_label,
    term = term_name,
    component = component_name,
    raw_support_n = raw_support,
    raw_nonzero_regressor_stations = raw_nonzero,
    n_stations_with_nonzero_residualized_target = nrow(station_q),
    residual_ss_total = total_ss,
    n_eff = n_eff,
    max_share = if (length(q)) q[1] else NA_real_,
    top2_share = if (length(q)) sum(utils::head(q, 2)) else NA_real_,
    top3_share = if (length(q)) sum(utils::head(q, 3)) else NA_real_,
    max_share_station_id = if (nrow(station_q)) station_q$station_id[1] else NA_character_,
    beta_fwl = beta_fwl,
    pooled_estimate = pooled_est,
    fwl_minus_pooled = beta_fwl - pooled_est,
    residualization_status = "ok",
    support_interpretation_note = note
  )
}

pooled_effective_support_13 <- function(stack_dt, model, mapping_name, radius_km) {
  data.table::rbindlist(
    lapply(dynamic_terms_13, function(tt) {
      effective_support_one_design_13(
        dt = stack_dt,
        term_name = tt,
        fe_rhs = "station_stack + month_stack",
        mapping_name = mapping_name,
        radius_km = radius_km,
        pooled_model = model
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
}

cohort_effective_support_13 <- function(stack_dt, mapping_name, radius_km) {
  keys <- unique(stack_dt[, .(stack_year, stack_label)])
  data.table::rbindlist(
    lapply(seq_len(nrow(keys)), function(k) {
      sy <- keys$stack_year[k]
      sl <- keys$stack_label[k]
      dd <- stack_dt[stack_year == sy]
      data.table::rbindlist(
        lapply(dynamic_terms_13, function(tt) {
          effective_support_one_design_13(
            dt = dd,
            term_name = tt,
            fe_rhs = "station_id + m_id",
            mapping_name = mapping_name,
            radius_km = radius_km,
            stack_year = sy,
            stack_label = sl
          )
        }),
        use.names = TRUE,
        fill = TRUE
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 3. Read station assignment and outcome data
# ------------------------------------------------------------

assign <- data.table::fread(ASSIGN_FILE)
data.table::setnames(assign, clean_names_simple(names(assign)))

if (!"code" %in% names(assign)) stop("Assignment file must contain code.")
assign[, station_id := toupper(trimws(as.character(code)))]
assign[, latitude := as.numeric(latitude)]
assign[, longitude := as.numeric(longitude)]
assign[, inside_ulez_2019 := to01(inside_ulez_2019)]
assign[, inside_ulez_2021 := to01(inside_ulez_2021)]
assign[, inside_ulez_2023 := to01(inside_ulez_2023)]

if (!"opening_date" %in% names(assign)) {
  stop(
    "Assignment file must contain opening_date for the predetermined-network robustness diagnostic."
  )
}
assign[, opening_date_raw := as.character(opening_date)]
assign[, opening_date_parsed := parse_monitor_date_13(opening_date_raw)]
assign[, pre_policy_monitor := as.integer(
  !is.na(opening_date_parsed) & opening_date_parsed < POLICY_FORMATION_CUTOFF
)]

if (anyDuplicated(assign$station_id)) {
  stop("Assignment station_id is not unique. Audit before continuing.")
}
if (any(!is.finite(assign$latitude)) || any(!is.finite(assign$longitude))) {
  stop("Assignment file contains missing/non-finite coordinates.")
}

assign[, cohort_date := as.Date(NA)]
assign[inside_ulez_2019 == 1L, cohort_date := as.Date("2019-04-01")]
assign[is.na(cohort_date) & inside_ulez_2021 == 1L,
       cohort_date := as.Date("2021-10-01")]
assign[is.na(cohort_date) & inside_ulez_2023 == 1L,
       cohort_date := as.Date("2023-08-01")]
assign[, cohort_id := month_id(cohort_date)]

assign_sf <- sf::st_as_sf(
  assign,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
assign_sf <- sf::st_transform(assign_sf, TARGET_CRS)

raw_panel <- data.table::fread(RAW_PANEL_FILE)
data.table::setnames(raw_panel, clean_names_simple(names(raw_panel)))
raw_panel[, station_id := toupper(trimws(as.character(code)))]
raw_panel[, month_date := as.Date(month)]
raw_panel[, m_id := month_id(month_date)]
raw_panel[, y := log(as.numeric(no2))]
raw_panel <- raw_panel[is.finite(y)]

receiver_ids <- sort(unique(raw_panel$station_id))
receiver_meta <- assign_sf[assign_sf$station_id %in% receiver_ids, ]
receiver_meta <- receiver_meta[match(receiver_ids, receiver_meta$station_id), ]
if (!identical(receiver_meta$station_id, receiver_ids)) {
  stop("Could not align receiver stations to assignment metadata.")
}

assign[, cohort_group_13 := data.table::fcase(
  cohort_date == as.Date("2019-04-01"), "2019 entrant",
  cohort_date == as.Date("2021-10-01"), "2021 entrant",
  cohort_date == as.Date("2023-08-01"), "2023 entrant",
  is.na(cohort_date), "Never treated",
  default = "Other"
)]

predetermined_receiver_ids <- receiver_meta$station_id[
  receiver_meta$pre_policy_monitor == 1L
]

network_composition_all <- assign[, .(
  scope = "assignment source universe",
  n_stations = .N,
  n_opening_date_known = sum(!is.na(opening_date_parsed)),
  n_pre_policy = sum(pre_policy_monitor == 1L),
  pre_policy_share = mean(pre_policy_monitor == 1L)
), by = .(cohort_group = cohort_group_13)]

receiver_assign <- assign[station_id %in% receiver_ids]
network_composition_receivers <- receiver_assign[, .(
  scope = "outcome receiver universe",
  n_stations = .N,
  n_opening_date_known = sum(!is.na(opening_date_parsed)),
  n_pre_policy = sum(pre_policy_monitor == 1L),
  pre_policy_share = mean(pre_policy_monitor == 1L)
), by = .(cohort_group = cohort_group_13)]

network_composition <- data.table::rbindlist(
  list(network_composition_all, network_composition_receivers),
  use.names = TRUE, fill = TRUE
)
network_composition[, `:=`(
  policy_formation_cutoff = POLICY_FORMATION_CUTOFF,
  rule = "opening_date < policy_formation_cutoff; unknown opening dates excluded"
)]

data.table::fwrite(
  network_composition,
  file.path(TAB_DIR, "table_13_16_predetermined_monitor_network_composition.csv")
)

message("Assignment station universe: ", nrow(assign_sf))
message("Outcome receiver universe: ", length(receiver_ids))
message(
  "Predetermined monitors (assignment / outcome receivers): ",
  sum(assign$pre_policy_monitor == 1L), " / ", length(predetermined_receiver_ids),
  "; cutoff = ", as.character(POLICY_FORMATION_CUTOFF)
)
message(
  "Unknown opening dates (assignment / outcome receivers): ",
  sum(is.na(assign$opening_date_parsed)), " / ",
  sum(is.na(receiver_assign$opening_date_parsed))
)

# ------------------------------------------------------------
# 4. Read and audit ULEZ polygons
# ------------------------------------------------------------

poly_2019 <- read_boundary_polygon(BOUNDARY_FILES[["2019"]])
poly_2021 <- read_boundary_polygon(BOUNDARY_FILES[["2021"]])
poly_2023 <- read_boundary_polygon(BOUNDARY_FILES[["2023"]])

# Entering-area polygons isolate the geographic increment associated with each
# ULEZ expansion. These are the boundary-based analogues of source cohorts.
enter_2019 <- poly_2019
enter_2021 <- polygon_difference(poly_2021, poly_2019)
enter_2023 <- polygon_difference(poly_2023, poly_2021)

# Prior-area polygons correspond to the area already regulated before cohort c.
prior_2019 <- empty_polygon_sf()
prior_2021 <- poly_2019
prior_2023 <- poly_2021

outside_2019_from_2021 <- polygon_difference(poly_2019, poly_2021)
outside_2021_from_2023 <- polygon_difference(poly_2021, poly_2023)

boundary_geometry_audit <- data.table(
  object = c(
    "ULEZ cumulative 2019",
    "ULEZ cumulative 2021",
    "ULEZ cumulative 2023",
    "Entering area 2019",
    "Entering area 2021",
    "Entering area 2023",
    "2019 area outside 2021 boundary",
    "2021 area outside 2023 boundary"
  ),
  area_km2 = c(
    polygon_area_km2(poly_2019),
    polygon_area_km2(poly_2021),
    polygon_area_km2(poly_2023),
    polygon_area_km2(enter_2019),
    polygon_area_km2(enter_2021),
    polygon_area_km2(enter_2023),
    polygon_area_km2(outside_2019_from_2021),
    polygon_area_km2(outside_2021_from_2023)
  )
)

data.table::fwrite(
  boundary_geometry_audit,
  file.path(TAB_DIR, "table_13_01_boundary_geometry_audit.csv")
)

# Validate that the point-in-polygon assignment implied by the boundary files
# agrees with the treatment-assignment file used by the existing pipeline.
assignment_validation <- data.table(station_id = assign_sf$station_id)
for (yr in c("2019", "2021", "2023")) {
  poly <- get(paste0("poly_", yr))
  inside_geo <- as.integer(lengths(sf::st_intersects(assign_sf, poly)) > 0)
  inside_file <- assign[[paste0("inside_ulez_", yr)]]

  assignment_validation[, (paste0("inside_boundary_", yr)) := inside_geo]
  assignment_validation[, (paste0("inside_assignment_", yr)) := inside_file]
  assignment_validation[, (paste0("mismatch_", yr)) := as.integer(inside_geo != inside_file)]
}

data.table::fwrite(
  assignment_validation,
  file.path(TAB_DIR, "table_13_02_boundary_assignment_validation.csv")
)

message("Boundary assignment mismatches:")
print(assignment_validation[, .(
  mismatch_2019 = sum(mismatch_2019),
  mismatch_2021 = sum(mismatch_2021),
  mismatch_2023 = sum(mismatch_2023)
)])

# ------------------------------------------------------------
# 5. Boundary-based exposure shares at receiver stations
# ------------------------------------------------------------

boundary_exposure_by_radius <- list()

for (r in RADII_KM) {
  message("Computing boundary-area shares for radius: ", r, " km")

  tmp <- data.table(station_id = receiver_meta$station_id)

  share_objects <- list(
    cum_2019   = poly_2019,
    cum_2021   = poly_2021,
    cum_2023   = poly_2023,
    enter_2019 = enter_2019,
    enter_2021 = enter_2021,
    enter_2023 = enter_2023,
    prior_2019 = prior_2019,
    prior_2021 = prior_2021,
    prior_2023 = prior_2023
  )

  for (nm in names(share_objects)) {
    sh <- compute_buffer_share(
      receiver_meta,
      share_objects[[nm]],
      radius_km = r,
      share_name = paste0("S_area_", nm)
    )
    tmp <- merge(tmp, sh, by = "station_id", all.x = TRUE, sort = FALSE)
  }

  tmp[, radius_km := r]
  boundary_exposure_by_radius[[as.character(r)]] <- tmp[]
}

boundary_exposure <- data.table::rbindlist(boundary_exposure_by_radius, use.names = TRUE)
saveRDS(boundary_exposure, file.path(RDS_DIR, "boundary_exposure_receiver_station_level.rds"))

# ------------------------------------------------------------
# 6. STATIC diagnostic: exact current benchmark vs boundary mapping
# ------------------------------------------------------------

static_panel <- data.table::fread(STATIC_PANEL_FILE)
data.table::setnames(static_panel, clean_names_simple(names(static_panel)))

required_static <- c("code", "month", "log_y", "d_active", "s_bench_0_3km")
if (!all(required_static %in% names(static_panel))) {
  stop(
    "Static benchmark panel is missing required variables: ",
    paste(setdiff(required_static, names(static_panel)), collapse = ", ")
  )
}

static_panel[, station_id := toupper(trimws(as.character(code)))]
static_panel[, month_date := as.Date(month)]
static_panel[, D := as.integer(d_active)]
static_panel[, S_station := as.numeric(s_bench_0_3km)]
static_panel[, y := as.numeric(log_y)]

for (r in RADII_KM) {
  # The current manuscript benchmark is 3 km. If more radii are later added,
  # station-based static counterparts should be supplied explicitly.
  if (r != 3) next

  bnd <- boundary_exposure[radius_km == r]
  static_r <- merge(
    static_panel,
    bnd[, .(station_id, S_area_cum_2019, S_area_cum_2021, S_area_cum_2023)],
    by = "station_id",
    all.x = TRUE,
    sort = FALSE
  )

  static_r[, S_boundary := 0]
  static_r[
    month_date >= STATIC_STAGE_START[1] & month_date < STATIC_STAGE_START[2],
    S_boundary := S_area_cum_2019
  ]
  static_r[
    month_date >= STATIC_STAGE_START[2] & month_date < STATIC_STAGE_START[3],
    S_boundary := S_area_cum_2021
  ]
  static_r[
    month_date >= STATIC_STAGE_START[3],
    S_boundary := S_area_cum_2023
  ]

  # Mapping correspondence on unique station x policy-phase cells, avoiding
  # mechanical over-weighting by the number of months in each phase.
  static_r[, phase := data.table::fcase(
    month_date < STATIC_STAGE_START[1], "pre",
    month_date < STATIC_STAGE_START[2], "2019",
    month_date < STATIC_STAGE_START[3], "2021",
    default = "2023"
  )]

  map_cells <- unique(
    static_r[phase != "pre", .(station_id, phase, S_station, S_boundary)]
  )

  static_corr <- map_cells[, .(
    radius_km = r,
    n_station_phase = .N,
    correlation = safe_cor(S_station, S_boundary),
    mean_station_exposure = safe_mean(S_station),
    mean_boundary_exposure = safe_mean(S_boundary),
    mean_abs_difference = safe_mean(abs(S_station - S_boundary)),
    n_positive_station = sum(S_station > EXPOSURE_EPS),
    n_positive_boundary = sum(S_boundary > EXPOSURE_EPS),
    binary_exposure_agreement = mean(
      (S_station > EXPOSURE_EPS) == (S_boundary > EXPOSURE_EPS)
    )
  ), by = phase]

  data.table::fwrite(
    static_corr,
    file.path(TAB_DIR, "table_13_03_static_mapping_correspondence.csv")
  )

  # Support under each mapping.
  support_list <- list()
  for (mapping_name in c("station_3km", "boundary_area_share_3km")) {
    S_var <- if (mapping_name == "station_3km") "S_station" else "S_boundary"
    dt <- data.table::copy(static_r)
    dt[, S_tmp := get(S_var)]
    dt[, G := as.integer(S_tmp > EXPOSURE_EPS)]
    dt[, cell := paste0("D", D, "/S", G)]

    support_list[[mapping_name]] <- dt[, .(
      radius_km = r,
      n_obs = .N,
      n_stations = data.table::uniqueN(station_id),
      mean_exposure = safe_mean(S_tmp),
      mean_positive_exposure = safe_mean(S_tmp[S_tmp > EXPOSURE_EPS])
    ), by = .(cell)][, mapping := mapping_name]
  }
  static_support <- data.table::rbindlist(support_list, use.names = TRUE, fill = TRUE)
  data.table::setcolorder(static_support, c("mapping", "radius_km", "cell"))
  data.table::fwrite(
    static_support,
    file.path(TAB_DIR, "table_13_04_static_support_station_vs_boundary.csv")
  )

  fit_static_mapping <- function(dt, S_var, mapping_name) {
    d <- data.table::copy(dt)
    d[, S := get(S_var)]
    d[, DS := D * S]
    d[, US := (1 - D) * S]

    mod <- fixest::feols(
      y ~ D + DS + US | station_id + month_date,
      cluster = ~ station_id,
      data = d
    )

    ct <- coef_dt(mod)
    labels <- data.table(
      term = c("D", "DS", "US"),
      component = c(
        "Direct treatment",
        "Treated-exposure slope",
        "Currently untreated exposure slope"
      )
    )
    out <- merge(labels, ct, by = "term", all.x = TRUE, sort = FALSE)
    out[, mapping := mapping_name]
    out[, radius_km := r]
    out[, percent_per_unit_exposure := pct_effect(estimate)]

    scales <- data.table(
      term = c("D", "DS", "US"),
      exposure_scale = c(
        1,
        safe_mean(d[D == 1 & S > EXPOSURE_EPS, S]),
        safe_mean(d[D == 0 & S > EXPOSURE_EPS, S])
      )
    )
    implied <- merge(out, scales, by = "term", all.x = TRUE, sort = FALSE)
    implied[, implied_percent := 100 * (exp(estimate * exposure_scale) - 1)]
    implied[, ci_low_percent := 100 * (
      exp((estimate - 1.96 * std_error) * exposure_scale) - 1
    )]
    implied[, ci_high_percent := 100 * (
      exp((estimate + 1.96 * std_error) * exposure_scale) - 1
    )]

    # For the boundary-area-share mapping, D=1/S=0 is structurally absent:
    # any treated station lies inside the policy polygon and therefore has a
    # positive treated-area share in its 3 km buffer. The standalone coefficient
    # on D is consequently an intercept at S=0 outside observed treated support.
    # Report the total treated effect at the observed mean S among treated units:
    # beta_D + s_bar * beta_DS, with covariance-aware uncertainty.
    b <- stats::coef(mod)
    V <- stats::vcov(mod)
    s_treated <- safe_mean(d[D == 1L, S])

    if (all(c("D", "DS") %in% names(b)) && is.finite(s_treated)) {
      est_total <- unname(b["D"] + s_treated * b["DS"])
      var_total <- unname(
        V["D", "D"] +
          s_treated^2 * V["DS", "DS"] +
          2 * s_treated * V["D", "DS"]
      )
      se_total <- sqrt(pmax(var_total, 0))
    } else {
      est_total <- NA_real_
      se_total <- NA_real_
    }

    total_treated <- data.table::data.table(
      mapping = mapping_name,
      radius_km = r,
      estimand = "Total treated effect at observed mean treated exposure",
      n_treated_stations = data.table::uniqueN(d[D == 1L, station_id]),
      mean_observed_treated_exposure = s_treated,
      log_effect = est_total,
      std_error = se_total,
      p_value = if (is.finite(est_total) && is.finite(se_total) && se_total > 0) {
        2 * stats::pnorm(-abs(est_total / se_total))
      } else NA_real_,
      implied_percent = if (is.finite(est_total)) 100 * (exp(est_total) - 1) else NA_real_,
      ci_low_percent = if (is.finite(est_total) && is.finite(se_total)) {
        100 * (exp(est_total - 1.96 * se_total) - 1)
      } else NA_real_,
      ci_high_percent = if (is.finite(est_total) && is.finite(se_total)) {
        100 * (exp(est_total + 1.96 * se_total) - 1)
      } else NA_real_
    )

    list(
      model = mod,
      coefficients = out[],
      implied = implied[],
      total_treated_observed_support = total_treated[]
    )
  }

  st_fit <- fit_static_mapping(static_r, "S_station", "station_3km")
  bd_fit <- fit_static_mapping(static_r, "S_boundary", "boundary_area_share_3km")

  static_coefs <- data.table::rbindlist(
    list(st_fit$coefficients, bd_fit$coefficients),
    use.names = TRUE,
    fill = TRUE
  )
  static_implied <- data.table::rbindlist(
    list(st_fit$implied, bd_fit$implied),
    use.names = TRUE,
    fill = TRUE
  )
  static_total_treated <- data.table::rbindlist(
    list(
      st_fit$total_treated_observed_support,
      bd_fit$total_treated_observed_support
    ),
    use.names = TRUE,
    fill = TRUE
  )

  data.table::fwrite(
    static_coefs,
    file.path(TAB_DIR, "table_13_05_static_coefficients_station_vs_boundary.csv")
  )
  data.table::fwrite(
    static_implied,
    file.path(TAB_DIR, "table_13_06_static_implied_effects_station_vs_boundary.csv")
  )
  data.table::fwrite(
    static_total_treated,
    file.path(TAB_DIR, "table_13_12_static_total_treated_effect_observed_support.csv")
  )

  # Scatter diagnostic at the station x policy-phase level.
  p_static <- ggplot2::ggplot(
    map_cells,
    ggplot2::aes(x = S_station, y = S_boundary)
  ) +
    ggplot2::geom_point(alpha = 0.65) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~ phase) +
    ggplot2::labs(
      x = "Station-based exposure (0--3 km)",
      y = "Boundary-area share within 3 km",
      title = "Static exposure mapping: monitoring-station versus policy-boundary representation"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  ggplot2::ggsave(
    file.path(FIG_DIR, "fig_13_01_static_station_vs_boundary_exposure.png"),
    p_static, width = 9, height = 4.8, dpi = 300
  )
  ggplot2::ggsave(
    file.path(FIG_DIR, "fig_13_01_static_station_vs_boundary_exposure.pdf"),
    p_static, width = 9, height = 4.8
  )

  saveRDS(st_fit$model, file.path(RDS_DIR, "static_model_station_3km.rds"))
  saveRDS(bd_fit$model, file.path(RDS_DIR, "static_model_boundary_3km.rds"))

  # ----------------------------------------------------------
  # 6.1 Predetermined monitoring-network robustness (static)
  # ----------------------------------------------------------
  # Measurement support is a maintained data restriction: restricting receivers
  # to monitors opened before policy formation changes the population over which
  # outcomes are observed. The station-based mapping has a second channel because
  # monitor locations also define exposure sources. We therefore distinguish:
  #   A. pre-policy receivers, original station-source universe;
  #   B. pre-policy receivers AND pre-policy station-source universe;
  #   C. pre-policy receivers, boundary-based exposure (no monitor sources).
  # This is not a claim that pre-policy monitors are spatially random.

  static_pre <- data.table::copy(
    static_r[station_id %in% predetermined_receiver_ids]
  )

  # Rebuild station exposure using only predetermined outcome monitors as both
  # receivers and sources, matching the historical static benchmark convention
  # within the restricted network.
  receiver_meta_pre <- receiver_meta[
    receiver_meta$station_id %in% predetermined_receiver_ids,
  ]
  receiver_meta_pre <- receiver_meta_pre[order(receiver_meta_pre$station_id), ]

  if (nrow(receiver_meta_pre) >= 2L) {
    W_static_pre <- build_W_cutoff(receiver_meta_pre, r)
    pre_order <- rownames(W_static_pre)
    pre_assign <- assign[match(pre_order, station_id)]

    static_pre_source_exposure <- data.table::data.table(station_id = pre_order)
    for (yy in c(2019L, 2021L, 2023L)) {
      z <- as.numeric(pre_assign[[paste0("inside_ulez_", yy)]] == 1L)
      static_pre_source_exposure[, (paste0("S_pre_sources_", yy)) :=
                                   as.numeric(W_static_pre %*% z)]
    }

    static_pre <- merge(
      static_pre,
      static_pre_source_exposure,
      by = "station_id", all.x = TRUE, sort = FALSE
    )
    static_pre[, S_station_pre_sources := 0]
    static_pre[
      month_date >= STATIC_STAGE_START[1] & month_date < STATIC_STAGE_START[2],
      S_station_pre_sources := S_pre_sources_2019
    ]
    static_pre[
      month_date >= STATIC_STAGE_START[2] & month_date < STATIC_STAGE_START[3],
      S_station_pre_sources := S_pre_sources_2021
    ]
    static_pre[
      month_date >= STATIC_STAGE_START[3],
      S_station_pre_sources := S_pre_sources_2023
    ]
  } else {
    static_pre[, S_station_pre_sources := NA_real_]
  }

  static_pred_variants <- list(
    station_full_sample = list(
      dt = static_r, S_var = "S_station", mapping = "station_full_sample"
    ),
    boundary_full_sample = list(
      dt = static_r, S_var = "S_boundary", mapping = "boundary_full_sample"
    ),
    station_pre2013_receivers_original_sources = list(
      dt = static_pre, S_var = "S_station",
      mapping = "station_pre2013_receivers_original_sources"
    ),
    station_pre2013_receivers_pre2013_sources = list(
      dt = static_pre, S_var = "S_station_pre_sources",
      mapping = "station_pre2013_receivers_pre2013_sources"
    ),
    boundary_pre2013_receivers = list(
      dt = static_pre, S_var = "S_boundary",
      mapping = "boundary_pre2013_receivers"
    )
  )

  static_pred_support <- data.table::rbindlist(lapply(static_pred_variants, function(v) {
    dd <- data.table::copy(v$dt)
    dd[, S_tmp := get(v$S_var)]
    dd <- dd[is.finite(S_tmp)]
    dd[, G := as.integer(S_tmp > EXPOSURE_EPS)]
    dd[, cell := paste0("D", D, "/S", G)]
    dd[, .(
      policy_formation_cutoff = POLICY_FORMATION_CUTOFF,
      n_obs = .N,
      n_stations = data.table::uniqueN(station_id),
      mean_exposure = safe_mean(S_tmp),
      mean_positive_exposure = safe_mean(S_tmp[S_tmp > EXPOSURE_EPS])
    ), by = cell][, variant := v$mapping]
  }), use.names = TRUE, fill = TRUE)
  data.table::setcolorder(
    static_pred_support,
    c("variant", "policy_formation_cutoff", "cell")
  )

  static_pred_fits <- lapply(static_pred_variants, function(v) {
    dd <- data.table::copy(v$dt)
    dd <- dd[is.finite(get(v$S_var))]
    if (nrow(dd) == 0L || data.table::uniqueN(dd$station_id) < 2L) return(NULL)
    tryCatch(
      fit_static_mapping(dd, v$S_var, v$mapping),
      error = function(e) {
        message("Static predetermined-network fit failed for ", v$mapping, ": ", e$message)
        NULL
      }
    )
  })

  static_pred_coefs <- data.table::rbindlist(lapply(static_pred_fits, function(ff) {
    if (is.null(ff)) return(NULL)
    x <- data.table::copy(ff$coefficients)
    x[, policy_formation_cutoff := POLICY_FORMATION_CUTOFF]
    x[]
  }), use.names = TRUE, fill = TRUE)

  static_pred_total <- data.table::rbindlist(lapply(static_pred_fits, function(ff) {
    if (is.null(ff)) return(NULL)
    x <- data.table::copy(ff$total_treated_observed_support)
    x[, policy_formation_cutoff := POLICY_FORMATION_CUTOFF]
    x[]
  }), use.names = TRUE, fill = TRUE)

  data.table::fwrite(
    static_pred_support,
    file.path(TAB_DIR, "table_13_17_static_predetermined_network_support.csv")
  )
  data.table::fwrite(
    static_pred_coefs,
    file.path(TAB_DIR, "table_13_18_static_predetermined_network_coefficients.csv")
  )
  data.table::fwrite(
    static_pred_total,
    file.path(TAB_DIR, "table_13_19_static_predetermined_network_total_treated_effects.csv")
  )
}

# ------------------------------------------------------------
# 7. DYNAMIC source-specific exposure: station vs boundary
# ------------------------------------------------------------

# Merge receiver cohort information into the outcome panel.
receiver_cohort <- assign[, .(station_id, cohort_date, cohort_id)]
dyn_panel <- merge(
  raw_panel,
  receiver_cohort,
  by = "station_id",
  all.x = TRUE,
  sort = FALSE
)

# Full assignment-universe W for the station-based dynamic benchmark, following
# script 12-v4. Sources therefore do not depend on whether they have outcome data.
for (r in RADII_KM) {
  message("Building dynamic source-specific exposures for radius: ", r, " km")

  W_full <- build_W_cutoff(assign_sf, r)
  station_order <- rownames(W_full)
  cohort_by_station <- assign$cohort_date[match(station_order, assign$station_id)]

  station_source_list <- list()
  boundary_source_list <- list()

  bnd <- boundary_exposure[radius_km == r]

  for (k in seq_along(EXPANSION_DATES)) {
    c_date <- EXPANSION_DATES[k]
    c_year <- EXPANSION_YEARS[k]
    c_label <- EXPANSION_LABELS[k]

    z_enter <- as.numeric(!is.na(cohort_by_station) & cohort_by_station == c_date)
    z_prior <- as.numeric(!is.na(cohort_by_station) & cohort_by_station < c_date)

    S_enter_station <- as.numeric(W_full %*% z_enter)
    S_prior_station <- as.numeric(W_full %*% z_prior)

    st_src <- data.table(
      station_id = station_order,
      stack_date = c_date,
      stack_year = c_year,
      stack_label = c_label,
      S_enter = S_enter_station,
      S_prior = S_prior_station
    )
    st_src <- st_src[station_id %in% receiver_ids]
    st_src[, mapping := paste0("station_", r, "km")]
    station_source_list[[k]] <- st_src[]

    enter_col <- paste0("S_area_enter_", c_year)
    prior_col <- paste0("S_area_prior_", c_year)

    bd_src <- bnd[, .(
      station_id,
      stack_date = c_date,
      stack_year = c_year,
      stack_label = c_label,
      S_enter = get(enter_col),
      S_prior = get(prior_col)
    )]
    bd_src[, mapping := paste0("boundary_area_share_", r, "km")]
    boundary_source_list[[k]] <- bd_src[]
  }

  station_sources <- data.table::rbindlist(station_source_list)
  boundary_sources <- data.table::rbindlist(boundary_source_list)

  # Receiver-status label for mapping correspondence. This status is relative to
  # the entering cohort and is distinct from calendar-time risk-set eligibility.
  receiver_status_for_stack <- function(station_ids, stack_date) {
    cd <- assign$cohort_date[match(station_ids, assign$station_id)]
    data.table::fcase(
      !is.na(cd) & cd == stack_date, "newly_treated",
      !is.na(cd) & cd < stack_date, "already_treated",
      !is.na(cd) & cd > stack_date, "future_treated",
      is.na(cd), "never_treated",
      default = "unclassified"
    )
  }

  corr_station_level <- merge(
    station_sources[, .(
      station_id, stack_date, stack_year, stack_label,
      S_enter_station = S_enter,
      S_prior_station = S_prior
    )],
    boundary_sources[, .(
      station_id, stack_date,
      S_enter_boundary = S_enter,
      S_prior_boundary = S_prior
    )],
    by = c("station_id", "stack_date"),
    all = TRUE
  )

  corr_station_level[, receiver_status := receiver_status_for_stack(station_id, stack_date)]
  corr_station_level[, radius_km := r]

  dynamic_corr <- corr_station_level[, .(
    n_stations = .N,
    correlation_enter = safe_cor(S_enter_station, S_enter_boundary),
    correlation_prior = safe_cor(S_prior_station, S_prior_boundary),
    mean_enter_station = safe_mean(S_enter_station),
    mean_enter_boundary = safe_mean(S_enter_boundary),
    mean_abs_diff_enter = safe_mean(abs(S_enter_station - S_enter_boundary)),
    n_positive_enter_station = sum(S_enter_station > EXPOSURE_EPS),
    n_positive_enter_boundary = sum(S_enter_boundary > EXPOSURE_EPS),
    binary_enter_agreement = mean(
      (S_enter_station > EXPOSURE_EPS) == (S_enter_boundary > EXPOSURE_EPS)
    )
  ), by = .(radius_km, stack_year, stack_label, receiver_status)]

  data.table::fwrite(
    dynamic_corr,
    file.path(TAB_DIR, "table_13_07_dynamic_mapping_correspondence_by_stack_status.csv")
  )
  data.table::fwrite(
    corr_station_level,
    file.path(TAB_DIR, "table_13_11_receiver_exposure_station_level.csv")
  )

  # ----------------------------------------------------------
  # 7.1 Build calendar-time-correct stacks
  # ----------------------------------------------------------

  build_stack <- function(panel_dt, source_dt, mapping_name) {
    out <- list()

    for (k in seq_along(EXPANSION_DATES)) {
      c_date <- EXPANSION_DATES[k]
      c_id <- month_id(c_date)
      c_year <- EXPANSION_YEARS[k]
      c_label <- EXPANSION_LABELS[k]

      dt <- data.table::copy(
        panel_dt[m_id >= c_id - L_PRE & m_id <= c_id + H_POST]
      )

      src <- source_dt[stack_date == c_date, .(station_id, S_enter, S_prior)]
      dt <- merge(dt, src, by = "station_id", all.x = TRUE, sort = FALSE)
      dt[is.na(S_enter), S_enter := 0]
      dt[is.na(S_prior), S_prior := 0]

      dt[, stack_date := c_date]
      dt[, stack_year := c_year]
      dt[, stack_label := c_label]
      dt[, stack_id := as.character(c_year)]
      dt[, mapping := mapping_name]
      dt[, c_id := c_id]
      dt[, rel_k := m_id - c_id]
      dt[, post := as.integer(m_id >= c_id)]

      dt[, new_cohort := as.integer(!is.na(cohort_id) & cohort_id == c_id)]
      dt[, future_cohort := as.integer(!is.na(cohort_id) & cohort_id > c_id)]
      dt[, prior_cohort := as.integer(!is.na(cohort_id) & cohort_id < c_id)]
      dt[, never_cohort := as.integer(is.na(cohort_id))]

      # Calendar-time risk-set rule from script 12-v4.
      dt[, future_at_risk := as.integer(
        future_cohort == 1L & !is.na(cohort_id) & m_id < cohort_id
      )]
      dt[, future_after_own_treatment := as.integer(
        future_cohort == 1L & !is.na(cohort_id) & m_id >= cohort_id
      )]
      dt[, risk_set_eligible := as.integer(future_after_own_treatment == 0L)]
      dt[, currently_untreated_receiver := as.integer(
        future_at_risk == 1L | never_cohort == 1L
      )]

      dt[, treated_entrant := new_cohort]
      dt[, already_treated_receiver := prior_cohort]
      dt[, untreated_receiver := currently_untreated_receiver]
      dt[, G_prior := as.integer(S_prior > EXPOSURE_EPS)]

      dt[, post_new_unexposed_prior := post * treated_entrant * (1L - G_prior)]
      dt[, post_new_exposed_prior := post * treated_entrant * G_prior]
      dt[, post_new_S_enter := post * treated_entrant * S_enter]
      dt[, post_already_S_enter := post * already_treated_receiver * S_enter]
      dt[, post_untreated_S_enter := post * untreated_receiver * S_enter]

      # Censor future-treated observations after own treatment.
      dt <- dt[risk_set_eligible == 1L]

      dt[, station_stack := paste0(stack_id, "::", station_id)]
      dt[, month_stack := paste0(stack_id, "::", m_id)]

      out[[k]] <- dt[]
    }

    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }

  stack_station <- build_stack(
    dyn_panel,
    station_sources,
    paste0("station_", r, "km")
  )
  stack_boundary <- build_stack(
    dyn_panel,
    boundary_sources,
    paste0("boundary_area_share_", r, "km")
  )

  # ----------------------------------------------------------
  # 7.2 Support comparison
  # ----------------------------------------------------------

  support_one_stack <- function(dt) {
    dt[, .(
      direct_new_unexposed_n = data.table::uniqueN(
        station_id[treated_entrant == 1L & G_prior == 0L]
      ),
      direct_new_exposed_n = data.table::uniqueN(
        station_id[treated_entrant == 1L & G_prior == 1L]
      ),
      entrant_exposure_n = data.table::uniqueN(
        station_id[treated_entrant == 1L & S_enter > EXPOSURE_EPS]
      ),
      already_treated_exposure_n = data.table::uniqueN(
        station_id[already_treated_receiver == 1L & S_enter > EXPOSURE_EPS]
      ),
      currently_untreated_exposure_n = data.table::uniqueN(
        station_id[untreated_receiver == 1L & S_enter > EXPOSURE_EPS]
      ),
      mean_S_enter_new = safe_mean(
        S_enter[treated_entrant == 1L & S_enter > EXPOSURE_EPS]
      ),
      mean_S_enter_already = safe_mean(
        S_enter[already_treated_receiver == 1L & S_enter > EXPOSURE_EPS]
      ),
      mean_S_enter_untreated = safe_mean(
        S_enter[untreated_receiver == 1L & S_enter > EXPOSURE_EPS]
      )
    ), by = .(mapping, stack_year, stack_label)]
  }

  dyn_support <- data.table::rbindlist(
    list(
      support_one_stack(stack_station),
      support_one_stack(stack_boundary)
    ),
    use.names = TRUE,
    fill = TRUE
  )
  dyn_support[, radius_km := r]
  data.table::setcolorder(dyn_support, c("mapping", "radius_km", "stack_year", "stack_label"))

  data.table::fwrite(
    dyn_support,
    file.path(TAB_DIR, "table_13_08_dynamic_support_station_vs_boundary.csv")
  )

  # ----------------------------------------------------------
  # 7.3 Pooled stacked models under each mapping
  # ----------------------------------------------------------

  fit_dynamic_mapping <- function(stack_dt, mapping_name) {
    mod <- fixest::feols(
      y ~ post_new_unexposed_prior +
        post_new_exposed_prior +
        post_new_S_enter +
        post_already_S_enter +
        post_untreated_S_enter |
        station_stack + month_stack,
      cluster = ~ station_id,
      data = stack_dt
    )

    ct <- coef_dt(mod)
    out <- merge(component_labels, ct, by = "term", all.x = TRUE, sort = FALSE)
    out[, mapping := mapping_name]
    out[, radius_km := r]

    scales <- data.table(
      term = component_labels$term,
      exposure_scale = c(
        1,
        1,
        safe_mean(stack_dt[treated_entrant == 1L & S_enter > EXPOSURE_EPS, S_enter]),
        safe_mean(stack_dt[already_treated_receiver == 1L & S_enter > EXPOSURE_EPS, S_enter]),
        safe_mean(stack_dt[untreated_receiver == 1L & S_enter > EXPOSURE_EPS, S_enter])
      )
    )

    implied <- merge(out, scales, by = "term", all.x = TRUE, sort = FALSE)
    implied[, implied_percent := 100 * (exp(estimate * exposure_scale) - 1)]
    implied[, ci_low_percent := 100 * (
      exp((estimate - 1.96 * std_error) * exposure_scale) - 1
    )]
    implied[, ci_high_percent := 100 * (
      exp((estimate + 1.96 * std_error) * exposure_scale) - 1
    )]

    # Direct-treatment coefficients are intercepts with respect to entering-cohort
    # exposure. Because most (and under the boundary mapping all) entrants have
    # positive S_enter, compare mappings using the total entrant effect evaluated
    # at the observed mean S_enter within each prior-exposure group.
    b <- stats::coef(mod)
    V <- stats::vcov(mod)
    eta_name <- "post_new_S_enter"

    combined_one <- function(g_prior_value, direct_name, group_label) {
      dd <- stack_dt[treated_entrant == 1L & G_prior == g_prior_value]
      sbar <- safe_mean(dd$S_enter)
      nst <- data.table::uniqueN(dd$station_id)

      if (!all(c(direct_name, eta_name) %in% names(b)) || !is.finite(sbar)) {
        return(data.table::data.table(
          mapping = mapping_name, radius_km = r, entrant_group = group_label,
          n_stations = nst, mean_observed_S_enter = sbar,
          log_effect = NA_real_, std_error = NA_real_, p_value = NA_real_,
          implied_percent = NA_real_, ci_low_percent = NA_real_, ci_high_percent = NA_real_
        ))
      }

      est <- unname(b[direct_name] + sbar * b[eta_name])
      vv <- unname(
        V[direct_name, direct_name] +
          sbar^2 * V[eta_name, eta_name] +
          2 * sbar * V[direct_name, eta_name]
      )
      se <- sqrt(pmax(vv, 0))

      data.table::data.table(
        mapping = mapping_name,
        radius_km = r,
        entrant_group = group_label,
        n_stations = nst,
        mean_observed_S_enter = sbar,
        log_effect = est,
        std_error = se,
        p_value = 2 * stats::pnorm(-abs(est / se)),
        implied_percent = 100 * (exp(est) - 1),
        ci_low_percent = 100 * (exp(est - 1.96 * se) - 1),
        ci_high_percent = 100 * (exp(est + 1.96 * se) - 1)
      )
    }

    combined_entrants <- data.table::rbindlist(
      list(
        combined_one(
          0L, "post_new_unexposed_prior",
          "Entrants not previously exposed"
        ),
        combined_one(
          1L, "post_new_exposed_prior",
          "Entrants previously exposed"
        )
      ),
      use.names = TRUE,
      fill = TRUE
    )

    list(
      model = mod,
      coefficients = out[],
      implied = implied[],
      combined_entrants = combined_entrants[]
    )
  }

  fit_station <- fit_dynamic_mapping(stack_station, paste0("station_", r, "km"))
  fit_boundary <- fit_dynamic_mapping(
    stack_boundary,
    paste0("boundary_area_share_", r, "km")
  )

  dyn_coefs <- data.table::rbindlist(
    list(fit_station$coefficients, fit_boundary$coefficients),
    use.names = TRUE,
    fill = TRUE
  )
  dyn_implied <- data.table::rbindlist(
    list(fit_station$implied, fit_boundary$implied),
    use.names = TRUE,
    fill = TRUE
  )
  dyn_combined_entrants <- data.table::rbindlist(
    list(fit_station$combined_entrants, fit_boundary$combined_entrants),
    use.names = TRUE,
    fill = TRUE
  )

  pooled_eff_station <- pooled_effective_support_13(
    stack_station, fit_station$model, paste0("station_", r, "km"), r
  )
  pooled_eff_boundary <- pooled_effective_support_13(
    stack_boundary, fit_boundary$model, paste0("boundary_area_share_", r, "km"), r
  )
  pooled_eff_13 <- data.table::rbindlist(
    list(pooled_eff_station, pooled_eff_boundary),
    use.names = TRUE,
    fill = TRUE
  )

  cohort_eff_station <- cohort_effective_support_13(
    stack_station, paste0("station_", r, "km"), r
  )
  cohort_eff_boundary <- cohort_effective_support_13(
    stack_boundary, paste0("boundary_area_share_", r, "km"), r
  )
  cohort_eff_13 <- data.table::rbindlist(
    list(cohort_eff_station, cohort_eff_boundary),
    use.names = TRUE,
    fill = TRUE
  )

  data.table::fwrite(
    dyn_coefs,
    file.path(TAB_DIR, "table_13_09_dynamic_coefficients_station_vs_boundary.csv")
  )
  data.table::fwrite(
    dyn_implied,
    file.path(TAB_DIR, "table_13_10_dynamic_implied_effects_station_vs_boundary.csv")
  )
  data.table::fwrite(
    dyn_combined_entrants,
    file.path(TAB_DIR, "table_13_13_dynamic_combined_entrant_effects_observed_support.csv")
  )
  data.table::fwrite(
    pooled_eff_13,
    file.path(TAB_DIR, "table_13_14_pooled_effective_identifying_support_station_vs_boundary.csv")
  )
  data.table::fwrite(
    cohort_eff_13,
    file.path(TAB_DIR, "table_13_15_cohort_effective_identifying_support_station_vs_boundary.csv")
  )

  # ----------------------------------------------------------
  # 7.4 Predetermined monitoring-network robustness (dynamic)
  # ----------------------------------------------------------
  # Three restricted designs isolate distinct roles of the monitoring network:
  #   1. boundary_pre2013_receivers:
  #        outcome support restricted to predetermined monitors; exposure is
  #        constructed from policy geography and does not use monitor sources.
  #   2. station_pre2013_receivers_all_sources:
  #        outcome support restricted, but exposure still uses the complete
  #        assignment-station source universe.
  #   3. station_pre2013_receivers_pre2013_sources:
  #        both receiver support and the station-based exposure-source network are
  #        restricted to monitors opened before policy formation.
  # Full-sample station and boundary estimates are retained as reference rows.

  dyn_panel_pre <- dyn_panel[station_id %in% predetermined_receiver_ids]

  assign_pre_sf <- assign_sf[assign_sf$pre_policy_monitor == 1L, ]
  station_sources_pre <- data.table::data.table()
  if (nrow(assign_pre_sf) >= 2L) {
    W_pre <- build_W_cutoff(assign_pre_sf, r)
    pre_station_order <- rownames(W_pre)
    pre_cohort_by_station <- assign$cohort_date[
      match(pre_station_order, assign$station_id)
    ]

    pre_source_list <- vector("list", length(EXPANSION_DATES))
    for (kk in seq_along(EXPANSION_DATES)) {
      cc_date <- EXPANSION_DATES[kk]
      cc_year <- EXPANSION_YEARS[kk]
      cc_label <- EXPANSION_LABELS[kk]
      zz_enter <- as.numeric(
        !is.na(pre_cohort_by_station) & pre_cohort_by_station == cc_date
      )
      zz_prior <- as.numeric(
        !is.na(pre_cohort_by_station) & pre_cohort_by_station < cc_date
      )
      pre_src <- data.table::data.table(
        station_id = pre_station_order,
        stack_date = cc_date,
        stack_year = cc_year,
        stack_label = cc_label,
        S_enter = as.numeric(W_pre %*% zz_enter),
        S_prior = as.numeric(W_pre %*% zz_prior)
      )
      pre_src <- pre_src[station_id %in% predetermined_receiver_ids]
      pre_src[, mapping := paste0("station_pre2013_sources_", r, "km")]
      pre_source_list[[kk]] <- pre_src[]
    }
    station_sources_pre <- data.table::rbindlist(
      pre_source_list, use.names = TRUE, fill = TRUE
    )
  }

  stack_station_pre_receivers_all_sources <- build_stack(
    dyn_panel_pre,
    station_sources,
    paste0("station_pre2013_receivers_all_sources_", r, "km")
  )
  stack_boundary_pre_receivers <- build_stack(
    dyn_panel_pre,
    boundary_sources,
    paste0("boundary_pre2013_receivers_", r, "km")
  )
  stack_station_pre_receivers_pre_sources <- if (nrow(station_sources_pre) > 0L) {
    build_stack(
      dyn_panel_pre,
      station_sources_pre,
      paste0("station_pre2013_receivers_pre2013_sources_", r, "km")
    )
  } else {
    data.table::data.table()
  }

  pred_stack_list <- list(
    station_full_sample = stack_station,
    boundary_full_sample = stack_boundary,
    station_pre2013_receivers_all_sources = stack_station_pre_receivers_all_sources,
    station_pre2013_receivers_pre2013_sources = stack_station_pre_receivers_pre_sources,
    boundary_pre2013_receivers = stack_boundary_pre_receivers
  )

  dynamic_pred_support <- data.table::rbindlist(lapply(names(pred_stack_list), function(nm) {
    dd <- pred_stack_list[[nm]]
    if (is.null(dd) || nrow(dd) == 0L) return(NULL)
    x <- support_one_stack(dd)
    x[, `:=`(
      variant = nm,
      policy_formation_cutoff = POLICY_FORMATION_CUTOFF
    )]
    x[]
  }), use.names = TRUE, fill = TRUE)

  pred_fit_list <- lapply(names(pred_stack_list), function(nm) {
    dd <- pred_stack_list[[nm]]
    if (is.null(dd) || nrow(dd) == 0L || data.table::uniqueN(dd$station_id) < 2L) {
      return(NULL)
    }
    tryCatch(
      fit_dynamic_mapping(dd, nm),
      error = function(e) {
        message("Dynamic predetermined-network fit failed for ", nm, ": ", e$message)
        NULL
      }
    )
  })
  names(pred_fit_list) <- names(pred_stack_list)

  dynamic_pred_coefs <- data.table::rbindlist(lapply(names(pred_fit_list), function(nm) {
    ff <- pred_fit_list[[nm]]
    if (is.null(ff)) return(NULL)
    x <- data.table::copy(ff$coefficients)
    x[, `:=`(
      variant = nm,
      policy_formation_cutoff = POLICY_FORMATION_CUTOFF
    )]
    x[]
  }), use.names = TRUE, fill = TRUE)

  dynamic_pred_implied <- data.table::rbindlist(lapply(names(pred_fit_list), function(nm) {
    ff <- pred_fit_list[[nm]]
    if (is.null(ff)) return(NULL)
    x <- data.table::copy(ff$implied)
    x[, `:=`(
      variant = nm,
      policy_formation_cutoff = POLICY_FORMATION_CUTOFF
    )]
    x[]
  }), use.names = TRUE, fill = TRUE)

  dynamic_pred_eff <- data.table::rbindlist(lapply(names(pred_fit_list), function(nm) {
    ff <- pred_fit_list[[nm]]
    dd <- pred_stack_list[[nm]]
    if (is.null(ff) || is.null(dd) || nrow(dd) == 0L) return(NULL)
    ee <- tryCatch(
      pooled_effective_support_13(dd, ff$model, nm, r),
      error = function(e) {
        message("Predetermined-network effective-support audit failed for ", nm, ": ", e$message)
        NULL
      }
    )
    if (is.null(ee)) return(NULL)
    ee[, policy_formation_cutoff := POLICY_FORMATION_CUTOFF]
    ee[]
  }), use.names = TRUE, fill = TRUE)

  data.table::fwrite(
    dynamic_pred_support,
    file.path(TAB_DIR, "table_13_20_dynamic_predetermined_network_support.csv")
  )
  data.table::fwrite(
    dynamic_pred_coefs,
    file.path(TAB_DIR, "table_13_21_dynamic_predetermined_network_coefficients.csv")
  )
  data.table::fwrite(
    dynamic_pred_implied,
    file.path(TAB_DIR, "table_13_22_dynamic_predetermined_network_implied_effects.csv")
  )
  data.table::fwrite(
    dynamic_pred_eff,
    file.path(TAB_DIR, "table_13_23_dynamic_predetermined_network_effective_support.csv")
  )

  saveRDS(stack_station, file.path(RDS_DIR, paste0("dynamic_stack_station_", r, "km.rds")))
  saveRDS(stack_boundary, file.path(RDS_DIR, paste0("dynamic_stack_boundary_", r, "km.rds")))
  saveRDS(fit_station$model, file.path(RDS_DIR, paste0("dynamic_model_station_", r, "km.rds")))
  saveRDS(fit_boundary$model, file.path(RDS_DIR, paste0("dynamic_model_boundary_", r, "km.rds")))

  # Dynamic mapping scatter by entering cohort.
  p_dynamic <- ggplot2::ggplot(
    corr_station_level,
    ggplot2::aes(x = S_enter_station, y = S_enter_boundary)
  ) +
    ggplot2::geom_point(alpha = 0.65) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~ stack_year) +
    ggplot2::labs(
      x = paste0("Station-based entering-cohort exposure (", r, " km)"),
      y = paste0("Boundary entering-area share within ", r, " km"),
      title = "Dynamic source-specific exposure: station versus policy-boundary representation"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  ggplot2::ggsave(
    file.path(FIG_DIR, "fig_13_02_dynamic_entering_station_vs_boundary_exposure.png"),
    p_dynamic, width = 9, height = 4.8, dpi = 300
  )
  ggplot2::ggsave(
    file.path(FIG_DIR, "fig_13_02_dynamic_entering_station_vs_boundary_exposure.pdf"),
    p_dynamic, width = 9, height = 4.8
  )
}

# ------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------

message("\n============================================================")
message("Pipeline 13 completed.")
message("Boundary-based ULEZ exposure diagnostic saved to:")
message("  ", OUT_DIR)
message("\nInterpretation rule:")
message("  Do not treat similarity as proof that either mapping is true.")
message("  Do not treat differences as automatic invalidation of the benchmark.")
message("  Compare mechanism, support, mapping correspondence, and estimates jointly.")

if (exists("static_total_treated")) {
  message("\nStatic total treated effect at observed mean treated exposure:")
  print(static_total_treated)
}
if (exists("dyn_combined_entrants")) {
  message("\nDynamic total entrant effects at observed mean entering-cohort exposure:")
  print(dyn_combined_entrants)
}
if (exists("pooled_eff_13")) {
  message("\nPooled effective identifying support: currently untreated exposure")
  print(
    pooled_eff_13[
      term == "post_untreated_S_enter",
      .(
        mapping, radius_km, raw_support_n, raw_nonzero_regressor_stations,
        n_stations_with_nonzero_residualized_target, n_eff,
        max_share, top2_share, max_share_station_id,
        beta_fwl, pooled_estimate, fwl_minus_pooled,
        residualization_status
      )
    ]
  )
}
if (exists("cohort_eff_13")) {
  message("\nCohort-specific effective identifying support: currently untreated exposure")
  print(
    cohort_eff_13[
      term == "post_untreated_S_enter",
      .(
        mapping, radius_km, stack_year, raw_support_n,
        raw_nonzero_regressor_stations,
        n_stations_with_nonzero_residualized_target,
        n_eff, max_share, top2_share, max_share_station_id,
        residualization_status
      )
    ][order(mapping, stack_year)]
  )
}
if (exists("network_composition")) {
  message("\nPredetermined monitoring-network composition:")
  print(network_composition)
}
if (exists("static_pred_total")) {
  message("\nStatic predetermined-network total treated effects:")
  print(static_pred_total)
}
if (exists("dynamic_pred_implied")) {
  message("\nDynamic predetermined-network: currently untreated exposure component:")
  print(dynamic_pred_implied[
    term == "post_untreated_S_enter",
    .(variant, mapping, radius_km, exposure_scale, estimate, std_error,
      implied_percent, ci_low_percent, ci_high_percent)
  ])
}
if (exists("dynamic_pred_eff")) {
  message("\nDynamic predetermined-network effective support: currently untreated exposure:")
  print(dynamic_pred_eff[
    term == "post_untreated_S_enter",
    .(mapping, radius_km, raw_support_n, raw_nonzero_regressor_stations,
      n_stations_with_nonzero_residualized_target, n_eff,
      max_share, top2_share, residualization_status)
  ])
}
message("============================================================")

# ============================================================
# Exportações canônicas do paper
# ============================================================

# 1. Correspondência entre station-based e boundary-based exposure.
map_paper <- data.table::copy(static_corr)
map_paper[, phase_order := match(phase, c("2019", "2021", "2023"))]
data.table::setorder(map_paper, phase_order)
map_paper[, phase_order := NULL]
agreement_txt <- paste(
  formatC(map_paper$binary_exposure_agreement, 3, format = "f"),
  collapse = ", "
)
writeLines(c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Correspondence between station-based and boundary-based 3 km exposure mappings}",
  "\\label{tab:ulez_boundary_mapping_correspondence}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabular}{lrrrrrr}", "\\toprule",
  "Policy phase & $N$ & Corr. & Mean $S^M$ & Mean $S^A$ & $S^M>0$ & $S^A>0$ \\\\",
  "\\midrule",
  paste0(
    map_paper$phase, " & ", map_paper$n_station_phase, " & ",
    formatC(map_paper$correlation, 3, format = "f"), " & ",
    formatC(map_paper$mean_station_exposure, 3, format = "f"), " & ",
    formatC(map_paper$mean_boundary_exposure, 3, format = "f"), " & ",
    map_paper$n_positive_station, " & ",
    map_paper$n_positive_boundary, " \\\\"
  ),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
  paste0(
    "\\item \\emph{Notes:} The table compares the station-based and policy-boundary exposure measures at each active ULEZ policy phase. ",
    "Binary exposed/unexposed agreement rates are ", agreement_txt,
    " in 2019, 2021, and 2023, respectively. ",
    "The comparison is a mapping diagnostic, not a causal-effect comparison."
  ),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_boundary_mapping_correspondence.tex"))

# 2. Support dinâmico da currently untreated exposure component.
cohort_support_paper <- cohort_eff_13[
  term == "post_untreated_S_enter",
  .(
    mapping,
    stack_year,
    raw_support_n,
    raw_nonzero_regressor_stations,
    n_eff
  )
]
st <- cohort_support_paper[grepl("^station_", mapping)]
bd <- cohort_support_paper[grepl("^boundary_", mapping)]
support_bnd_paper <- merge(st, bd, by = "stack_year", suffixes = c("_station", "_boundary"))
data.table::setorder(support_bnd_paper, stack_year)

pooled_support_paper <- pooled_eff_13[
  term == "post_untreated_S_enter",
  .(mapping, raw_support_n, raw_nonzero_regressor_stations, n_eff)
]
pooled_station <- pooled_support_paper[grepl("^station_", mapping)]
pooled_boundary <- pooled_support_paper[grepl("^boundary_", mapping)]
if (nrow(pooled_station) != 1L || nrow(pooled_boundary) != 1L) {
  stop("Expected one pooled effective-support row for each 3 km mapping.")
}

writeLines(c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Dynamic identifying support for the currently untreated exposure component}",
  "\\label{tab:ulez_boundary_dynamic_support}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabular}{lrrrrrr}", "\\toprule",
  "& \\multicolumn{2}{c}{Raw exposed receivers} & \\multicolumn{2}{c}{Nonzero raw target} & \\multicolumn{2}{c}{$N_{\\mathrm{eff}}$} \\\\",
  "Cohort & Station & Boundary & Station & Boundary & Station & Boundary \\\\",
  "\\midrule",
  paste0(
    support_bnd_paper$stack_year, " & ",
    support_bnd_paper$raw_support_n_station, " & ", support_bnd_paper$raw_support_n_boundary, " & ",
    support_bnd_paper$raw_nonzero_regressor_stations_station, " & ",
    support_bnd_paper$raw_nonzero_regressor_stations_boundary, " & ",
    formatC(support_bnd_paper$n_eff_station, 2, format = "f"), " & ",
    formatC(support_bnd_paper$n_eff_boundary, 2, format = "f"), " \\\\"
  ),
  "\\midrule",
  paste0(
    "Pooled & ",
    pooled_station$raw_support_n, " & ", pooled_boundary$raw_support_n, " & ",
    pooled_station$raw_nonzero_regressor_stations, " & ",
    pooled_boundary$raw_nonzero_regressor_stations, " & ",
    formatC(pooled_station$n_eff, 2, format = "f"), " & ",
    formatC(pooled_boundary$n_eff, 2, format = "f"), " \\\\"
  ),
  "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
  "\\item \\emph{Notes:} ``Raw exposed receivers'' counts currently untreated stations with positive entering-cohort exposure; ``Nonzero raw target'' counts stations for which the corresponding post-treatment target regressor is nonzero. $N_{\\mathrm{eff}}$ is computed from station shares of residualized identifying variation after partialling out fixed effects and the remaining decomposition terms. The boundary mapping expands pooled support but does not resolve the very weak 2023 cohort-specific support.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_boundary_dynamic_support.tex"))

# 3. Efeitos avaliados no suporte observado.
untreated_effect <- dyn_implied[term == "post_untreated_S_enter", .(
  mapping,
  untreated_pct = implied_percent,
  untreated_lo = ci_low_percent,
  untreated_hi = ci_high_percent
)]
entrant0 <- dyn_combined_entrants[entrant_group == "Entrants not previously exposed", .(
  mapping,
  entrant0_pct = implied_percent,
  entrant0_lo = ci_low_percent,
  entrant0_hi = ci_high_percent
)]
entrant1 <- dyn_combined_entrants[entrant_group == "Entrants previously exposed", .(
  mapping,
  entrant1_pct = implied_percent,
  entrant1_lo = ci_low_percent,
  entrant1_hi = ci_high_percent
)]
static_eff <- static_total_treated[, .(
  mapping,
  static_pct = implied_percent,
  static_lo = ci_low_percent,
  static_hi = ci_high_percent
)]
effects_paper <- Reduce(function(x, y) merge(x, y, by = "mapping"),
                        list(static_eff, entrant0, entrant1, untreated_effect))
effects_paper[, mapping_label := data.table::fifelse(
  grepl("^station", mapping), "Station-based", "Boundary-based"
)]
fmt_effect_ci <- function(est, lo, hi) {
  paste0(
    "\\shortstack{\\(", formatC(est, 2, format = "f"), "\\%\\)\\\\{\\scriptsize \\([",
    formatC(lo, 2, format = "f"), ",", formatC(hi, 2, format = "f"), "]\\)}}"
  )
}
writeLines(c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Effects evaluated at observed exposure support under alternative 3 km mappings}",
  "\\label{tab:ulez_boundary_effects}", "\\small", "\\begin{threeparttable}",
  "\\begin{tabularx}{0.96\\textwidth}{>{\\raggedright\\arraybackslash}X >{\\centering\\arraybackslash}m{3.4cm} >{\\centering\\arraybackslash}m{3.4cm}}",
  "\\toprule", "Estimand & Station-based & Boundary-based \\\\", "\\midrule",
  paste0("Static total treated effect & ",
         fmt_effect_ci(effects_paper[mapping_label == "Station-based", static_pct], effects_paper[mapping_label == "Station-based", static_lo], effects_paper[mapping_label == "Station-based", static_hi]), " & ",
         fmt_effect_ci(effects_paper[mapping_label == "Boundary-based", static_pct], effects_paper[mapping_label == "Boundary-based", static_lo], effects_paper[mapping_label == "Boundary-based", static_hi]), " \\\\"),
  paste0("Dynamic entrants, not previously exposed & ",
         fmt_effect_ci(effects_paper[mapping_label == "Station-based", entrant0_pct], effects_paper[mapping_label == "Station-based", entrant0_lo], effects_paper[mapping_label == "Station-based", entrant0_hi]), " & ",
         fmt_effect_ci(effects_paper[mapping_label == "Boundary-based", entrant0_pct], effects_paper[mapping_label == "Boundary-based", entrant0_lo], effects_paper[mapping_label == "Boundary-based", entrant0_hi]), " \\\\"),
  paste0("Dynamic entrants, previously exposed & ",
         fmt_effect_ci(effects_paper[mapping_label == "Station-based", entrant1_pct], effects_paper[mapping_label == "Station-based", entrant1_lo], effects_paper[mapping_label == "Station-based", entrant1_hi]), " & ",
         fmt_effect_ci(effects_paper[mapping_label == "Boundary-based", entrant1_pct], effects_paper[mapping_label == "Boundary-based", entrant1_lo], effects_paper[mapping_label == "Boundary-based", entrant1_hi]), " \\\\"),
  paste0("Dynamic currently untreated exposure & ",
         fmt_effect_ci(effects_paper[mapping_label == "Station-based", untreated_pct], effects_paper[mapping_label == "Station-based", untreated_lo], effects_paper[mapping_label == "Station-based", untreated_hi]), " & ",
         fmt_effect_ci(effects_paper[mapping_label == "Boundary-based", untreated_pct], effects_paper[mapping_label == "Boundary-based", untreated_lo], effects_paper[mapping_label == "Boundary-based", untreated_hi]), " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\begin{tablenotes}", "\\footnotesize",
  "\\item \\emph{Notes:} Entries are implied percentage effects with 95 percent confidence intervals in brackets. Effects are evaluated at the observed mean exposure relevant to each estimand. Because the two mappings generate different exposure distributions, the magnitudes are mapping-specific and should not be interpreted as estimates of a single invariant parameter.",
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_boundary_effects.tex"))

# 4. Predetermined-monitoring-network diagnostic.
# Combina total treated effect (static) com currently untreated exposure (dynamic).
pred_u <- dynamic_pred_implied[term == "post_untreated_S_enter", .(
  variant,
  untreated_pct = implied_percent,
  untreated_lo = ci_low_percent,
  untreated_hi = ci_high_percent
)]
pred_s <- static_pred_total[, .(
  variant = mapping,
  treated_pct = implied_percent,
  treated_lo = ci_low_percent,
  treated_hi = ci_high_percent
)]
# A denominação do variant estático usa "original_sources"; harmonizar antes do merge.
pred_s[variant == "station_pre2013_receivers_original_sources",
       variant := "station_pre2013_receivers_all_sources"]
pred_eff <- dynamic_pred_eff[term == "post_untreated_S_enter", .(
  variant = mapping,
  raw_support = raw_support_n,
  n_eff
)]
pred <- merge(pred_s, pred_u, by = "variant", all = TRUE)
pred <- merge(pred, pred_eff, by = "variant", all = TRUE)
receiver_n <- c(
  station_full_sample = data.table::uniqueN(raw_panel$station_id),
  boundary_full_sample = data.table::uniqueN(raw_panel$station_id),
  station_pre2013_receivers_all_sources = length(predetermined_receiver_ids),
  station_pre2013_receivers_pre2013_sources = length(predetermined_receiver_ids),
  boundary_pre2013_receivers = length(predetermined_receiver_ids)
)
variant_label <- c(
  station_full_sample = "Station, full",
  boundary_full_sample = "Boundary, full",
  station_pre2013_receivers_original_sources = "Station, pre-2013 receivers, all sources",
  station_pre2013_receivers_all_sources = "Station, pre-2013 receivers, all sources",
  station_pre2013_receivers_pre2013_sources = "Station, pre-2013 receivers, pre-2013 sources",
  boundary_pre2013_receivers = "Boundary, pre-2013 receivers"
)
pred[, specification := unname(variant_label[variant])]
pred[, receivers := unname(receiver_n[variant])]
order_variants <- c(
  "station_full_sample", "boundary_full_sample",
  "station_pre2013_receivers_all_sources",
  "station_pre2013_receivers_pre2013_sources",
  "boundary_pre2013_receivers"
)
pred[, ord := match(variant, order_variants)]
data.table::setorder(pred, ord)
writeLines(c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Predetermined-monitoring-network diagnostic}",
  "\\label{tab:ulez_predetermined_network}",
  "\\small",
  "\\resizebox{0.96\\textwidth}{!}{%",
  "\\begin{tabular}{>{\\raggedright\\arraybackslash}m{5.0cm} >{\\centering\\arraybackslash}m{1.8cm} >{\\centering\\arraybackslash}m{3.0cm} >{\\centering\\arraybackslash}m{3.0cm} >{\\centering\\arraybackslash}m{1.7cm} >{\\centering\\arraybackslash}m{1.7cm}}",
  "\\toprule",
  "Specification & Receivers & Total treated & Currently untreated exposure & Raw support & $N_{\\mathrm{eff}}$ \\\\",
  "\\midrule",
  paste0(
    pred$specification, " & ", pred$receivers, " & ",
    "\\shortstack{\\(", formatC(pred$treated_pct, 2, format = "f"), "\\%\\)\\\\{\\scriptsize \\([",
    formatC(pred$treated_lo, 2, format = "f"), ",", formatC(pred$treated_hi, 2, format = "f"), "]\\)}} & ",
    "\\shortstack{\\(", formatC(pred$untreated_pct, 2, format = "f"), "\\%\\)\\\\{\\scriptsize \\([",
    formatC(pred$untreated_lo, 2, format = "f"), ",", formatC(pred$untreated_hi, 2, format = "f"), "]\\)}} & ",
    pred$raw_support, " & ", formatC(pred$n_eff, 2, format = "f"), " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "",
  "\\vspace{0.4em}",
  "\\parbox{0.96\\textwidth}{\\footnotesize \\emph{Notes:} Entries in the two effect columns report implied percentage effects, with 95 percent confidence intervals in brackets. Raw support and $N_{\\mathrm{eff}}$ refer to the currently untreated exposure component. Unknown opening dates are excluded from the pre-policy network. The pre-policy restriction is a diagnostic for policy-responsive monitor placement; it does not assume that historical station placement is spatially random or exogenous to pollution and traffic conditions.}",
  "\\end{table}"
), file.path(TABLE_DIR, "table_ulez_predetermined_network.tex"))

saveRDS(
  list(
    static_mapping_correspondence = map_paper,
    dynamic_support = support_bnd_paper,
    effects = effects_paper,
    predetermined_network = pred,
    pooled_effective_support = pooled_eff_13,
    cohort_effective_support = cohort_eff_13,
    network_composition = network_composition
  ),
  BOUNDARY_RESULTS_RDS
)
message("[04] Diagnósticos de boundary/network concluídos.")
