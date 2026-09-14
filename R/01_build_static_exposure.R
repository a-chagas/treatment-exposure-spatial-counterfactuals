# ============================================================
# 01_build_static_exposure.R
# Reconstrói o painel estático e todas as exposições usadas no paper
# ============================================================

message("[01] Construindo painel estático e matrizes de exposição...")

raw_panel <- data.table::fread(RAW_PANEL_FILE)
assign <- data.table::fread(ASSIGN_FILE)
data.table::setnames(raw_panel, clean_names_simple(names(raw_panel)))
data.table::setnames(assign, clean_names_simple(names(assign)))

raw_panel[, code := toupper(trimws(as.character(code)))]
raw_panel[, month := as.Date(month)]
assign[, code := toupper(trimws(as.character(code)))]
assign[, latitude := as.numeric(latitude)]
assign[, longitude := as.numeric(longitude)]
assign[, inside_ulez_2019 := to01(inside_ulez_2019)]
assign[, inside_ulez_2021 := to01(inside_ulez_2021)]
assign[, inside_ulez_2023 := to01(inside_ulez_2023)]
assign <- unique(assign, by = "code")

receiver_ids <- sort(unique(raw_panel$code))
receiver_meta <- assign[match(receiver_ids, code)]
if (anyNA(receiver_meta$code) || any(!is.finite(receiver_meta$latitude)) ||
    any(!is.finite(receiver_meta$longitude))) {
  stop("Não foi possível associar coordenadas válidas a todas as estações do painel.")
}

# Convenção estática do artigo: W é definido no universo das estações com outcome.
receiver_sf <- sf::st_as_sf(
  receiver_meta,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
receiver_sf <- sf::st_transform(receiver_sf, TARGET_CRS)
dist_km <- as.matrix(units::drop_units(sf::st_distance(receiver_sf))) / 1000
rownames(dist_km) <- receiver_ids
colnames(dist_km) <- receiver_ids

row_standardize <- function(W) {
  rs <- rowSums(W, na.rm = TRUE)
  ok <- rs > 0
  W[ok, ] <- W[ok, , drop = FALSE] / rs[ok]
  W[!ok, ] <- 0
  W
}

make_W_threshold <- function(upper_km) {
  W <- 1 * (dist_km > 0 & dist_km <= upper_km)
  dimnames(W) <- dimnames(dist_km)
  row_standardize(W)
}

make_W_ring <- function(lower_km, upper_km) {
  W <- 1 * (dist_km > lower_km & dist_km <= upper_km)
  dimnames(W) <- dimnames(dist_km)
  row_standardize(W)
}

W_list <- list(
  bench_0_3km = make_W_threshold(3),
  near_2_5km = make_W_ring(2, 5),
  far_10_20km = make_W_ring(10, 20),
  bench_0_2km = make_W_threshold(2),
  bench_0_5km = make_W_threshold(5),
  near_0_2km = make_W_ring(0, 2),
  medium_5_10km = make_W_ring(5, 10),
  far_20_40km = make_W_ring(20, 40)
)

months <- sort(unique(raw_panel$month))
D19 <- outer(receiver_meta$inside_ulez_2019, months >= STATIC_STAGE_START[1], `*`)
D21 <- outer(receiver_meta$inside_ulez_2021, months >= STATIC_STAGE_START[2], `*`)
D23 <- outer(receiver_meta$inside_ulez_2023, months >= STATIC_STAGE_START[3], `*`)
Dmat <- 1 * ((D19 + D21 + D23) > 0)
rownames(Dmat) <- receiver_ids
colnames(Dmat) <- as.character(months)

exposure_long <- lapply(names(W_list), function(nm) {
  S <- W_list[[nm]] %*% Dmat
  data.table::data.table(
    code = rep(receiver_ids, times = length(months)),
    month = rep(months, each = length(receiver_ids)),
    value = as.vector(S)
  )[, variable := paste0("S_", nm)]
})
exposure_long <- data.table::rbindlist(exposure_long)
exposure_wide <- data.table::dcast(exposure_long, code + month ~ variable, value.var = "value")

D_long <- data.table::data.table(
  code = rep(receiver_ids, times = length(months)),
  month = rep(months, each = length(receiver_ids)),
  D_ulez_2019 = as.integer(as.vector(D19)),
  D_ulez_2021 = as.integer(as.vector(D21)),
  D_ulez_2023 = as.integer(as.vector(D23)),
  D_active = as.integer(as.vector(Dmat))
)

panel <- merge(raw_panel, D_long, by = c("code", "month"), all.x = TRUE, sort = FALSE)
panel <- merge(panel, exposure_wide, by = c("code", "month"), all.x = TRUE, sort = FALSE)
panel[, log_y := log(as.numeric(no2))]
panel[, `:=`(
  S_bench = S_bench_0_3km,
  S_near = S_near_2_5km,
  S_far = S_far_10_20km
)]
panel[, `:=`(
  U_bench = (1 - D_active) * S_bench,
  U_near = (1 - D_active) * S_near,
  U_far = (1 - D_active) * S_far,
  T_bench = D_active * S_bench,
  T_near = D_active * S_near,
  T_far = D_active * S_far
)]

# Ordem preservada para facilitar a comparação com os arquivos usados na submissão.
col_order <- c(
  "code", "month", "no2", "n_hours", "D_active",
  "D_ulez_2019", "D_ulez_2021", "D_ulez_2023",
  "S_bench_0_3km", "S_near_2_5km", "S_far_10_20km",
  "S_bench_0_2km", "S_bench_0_5km", "S_near_0_2km",
  "S_medium_5_10km", "S_far_20_40km", "log_y",
  "S_bench", "S_near", "S_far", "U_bench", "U_near", "U_far",
  "T_bench", "T_near", "T_far"
)
data.table::setcolorder(panel, col_order)
data.table::setorder(panel, code, month)

# Duas verificações estruturais mínimas, usadas apenas para impedir uma reprodução silenciosamente divergente.
stopifnot(nrow(panel) == 10138L)
stopifnot(data.table::uniqueN(panel$code) == 129L)
stopifnot(data.table::uniqueN(panel$month) == 108L)
stopifnot(all(panel$D_active %in% c(0L, 1L)))

data.table::fwrite(panel, STATIC_PANEL_FILE)
saveRDS(W_list, file.path(DERIVED_DIR, "static_W_matrices.rds"))
message("[01] Painel estático salvo em: ", STATIC_PANEL_FILE)
