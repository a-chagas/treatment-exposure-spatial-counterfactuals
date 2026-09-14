# ============================================================
# 00_config.R
# Configuração comum do pacote de reprodução
# ============================================================

# O repositório deve ser executado a partir de sua raiz.
REPRO_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(REPRO_ROOT, "run_all.R"))) {
  stop("Execute run_all.R a partir da raiz do repositório de reprodução.")
}

DATA_DIR    <- file.path(REPRO_ROOT, "data", "processed", "london_ulez_air")
DERIVED_DIR <- file.path(REPRO_ROOT, "data", "derived")
TABLE_DIR   <- file.path(REPRO_ROOT, "outputs", "tables")
FIGURE_DIR  <- file.path(REPRO_ROOT, "outputs", "figures")
LOG_DIR     <- file.path(REPRO_ROOT, "outputs", "logs")

dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

REQUIRED_PACKAGES <- c(
  "data.table", "dplyr", "tidyr", "tibble", "readr", "stringr",
  "ggplot2", "fixest", "sf", "lubridate", "units", "glue", "janitor"
)
missing_packages <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Pacotes R ausentes: ", paste(missing_packages, collapse = ", "),
    ". Instale-os antes de executar a reprodução."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(fixest)
  library(sf)
  library(lubridate)
  library(units)
  library(glue)
  library(janitor)
})

fixest::setFixest_notes(FALSE)
sf::sf_use_s2(FALSE)

EXPOSURE_EPS <- 1e-10
TARGET_CRS <- 27700L
BENCHMARK_RADIUS_KM <- 3
DYNAMIC_RADII_KM <- c(3)
POLICY_FORMATION_CUTOFF <- as.Date("2013-02-13")

# Datas efetivas da política e primeiro mês-calendário completo usado no painel estático.
EXPANSION_DATES <- as.Date(c("2019-04-01", "2021-10-01", "2023-08-01"))
STATIC_STAGE_START <- as.Date(c("2019-05-01", "2021-11-01", "2023-09-01"))
EXPANSION_LABELS <- c(
  "Central ULEZ, Apr. 2019",
  "Inner London expansion, Oct. 2021",
  "London-wide expansion, Aug. 2023"
)

RAW_PANEL_FILE <- file.path(DATA_DIR, "laqn_full_monthly_no2_2016_2024.csv")
ASSIGN_FILE <- file.path(DATA_DIR, "laqn_station_ulez_assignment.csv")
BOUNDARY_FILES <- c(
  `2019` = file.path(DATA_DIR, "ulez_2019.gpkg"),
  `2021` = file.path(DATA_DIR, "ulez_2021.gpkg"),
  `2023` = file.path(DATA_DIR, "ulez_2023.gpkg")
)
STATIC_PANEL_FILE <- file.path(DERIVED_DIR, "analysis_panel_full_with_exposures.csv")
STATIC_RESULTS_RDS <- file.path(DERIVED_DIR, "static_results.rds")
DYNAMIC_RESULTS_RDS <- file.path(DERIVED_DIR, "dynamic_results.rds")
BOUNDARY_RESULTS_RDS <- file.path(DERIVED_DIR, "boundary_network_results.rds")

required_files <- c(RAW_PANEL_FILE, ASSIGN_FILE, BOUNDARY_FILES)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Arquivos de dados ausentes:\n", paste(missing_files, collapse = "\n"))
}

clean_names_simple <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "true", "t", "yes", "y", "sim")] <- 1L
  out[z %in% c("0", "false", "f", "no", "n", "nao", "não")] <- 0L
  out
}

month_id <- function(date) {
  lubridate::year(date) * 12L + lubridate::month(date)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

pct_effect <- function(beta) 100 * (exp(beta) - 1)

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}
fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "$<0.001$", formatC(p, digits = 3, format = "f")))
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
