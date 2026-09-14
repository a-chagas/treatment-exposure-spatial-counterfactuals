# Reproduction package — Treatment, Exposure, and Spatial Counterfactuals

This repository reproduces the empirical tables and figures reported in the London ULEZ application of **“Treatment, Exposure, and Spatial Counterfactuals in Regional Policy Evaluation: A Framework for Impact Analysis with Spatial Spillovers.”**

## One-command reproduction

From the repository root, run:

```r
source("run_all.R")
```

The pipeline reconstructs treatment and exposure variables from the processed monitoring data and ULEZ boundary files, estimates the paper specifications, and writes the paper-ready empirical artifacts to:

- `outputs/tables/`
- `outputs/figures/`

`outputs/manifest.csv` lists the canonical empirical paper artifacts. The final export step fails if an expected artifact is missing, if an unexpected artifact is present in either paper-output directory, or if an expected artifact is stale relative to the current run.

A complete release-environment run has been verified successfully. The terminal validation message is:

```text
[06] 24 paper outputs verified.
```

The canonical set contains 18 tables and 6 figures.

## Repository structure

- `run_all.R` — single entry point for the complete reproduction.
- `R/00_config.R` — paths, package checks, common constants, and helpers.
- `R/01_build_static_exposure.R` — reconstructs treatment and the eight station-based static exposure mappings from the processed data.
- `R/02_static_results.R` — benchmark/static specifications, support and design diagnostics, pre-trend figures, restricted-vs-unrestricted comparison, and across-matrix results.
- `R/03_dynamic_results.R` — stacked dynamic decomposition, calendar-time risk set, effective identifying support, event-study figures, and announcement/timing diagnostics.
- `R/04_boundary_network_results.R` — boundary-area exposure mapping and the pre-2013 predetermined-monitoring-network robustness exercise.
- `R/05_appendix_robustness.R` — descriptive statistics, implied-effect tables, conventional robustness, spatial-inference diagnostics, and the exposure-bin diagnostic.
- `R/06_export_paper_outputs.R` — final export and canonical-output validation.
- `data/processed/` — public processed inputs required by the reproduction pipeline.
- `data/derived/` — regenerated intermediate files; excluded from version control except for `.gitkeep`.
- `outputs/` — canonical paper tables, figures, and manifest.
- `VALIDATION.md` — detailed audit and validation record.
- `data/README_data.md` — data provenance, licensing, and attribution.

## Important implementation conventions

The static benchmark uses a row-standardized binary 0–3 km matrix defined on the monitoring stations observed in the outcome panel. Alternative static mappings are 0–2 and 0–5 km thresholds and 0–2, 2–5, 5–10, 10–20, and 20–40 km rings.

The dynamic implementation defines exposure sources on the complete assignment-station universe. For an expansion cohort `c`, future-treated receivers remain in the currently untreated risk set only while `t < C_i`; observations are censored from an earlier stack once the receiver's own treatment begins.

The boundary diagnostic defines local exposure as the fraction of a 3 km buffer around an outcome station that overlaps the active ULEZ area. The predetermined-network diagnostic classifies a monitor as predetermined when its opening date is strictly before 13 February 2013; monitors with unknown opening dates are excluded from that restricted set.

For the exposure-response functional-form diagnostic, the implementation first attempts quartile cutpoints among positive benchmark exposures. Because tied positive-exposure quantiles make those cutpoints non-unique in the ULEZ sample, the prespecified fallback of four equal-width bins is used. The realized bins are reported as `B1`–`B4`.

For spatial inference, Conley covariance estimates use the triangular distance convention. The eight-block spatial bootstrap uses separate deterministic seeds for geographic blocking and resampling, canonicalizes the arbitrary k-means cluster labels, and retains 499 valid bootstrap replications. Non-estimable resamples are discarded and replaced by the next draw from the same deterministic RNG stream.

## Dependencies

The scripts require R and the following packages:

`data.table`, `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`, `ggplot2`, `fixest`, `sf`, `lubridate`, `units`, `glue`, `janitor`, and `purrr`.

The scripts deliberately do **not** install packages automatically.

For an environment-pinned reproduction, restore the committed `renv.lock` when available:

```r
install.packages("renv")
renv::restore()
```

## Data and licensing

The processed monitoring and ULEZ boundary inputs redistributed under `data/` are derived from sources made available under the UK Open Government Licence v2.0. See `data/README_data.md` for source-specific details and attribution.

The data licence does not automatically license the R code. See `LICENSE_NOTE.md` and the repository `LICENSE` file once the software licence is finalized.

## Validation

The release pipeline has been executed end-to-end from the repository root and verified to regenerate the complete canonical set of 24 empirical paper artifacts. The static exposure reconstruction, dynamic risk-set implementation, policy-boundary diagnostic, deterministic spatial block bootstrap, triangular Conley covariance estimates, and equal-width exposure-bin fallback have each been audited separately.

See `VALIDATION.md` for the detailed record.
