# Release validation

## Final end-to-end reproduction

The release candidate was executed from the repository root with:

```r
source("run_all.R")
```

The run completed successfully with:

```text
[06] 24 paper outputs verified.

Reprodução concluída. Ver outputs/tables e outputs/figures.
```

The canonical paper-output set contains 18 tables and 6 figures. The final export step checks for missing, unexpected, and stale paper artifacts before rewriting `outputs/manifest.csv`.

## Static reconstruction

The static exposure reconstruction was independently checked against the analysis panel used in the revision. Starting only from the monthly NO2 panel and station assignment file, the reconstructed treatment indicator and all eight reported station-based exposure variables match the analysis panel up to floating-point precision:

- maximum absolute difference in `D_active`: `0`;
- maximum absolute difference across exposure variables in an independent Python/PROJ reconstruction: below `2.0e-15`.

The reconstructed sample matches the manuscript panel dimensions: 10,138 station-month observations, 129 monitoring stations, and 108 months.

The GeoPackage inputs were verified as valid OGC GeoPackage/SQLite files. The active R scripts contain no machine-specific `G:/...` or `C:/...` paths, no `setwd()`, no network downloads, and no automatic package installation.

## Spatial inference

The historical spatial-inference implementation was recovered and reconciled with the consolidated pipeline. The manuscript Conley standard errors are reproduced by `fixest::vcov_conley(..., distance = "triangular")`; the spherical-distance option does not reproduce them.

The public spatial block bootstrap uses a prospectively deterministic implementation:

- eight station blocks from k-means on scaled raw longitude/latitude;
- `nstart = 50`, `algorithm = "Hartigan-Wong"`;
- separate fixed seeds for k-means and bootstrap resampling;
- k-means labels canonicalized by centroid longitude and then latitude;
- target of 499 valid bootstrap replications;
- non-estimable resamples discarded and deterministically replaced by the next draw from the same RNG stream.

The validated canonical row is:

- direct-treatment component: `-6.60 (4.67)`;
- treated-exposure component: `-3.31 (4.84)`;
- currently untreated exposure component: `-9.85*** (3.32)`.

These replace the previous finite Monte Carlo realization. Point estimates and substantive inference are unchanged.

## Exposure-bin diagnostic

The recovered implementation first attempts quartile cutpoints among observations with positive 0–3 km exposure. In the benchmark panel, tied quantiles make the five candidate cutpoints non-unique, activating the prespecified four-bin equal-width fallback.

The realized bins are labeled `B1`–`B4`. Zero exposure remains the explicit reference condition, and fail-fast checks guard against missing bin indicators or unintended observation deletion.

Validated support and coefficients are:

- B1: treated support 0; currently untreated support 108; currently untreated exposure `-3.70 (3.99)`;
- B2: treated support 22, `0.74 (3.04)`; currently untreated support 72, `-7.18** (2.74)`;
- B3: treated support 296, `-6.97 (4.16)`; currently untreated support 55, `-13.07* (7.24)`;
- B4: treated support 2,264, `-4.55 (3.07)`; currently untreated support 151, `-6.95*** (2.30)`.

## Dynamic and boundary diagnostics

The stacked dynamic implementation uses calendar-time risk sets: future-treated receivers remain currently untreated only while `t < C_i`. The final outputs include pooled and cohort-specific raw and residualized support diagnostics.

The alternative policy-boundary representation and the predetermined-monitoring-network diagnostic are generated from the public processed inputs. The final table generators include the pooled support rows and explanatory notes used in the revised manuscript.

## Manuscript-output reconciliation

The 18 generated empirical tables and 6 generated figures were reconciled with the final clean manuscript assets. The manuscript uses renamed copies whose filenames begin with their final paper numbers; the reproduction repository retains descriptive canonical filenames so that the computational pipeline remains stable.

The remaining manuscript tables are conceptual/methodological tables and are not outputs of the empirical reproduction pipeline.

## Release items

The computational and empirical audit is complete. Before making the repository public, finalize the software licence and commit an `renv.lock` created from the verified release R environment. Source-data redistribution terms are documented in `data/README_data.md`.
