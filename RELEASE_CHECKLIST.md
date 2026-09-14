# GitHub release checklist

- [x] One-command entry point: `source("run_all.R")`.
- [x] No machine-specific paths, `setwd()`, automatic installs, or network downloads in active scripts.
- [x] Public processed inputs included and documented.
- [x] Source-data redistribution terms documented.
- [x] Final end-to-end run completed successfully.
- [x] Canonical validation: 18 tables + 6 figures = 24 outputs.
- [x] Missing, unexpected, and stale output checks active.
- [x] Deterministic spatial bootstrap validated.
- [x] Triangular Conley convention validated.
- [x] Equal-width exposure-bin fallback validated and documented.
- [x] Final empirical manuscript assets reconciled with repository outputs.
- [ ] Create `renv.lock` from the verified release R environment and commit it.
- [ ] Add the selected software `LICENSE`.
- [ ] Create GitHub repository and push the release candidate.
- [ ] Run `source("run_all.R")` once from a fresh clone and confirm `[06] 24 paper outputs verified.`
- [ ] Tag the validated release (suggested: `v1.0.0`).
