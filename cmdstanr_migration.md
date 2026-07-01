# Migrating multip2 from rstan to cmdstanr

*Assessment written 2026-07-01. Based on the package state on branch `features/random_slopes` (post-merge with master).*

## Verdict

Very feasible — the rstan coupling is much narrower than typical. There are only 5 direct rstan call sites in `R/`, and the package's own wrapper functions (`extract_draws`, `create_summary`) already funnel most access through choke points. Realistic effort: **roughly a week of focused work**, most of it validation rather than rewriting. Two structural snags need to be resolved (see below).

## Why cmdstanr is genuinely better for this package

- **Current Stan.** CRAN rstan is stuck around Stan 2.32; CmdStan is at 2.36+. Years of speedups, better diagnostics, and `reduce_sum` threading (relevant for the dyad loop, which is embarrassingly parallel).
- **Stability.** No StanHeaders/R-version compatibility breakage; sampling runs in a separate process, so a model crash cannot kill the R session.
- **Data provenance.** cmdstanr persists the data JSON and output CSVs next to every run. The "can I recover stan_data from a saved fit?" problem (hit during the PSIS-LOO work, 2026-06) would not exist.
- **Standalone generated quantities** are a first-class citizen: `mod$generate_quantities(fitted_params = fit, data = ...)` — cleaner than `rstan::gqs()` as used in `inst/replication_code/applied_multilevel_bully/PSIS-LOO.R`.
- **Lower memory.** Draws stream to CSV instead of accumulating in RAM (current fits are 265MB+ in-memory objects).
- **Repo hygiene bonus:** the five 20–40MB `.rds` files in `inst/stan/` (~180MB total) are rstan `auto_write` caches. They are rstan-specific and become deletable, shrinking the repo dramatically.

## The two real snags

### 1. `include = FALSE` has no cmdstanr equivalent

`fit.Mp2Model()` (R/classes.R) and the replication scripts rely on `par = "x_beta", include = FALSE` to avoid saving x_beta. CmdStan writes **all** transformed parameters to its output CSVs — no exclusion option. x_beta is `sum(N) × K` (thousands of dyads × 16 outcomes), so the CSVs would balloon to gigabytes and I/O would dominate.

**Fix:** refactor the model — move the x_beta computation into a user-defined function, make x_beta a *local* variable in the `model` block, and call the same function in generated quantities. Mechanical but touches the heart of every model file; needs numerical validation afterward. **Estimate: ~1 day including checks.**

### 2. Old array syntax

All active models use pre-2.26 syntax (`int y[N]`, `real x[2]` — 90 such declarations in `multilevel_multiplex_p2_uncentered.stan` alone). Stan 2.33+ **removed** this syntax, so current CmdStan will not compile the files at all.

**Fix:** automated — `stanc --canonicalize deprecations` rewrites files in place (cmdstanr can invoke it). Diff and recompile-check each file. Only 2–3 of the 12 Stan files are actively used; migrate those, leave the rest. **Estimate: half a day.**

## The mechanical rest (small)

| Site | Change |
|---|---|
| `R/classes.R:268` `rstan::stan(...)` | `cmdstan_model(fpath)$sample(...)`; remap `iter`/`warmup` → `iter_sampling`/`iter_warmup`, `mc.cores` → `parallel_chains` |
| `R/classes.R:392` `rstan::summary(fit)$summary` | `fit$summary()` via posterior — column names/order differ (`n_eff`/`Rhat` → `ess_bulk`/`rhat`, different quantiles), and `extract_model_info` slices columns `1:10` by position, so the summary-building functions need a careful pass |
| `R/classes.R:554` `rstan::extract(fit)` | `posterior::extract_variable_array(fit$draws(par))` — one wrapper, everything else goes through it |
| `R/simulations.R:41` `expose_stan_functions` | `cmdstan_model(..., compile_standalone = TRUE)$functions$...` — direct equivalent |
| `R/convergence_check.R` `ggmcmc::ggs(stanfit)` | ggs has no CmdStanMCMC method; convert with `posterior::as_draws_df` → coda `mcmc.list` first (small shim) |

Also update `DESCRIPTION` (drop `rstan` from Imports; add `cmdstanr` — see CRAN note below — plus `posterior`) and the `@import rstan` roxygen tags.

## Cons to weigh

- **cmdstanr is not on CRAN** (R-universe only), and each user must run `install_cmdstan()` once. If multip2 should ever go to CRAN, cmdstanr can only be an optional backend (the brms pattern), which roughly doubles the migration work. For a GitHub research package, a non-issue.
- **Fits become file-backed.** A `CmdStanMCMC` object references CSV files; naive `saveRDS()` breaks unless `fit$save_object()` is called first. The current workflow saves fits to `.Rds` constantly, so this is a habit change (arguably an improvement — `save_object` is explicit).
- **Old saved fits don't interoperate.** Existing `stanfit` objects (all the bully/gossip runs) can still be read by old code and by `posterior::as_draws`, but package functions rewritten for cmdstanr won't accept them without a small compatibility branch. Keep `extract_draws`/`create_summary` tolerant of both classes during the transition.

## Recommended sequence

Clean switch, not the dual-backend pattern (overkill for a research package):

1. Canonicalize the 2–3 active Stan files (`stanc --canonicalize deprecations`), diff, recompile-check.
2. x_beta refactor (local variable + shared function), with a numerical equivalence check against an existing rstan fit.
3. Swap the 5 rstan call sites in `R/`.
4. Fix summary column handling (`extract_model_info`, `make_fixed_summary`, `make_random_summary`).
5. ggmcmc shim in `convergence_check.R`.
6. Delete the `inst/stan/*.rds` auto_write caches; update `DESCRIPTION`/roxygen.
7. Run `tests/testthat/` (esp. `test-fit_stan.R`) as the safety net.

**Timing:** hold off until the current LOO analysis is done — mid-analysis is the wrong time to swap inference backends, and `PSIS-LOO.R` as written depends on `rstan::gqs()` against rstan-saved fits.
