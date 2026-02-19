
# MSSTNB Algorithm 2 sampler (R)

This folder contains an R implementation of **Algorithm 2** (one blocked Gibbs iteration) for the MSSTNB model with regression and ICAR.

It is designed to work directly with the `.rds` files written by your **01_generate_all_data.R** driver:
`data_sim/<scenario_id>/rep0001.rds`, etc.

## Folder layout

- `R/` core functions
  - `sampler_algo2.R` the blocked Gibbs driver (Algorithm 2)
  - `update_kappa.R` Eq. (27)
  - `update_beta.R` ESS targeting Eq. (32)
  - `update_phi.R` ICAR-aware ESS targeting Eq. (33)
  - `update_tau_phi.R` Eq. (31)
  - `ffbs_lambda.R` filtering Eq. (29) and smoothing Eq. (34)
  - `ffbs_omega.R` filtering Eq. (30) and smoothing Eq. (35)–(37)
  - `update_r.R` MH on log r targeting Appendix A.10
  - `update_discounts.R` grid SIR for gamma and delta (Section 8.6)
- `scripts/`
  - `02_fit_one_replicate.R` fit one replicate file
  - `03_fit_all_replicates.R` fit all scenarios and replicates in `data_sim/`

## Dependencies

- R base packages
- `Matrix` (required)

Install if needed:
```r
install.packages("Matrix")
```

## How to run

From your project root (the same directory that contains `data_sim/`):

### Fit one replicate
```bash
Rscript scripts/02_fit_one_replicate.R \
  --in  data_sim/S1_simple/rep0001.rds \
  --out mcmc_out/S1_simple/rep0001_fit.rds \
  --n_iter 2000 --burn 1000 --thin 1 --seed 123 \
  --store_gamma 1 --store_delta 1
```

### Fit everything
```bash
Rscript scripts/03_fit_all_replicates.R \
  --data_dir data_sim \
  --out_dir  mcmc_out \
  --n_iter 2000 --burn 1000 --thin 1 --seed0 1000
```

A `fit_summary.csv` with basic diagnostics is written to `mcmc_out/`.

## What the output contains

Each `*_fit.rds` is a list:
- `draws`: stored posterior draws (after burn-in and thinning)
- `last_state`: final MCMC state (useful for restarts)

By default, the driver scripts store:
- `beta0`, `beta`, `tau_phi`
- `gamma` (per coarsest region)
- `delta` (per internal node, flattened; see `draws$delta_nodes` for mapping)

You can turn on large objects by passing:
- `--store_phi 1` to store `phi`
- `--store_lambda 1` to store the full `lambda_tilde` path
- `--store_r 1` to store `r1`

## Hyperparameters and tuning

Defaults are set in `R/hyper.R` via `default_hyper(p)`.  
In `scripts/02_fit_one_replicate.R`, you can override:
- regression prior `(m0,s0,m_beta,V_beta)`
- ICAR precision prior `(a_phi,b_phi)`
- lambda initial prior `(a0,b0)`
- split concentration `alpha_dir`
- dispersion prior `(a_r,b_r)`
- discount priors `(a_gamma,b_gamma,a_delta,b_delta)`

Tuning knobs (in the `mcmc` list) include:
- `prop_sd_r`: proposal SD for the log r MH step
- `update_discounts_every`: how often to update discount factors
- `gamma_grid`, `delta_grid`: grids for SIR

## Notes on correctness vs speed

This is a “correctness-first” reference implementation:
- Everything is written in R.
- Sparse ops use `Matrix` where helpful, but `B` is dense (as discussed previously), so some products become dense.
- Once you are satisfied with correctness, the most impactful accelerations usually are:
  - moving the coarsest likelihood evaluations into Rcpp
  - exploiting sparsity in `H` more aggressively for larger `n1`
  - parallelizing across replicates and scenarios
