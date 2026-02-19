
#!/usr/bin/env Rscript
# 02_fit_one_replicate.R
# Fit MSSTNB with Algorithm 2 to a single simulated dataset (.rds created by 01_generate_all_data.R).

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value after", flag))
  args[hit + 1]
}

in_file  <- get_arg("--in",  default = NULL)
out_file <- get_arg("--out", default = NULL)

if (is.null(in_file) || is.null(out_file)) {
  cat("Usage:\n")
  cat("  Rscript scripts/02_fit_one_replicate.R --in data_sim/<scenario>/rep0001.rds --out mcmc_out/<scenario>/rep0001_fit.rds\n")
  cat("Optional:\n")
  cat("  --n_iter 2000 --burn 1000 --thin 1 --seed 123 --store_phi 0/1 --store_lambda 0/1 --store_r 0/1 --store_gamma 0/1 --store_delta 0/1\n")
  quit(status = 1)
}

n_iter <- as.integer(get_arg("--n_iter", 2000))
burn   <- as.integer(get_arg("--burn",   floor(n_iter/2)))
thin   <- as.integer(get_arg("--thin",   1))
seed   <- as.integer(get_arg("--seed",   1))

store_phi    <- as.integer(get_arg("--store_phi", 0)) == 1
store_lambda <- as.integer(get_arg("--store_lambda", 0)) == 1
store_r      <- as.integer(get_arg("--store_r", 0)) == 1
store_gamma  <- as.integer(get_arg("--store_gamma", 1)) == 1
store_delta  <- as.integer(get_arg("--store_delta", 1)) == 1

# Source sampler
source("R/load_all.R")
src_msstnb_algo2(".")

ds <- readRDS(in_file)
data <- ds$data

# Hyperparameters: start with defaults; can override here
p <- dim(data$x)[3]
hyper <- default_hyper(p = p)

mcmc <- list(
  n_iter = n_iter,
  burn = burn,
  thin = thin,
  prop_sd_r = 0.2,
  update_discounts_every = 1,
  store_phi = store_phi,
  store_lambda = store_lambda,
  store_r = store_r,
  store_gamma = store_gamma,
  store_delta = store_delta
)

fit <- msstnb_algo2_mcmc(data = data, hyper = hyper, mcmc = mcmc, seed = seed, verbose = TRUE)

# Ensure output directory exists
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(fit, out_file)
cat("Saved fit to:", out_file, "\n")
