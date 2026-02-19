
#!/usr/bin/env Rscript
# 03_fit_all_replicates.R
# Fit MSSTNB with Algorithm 2 across all scenarios and replicates in data_sim/.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value after", flag))
  args[hit + 1]
}

data_dir <- get_arg("--data_dir", "data_sim")
out_dir  <- get_arg("--out_dir",  "mcmc_out")
n_iter   <- as.integer(get_arg("--n_iter", 2000))
burn     <- as.integer(get_arg("--burn", floor(n_iter/2)))
thin     <- as.integer(get_arg("--thin", 1))
seed0    <- as.integer(get_arg("--seed0", 1))

if (!dir.exists(data_dir)) stop(paste("Missing data_dir:", data_dir))

source("R/load_all.R")
src_msstnb_algo2(".")

scenarios <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
if (length(scenarios) == 0) stop("No scenarios found in data_dir.")

summary_rows <- list()
k <- 1

for (scen in scenarios) {
  scen_path <- file.path(data_dir, scen)
  rep_files <- list.files(scen_path, pattern = "^rep[0-9]{4}\\.rds$", full.names = TRUE)
  if (length(rep_files) == 0) next

  cat("\nScenario:", scen, "reps:", length(rep_files), "\n")
  out_scen_dir <- file.path(out_dir, scen)
  if (!dir.exists(out_scen_dir)) dir.create(out_scen_dir, recursive = TRUE)

  for (rf in rep_files) {
    rep_id <- sub("^rep([0-9]{4})\\.rds$", "\\1", basename(rf))
    out_file <- file.path(out_scen_dir, paste0("rep", rep_id, "_fit.rds"))
    if (file.exists(out_file)) {
      cat("  skip existing:", out_file, "\n")
      next
    }
    ds <- readRDS(rf)
    data <- ds$data
    p <- dim(data$x)[3]
    hyper <- default_hyper(p = p)

    mcmc <- list(
      n_iter = n_iter,
      burn = burn,
      thin = thin,
      prop_sd_r = 0.2,
      update_discounts_every = 1,
      store_phi = FALSE,
      store_lambda = FALSE,
      store_r = FALSE,
      store_gamma = TRUE,
      store_delta = TRUE
    )

    seed <- seed0 + as.integer(rep_id)
    fit <- msstnb_algo2_mcmc(data = data, hyper = hyper, mcmc = mcmc, seed = seed, verbose = FALSE)
    saveRDS(fit, out_file)

    # collect summary
    acc_r <- mean(fit$draws$diagnostics$acc_r)
    ess_b <- mean(fit$draws$diagnostics$ess_beta_steps)
    ess_p <- mean(fit$draws$diagnostics$ess_phi_steps)
    summary_rows[[k]] <- data.frame(scenario = scen, rep = rep_id,
                                    acc_r = acc_r, ess_beta = ess_b, ess_phi = ess_p,
                                    out_file = out_file, stringsAsFactors = FALSE)
    k <- k + 1
    cat("  wrote:", out_file, "\n")
  }
}

if (length(summary_rows) > 0) {
  sumdf <- do.call(rbind, summary_rows)
  sum_file <- file.path(out_dir, "fit_summary.csv")
  write.csv(sumdf, sum_file, row.names = FALSE)
  cat("\nWrote summary:", sum_file, "\n")
}
