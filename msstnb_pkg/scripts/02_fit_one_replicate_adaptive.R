# 02_fit_one_replicate_interactive.R
# Run MSSTNB Algorithm 2 for one replicate from an interactive R session.
#
# How to use
#   1. Open this file in RStudio
#   2. Edit the user settings section below
#   3. Run the whole script


###############################################################################
# User settings
###############################################################################

# Project root folder that contains the R folder
project_root <- "."

# Path to the msstnb_pkg folder (code lives under msstnb_pkg/R)
# Default assumes your folder structure is:
#   project_root/
#     data_sim/
#     mcmc_out/
#     msstnb_pkg/
pkg_root <- file.path(project_root, "msstnb_pkg")

# Choose one of the two input methods

# Method A: specify scenario name and replicate id
scenario_name <- "S0"
rep_id <- 1L

# Method B: directly specify the full input file path
# If you set in_file manually, it will override scenario_name and rep_id.
in_file <- NULL

# Output file path
# If NULL, the script will save next to the input file with suffix _fit.rds
out_file <- NULL

# MCMC controls
n_iter <- 2000L
burn <- 1000L
thin <- 1L
seed <- 1L
verbose <- TRUE

# Storage controls
store_phi <- FALSE
store_lambda <- FALSE
store_gamma <- TRUE
store_delta <- TRUE
store_r <- FALSE

# Sampler tuning
prop_sd_r <- 0.2
update_discounts_every <- 1L

# Adaptive Metropolis-Hastings for r (global prop_sd_r)
adapt_r <- TRUE
adapt_r_window <- 100L
adapt_r_target <- 0.30
adapt_r_gamma0 <- 0.05
adapt_r_t0 <- 10
adapt_r_power <- 0.5
adapt_r_min_sd <- 0.01
adapt_r_max_sd <- 2.0

# By default, adapt only during burn-in. Set to NULL to use burn.
adapt_r_until <- NULL

# Optional discount grids
# Set to NULL to use the defaults inside msstnb_algo2_mcmc
# gamma_grid <- seq(0.05, 0.999, length.out = 200L)
# delta_grid <- seq(0.05, 0.999, length.out = 200L)
gamma_grid <- NULL
delta_grid <- NULL

# Optional split prior override
# q0_splits is used in ffbs_omega and update_delta_grid
q0_splits <- NULL

# Optional time index for discount updates
# Use NULL for all times
# Example: discount_time_idx <- 10:TT
discount_time_idx <- NULL


###############################################################################
# Helper functions
###############################################################################

build_rep_paths <- function(project_root, scenario_name, rep_id) {
    rep_tag <- sprintf("rep%04d.rds", as.integer(rep_id))
    in_path <- file.path(project_root, "data_sim", scenario_name, rep_tag)

    out_tag <- sprintf("rep%04d_fit.rds", as.integer(rep_id))
    out_path <- file.path(project_root, "mcmc_out", scenario_name, out_tag)

    list(in_file = in_path, out_file = out_path)
}

source_msstnb_algo2 <- function(project_root) {
    r_dir <- file.path(project_root, "R")

    files <- c(
        "cache.R",
        "utils.R",
        "tree_utils.R",
        "ess.R",
        "loglik.R",
        "update_kappa.R",
        "update_beta.R",
        "update_phi.R",
        "update_tau_phi.R",
        "ffbs_lambda.R",
        "ffbs_omega.R",
        "update_r.R",
        "update_discounts.R",
        "hyper.R",
        "init_state.R",
        "sampler_algo2_reordered.R"
    )

    missing <- files[!file.exists(file.path(r_dir, files))]
    if (length(missing) > 0) {
        stop(
            paste(
                "Missing files under", shQuote(r_dir), ":\n",
                paste(missing, collapse = "\n")
            )
        )
    }

    for (f in files) {
        source(file.path(r_dir, f))
    }

    if (!exists("msstnb_algo2_mcmc", mode = "function")) {
        stop("msstnb_algo2_mcmc was not found after sourcing files")
    }

    invisible(TRUE)
}


###############################################################################
# Main script
###############################################################################

# 1. Resolve input and output paths
if (is.null(in_file)) {
    paths <- build_rep_paths(project_root, scenario_name, rep_id)
    in_file <- paths$in_file
    if (is.null(out_file)) {
        out_file <- paths$out_file
    }
} else {
    if (is.null(out_file)) {
        out_file <- sub("\\.rds$", "_fit.rds", in_file)
    }
}

if (!file.exists(in_file)) {
    stop(paste("Input file does not exist:", shQuote(in_file)))
}

# 2. Source the sampler and all dependencies
source_msstnb_algo2(pkg_root)

# 3. Load one replicate
ds <- readRDS(in_file)
if (is.null(ds$data)) {
    stop("The input rds does not contain ds$data")
}

data <- ds$data

# 4. Build hyperparameters
# data$x is TT by n1 by p
p <- dim(data$x)[3]
hyper <- default_hyper(p = p)

# 5. Build MCMC list
if (is.null(adapt_r_until)) {
    adapt_r_until <- as.integer(burn)
} else {
    adapt_r_until <- as.integer(adapt_r_until)
}

mcmc <- list(
    n_iter = as.integer(n_iter),
    burn = as.integer(burn),
    thin = as.integer(thin),
    prop_sd_r = prop_sd_r,
    adapt_r = isTRUE(adapt_r),
    adapt_r_window = as.integer(adapt_r_window),
    adapt_r_target = as.numeric(adapt_r_target),
    adapt_r_gamma0 = as.numeric(adapt_r_gamma0),
    adapt_r_t0 = as.numeric(adapt_r_t0),
    adapt_r_power = as.numeric(adapt_r_power),
    adapt_r_min_sd = as.numeric(adapt_r_min_sd),
    adapt_r_max_sd = as.numeric(adapt_r_max_sd),
    adapt_r_until = as.integer(adapt_r_until),
    update_discounts_every = as.integer(update_discounts_every),
    store_phi = isTRUE(store_phi),
    store_lambda = isTRUE(store_lambda),
    store_gamma = isTRUE(store_gamma),
    store_delta = isTRUE(store_delta),
    store_r = isTRUE(store_r),
    gamma_grid = gamma_grid,
    delta_grid = delta_grid,
    q0_splits = q0_splits,
    discount_time_idx = discount_time_idx
)

# 6. Run sampler
cat("Running msstnb_algo2_mcmc\n")
cat("Input :", in_file, "\n")
cat("Output:", out_file, "\n")

start_time <- Sys.time()
fit <- msstnb_algo2_mcmc(
    data = data,
    hyper = hyper,
    mcmc = mcmc,
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
)
end_time <- Sys.time()

cat("Finished\n")
cat("Elapsed:", as.numeric(difftime(end_time, start_time, units = "secs")), "seconds\n")

# 7. Save results
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

saveRDS(fit, out_file)
cat("Saved fit to:", out_file, "\n")

# 8. Quick diagnostics summary
if (!is.null(fit$diag) && is.list(fit$diag)) {
    if (!is.null(fit$diag$accept_r)) {
        cat("Mean accept fraction across regions (all iters): ", mean(fit$diag$accept_r), "\n", sep = "")
        if (burn < length(fit$diag$accept_r)) {
            cat("Mean accept fraction across regions (post burn): ", mean(fit$diag$accept_r[(burn + 1L):length(fit$diag$accept_r)]), "\n", sep = "")
        }
    }
    if (!is.null(fit$diag$accept_r_by_region_postburn)) {
        x <- as.numeric(fit$diag$accept_r_by_region_postburn)
        cat("Acceptance by region post burn summary:\n")
        cat("  min : ", min(x), "\n", sep = "")
        cat("  mean: ", mean(x), "\n", sep = "")
        cat("  max : ", max(x), "\n", sep = "")
    }
    if (!is.null(fit$diag$prop_sd_r)) {
        cat("prop sd r last value: ", tail(fit$diag$prop_sd_r, 1), "\n", sep = "")
    }
}
