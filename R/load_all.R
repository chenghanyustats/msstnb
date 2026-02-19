
# load_all.R
# Source all sampler components (assumes working directory set to project root).

src_msstnb_algo2 <- function(root = ".") {
  rdir <- file.path(root, "R")
  files <- c("utils.R","ess.R","loglik.R",
             "update_kappa.R","update_beta.R","update_phi.R","update_tau_phi.R",
             "ffbs_lambda.R","ffbs_omega.R",
             "update_r.R","update_discounts.R",
             "hyper.R","init_state.R","sampler_algo2.R")
  for (f in files) source(file.path(rdir, f))
  invisible(TRUE)
}
