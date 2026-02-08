# generate_dataset.R
# Main dataset generator for MSSTNB simulation studies
#
# Expected that the driver script has already sourced these files
# R/utils_io.R
# R/make_graph.R
# R/sim_icar.R
# R/make_tree.R
# R/sim_lambda_tilde.R
# R/sim_splits.R

compute_xbeta <- function(x, beta) {
    # x is TT by n1 by p, beta is p
    d <- dim(x)
    TT <- d[1]; n1 <- d[2]; p <- d[3]
    assert_true(length(beta) == p, "beta length must match p.")
    out <- matrix(0, nrow = TT, ncol = n1)
    for (k in seq_len(p)) out <- out + x[, , k] * beta[k]
    out
}

allocate_counts <- function(tree, y1, omega) {
  # y1 is TT by n1, omega is list omega[[l]][[j]] matrices TT by K
  TT <- nrow(y1)
  y_levels <- vector("list", tree$L)
  y_levels[[1]] <- y1

  if (tree$L == 1) return(y_levels)

  for (l in seq_len(tree$L - 1)) {
    n_child <- tree$n[l + 1]
    y_next <- matrix(0L, nrow = TT, ncol = n_child)
    # y_next <- matrix(0, nrow = T, ncol = n_child)

    for (j in seq_len(tree$n[l])) {
      ch <- tree$children[[l]][[j]]
      K <- length(ch)
      if (K == 0) next
      omega_mat <- omega[[l]][[j]]
      assert_true(is.matrix(omega_mat), "omega_mat must be a matrix.")
      for (t in seq_len(TT)) {
        y_next[t, ch] <- rmultinom1(y_levels[[l]][t, j], omega_mat[t, ])
      }
      # for (t in seq_len(T)) {
      #   size_tj <- y_levels[[l]][t, j]
      #   assert_true(is.finite(size_tj) && size_tj >= 0, "Parent count must be finite and nonnegative.")
      #   y_next[t, ch] <- rmultinom1(size_tj, omega_mat[t, ])
      # }

    }
    y_levels[[l + 1]] <- y_next
  }
  y_levels
}

generate_dataset <- function(scenario_row, rep_id = 1L) {
  # scenario_row can be a one row data.frame or a named list
  scen <- scenario_row
  if (is.data.frame(scen)) scen <- as.list(scen[1, , drop = TRUE])

  required <- c("scenario_id", "TT", "n1", "L", "branching", "refine_prob",
                "graph_type", "tau_phi", "r", "alpha_dir", "gamma_lambda",
                "lambda_a0", "lambda_b0",
                "p", "beta0", "beta_vals", "x_sd", "loge_mean", "loge_sd",
                "zi_prob", "comp_mode", "cp_time", "cp_strength", "q_drift_sd",
                "q_conc",
                "seed_base")
  missing <- setdiff(required, names(scen))
  assert_true(length(missing) == 0, paste("Missing scenario fields:", paste(missing, collapse = ", ")))

  TT <- as.integer(scen$TT)
  n1 <- as.integer(scen$n1)
  L <- as.integer(scen$L)
  branching <- parse_int_vec(scen$branching)
  refine_prob <- as.numeric(scen$refine_prob)
  graph_type <- as.character(scen$graph_type)
  tau_phi <- as.numeric(scen$tau_phi)
  r <- as.numeric(scen$r)
  alpha_dir <- as.numeric(scen$alpha_dir)
  gamma_lambda <- as.numeric(scen$gamma_lambda)
  lambda_a0 <- as.numeric(scen$lambda_a0)
  lambda_b0 <- as.numeric(scen$lambda_b0)
  p <- as.integer(scen$p)
  beta0 <- as.numeric(scen$beta0)
  beta_vals <- parse_num_vec(scen$beta_vals)
  x_sd <- as.numeric(scen$x_sd)
  loge_mean <- as.numeric(scen$loge_mean)
  loge_sd <- as.numeric(scen$loge_sd)
  zi_prob <- as.numeric(scen$zi_prob)
  comp_mode <- as.character(scen$comp_mode)
  cp_time <- as.integer(scen$cp_time)
  cp_strength <- as.numeric(scen$cp_strength)
  q_drift_sd <- as.numeric(scen$q_drift_sd)
  q_conc <- as.numeric(scen$q_conc)
  seed_base <- as.integer(scen$seed_base)

  assert_true(length(beta_vals) == p, "beta_vals must have length p.")
  assert_true(zi_prob >= 0 && zi_prob < 1, "zi_prob must be in [0,1).")

  set_seed(seed_base + as.integer(rep_id))

  # Graph and ICAR
  graph <- make_graph(n1, graph_type = graph_type)
  phi <- sim_phi_icar(graph$H, tau_phi = tau_phi)

  # Exposures and covariates
  e <- matrix(exp(rnorm(TT * n1, mean = loge_mean, sd = loge_sd)), nrow = TT, ncol = n1)
  x <- array(rnorm(TT * n1 * p, mean = 0, sd = x_sd), dim = c(TT, n1, p))
  xbeta <- compute_xbeta(x, beta_vals)

  # Dynamic residual risk lambda tilde
  lambda_tilde <- sim_lambda_tilde(
    TT, n1,
    gamma = gamma_lambda,
    a0 = lambda_a0, b0 = lambda_b0,
    mode = "discount",
    seed = seed_base + 91000 + rep_id
  )

  # Overdispersion parameters and kappa
  r1 <- rep(r, n1)
  kappa <- matrix(NA_real_, nrow = TT, ncol = n1)
  for (j in seq_len(n1)) kappa[, j] <- rgamma(TT, shape = r1[j], rate = r1[j])

  # Coarsest mean and counts
  linpred <- beta0 + xbeta + matrix(phi, nrow = TT, ncol = n1, byrow = TRUE)
  mu <- e * exp(linpred) * lambda_tilde
  y1 <- matrix(rpois(TT * n1, lambda = as.numeric(mu * kappa)), nrow = TT, ncol = n1)

  # Optional zero inflation misspecification
  if (zi_prob > 0) {
    zmask <- matrix(rbinom(TT * n1, size = 1, prob = zi_prob) == 1,
                    nrow = TT, ncol = n1)
    y1[zmask] <- 0L
  }

  # Tree and splits
  tree <- make_tree(
    n1, L,
    branching = branching,
    refine_prob = refine_prob,
    seed = seed_base + 42000 + rep_id
  )
  splits <- sim_splits(
    tree, TT,
    alpha_dir = alpha_dir,
    comp_mode = comp_mode,
    cp_time = cp_time,
    cp_strength = cp_strength,
    q_drift_sd = q_drift_sd,
    q_conc = q_conc,
    seed = seed_base + 77000 + rep_id
  )

  y_levels <- allocate_counts(tree, y1, splits$omega)

  list(
    data = list(
      y_levels = y_levels,
      e = e,
      x = x,
      tree = tree,
      graph = list(edges = graph$edges, graph_type = graph$graph_type),
      level1 = list(n1 = n1, TT = TT, p = p)
    ),
    truth = list(
      beta0 = beta0,
      beta = beta_vals,
      phi = phi,
      tau_phi = tau_phi,
      r1 = r1,
      kappa = kappa,
      lambda_tilde = lambda_tilde,
      omega = splits$omega,
      q = splits$q,
      alpha_dir = alpha_dir,
      scenario = scen
    )
  )
}

