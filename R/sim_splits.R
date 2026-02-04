\
# sim_splits.R
# Simulate baseline compositions q and split probabilities omega down the tree

make_initial_q <- function(K, q_conc = 30) {
  # Random q near uniform: Dirichlet(q_conc * 1/K)
  rdirichlet1(rep(q_conc / K, K))
}

perturb_simplex <- function(q, sd = 0.05) {
  # Multiplicative log normal perturbation, stays in simplex after renormalization
  g <- rnorm(length(q), mean = 0, sd = sd)
  z <- q * exp(g)
  z / sum(z)
}

spike_simplex <- function(q, strength = 3.0) {
  K <- length(q)
  kstar <- sample.int(K, 1)
  z <- q
  z[kstar] <- z[kstar] * strength
  z / sum(z)
}

sim_splits <- function(tree, T,
                       alpha_dir = 50,
                       comp_mode = c("static", "drift", "changepoint"),
                       cp_time = 0L,
                       cp_strength = 3.0,
                       q_drift_sd = 0.05,
                       q_conc = 30,
                       seed = NULL) {
  set_seed(seed)
  comp_mode <- match.arg(comp_mode)

  assert_true(alpha_dir > 0, "alpha_dir must be positive.")
  assert_true(T >= 1, "T must be at least 1.")

  L <- tree$L
  if (L == 1) {
    return(list(omega = list(), q = list(), alpha_dir = alpha_dir))
  }

  omega <- vector("list", L - 1)
  q_time <- vector("list", L - 1)

  for (l in seq_len(L - 1)) {
    n_parent <- tree$n[l]
    omega_l <- vector("list", n_parent)
    q_l <- vector("list", n_parent)

    for (j in seq_len(n_parent)) {
      ch <- tree$children[[l]][[j]]
      K <- length(ch)
      if (K == 0) {
        omega_l[[j]] <- NULL
        q_l[[j]] <- NULL
        next
      }

      qmat <- matrix(NA_real_, nrow = T, ncol = K)
      qmat[1, ] <- make_initial_q(K, q_conc = q_conc)

      if (comp_mode == "static") {
        for (t in 2:T) qmat[t, ] <- qmat[1, ]
      }

      if (comp_mode == "drift") {
        for (t in 2:T) qmat[t, ] <- perturb_simplex(qmat[t - 1, ], sd = q_drift_sd)
      }

      if (comp_mode == "changepoint") {
        for (t in 2:T) qmat[t, ] <- qmat[1, ]
        if (!is.null(cp_time) && cp_time > 0 && cp_time < T) {
          q_new <- spike_simplex(qmat[1, ], strength = cp_strength)
          qmat[(cp_time + 1):T, ] <- matrix(q_new, nrow = T - cp_time, ncol = K, byrow = TRUE)
        }
      }

      # Omega draws
      omega_mat <- matrix(NA_real_, nrow = T, ncol = K)
      for (t in seq_len(T)) {
        omega_mat[t, ] <- rdirichlet1(alpha_dir * qmat[t, ])
      }

      omega_l[[j]] <- omega_mat
      q_l[[j]] <- qmat
    }

    omega[[l]] <- omega_l
    q_time[[l]] <- q_l
  }

  list(omega = omega, q = q_time, alpha_dir = alpha_dir)
}

