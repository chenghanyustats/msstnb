
# init_state.R
# Initialize parameters for the blocked Gibbs sampler (Algorithm 2).

init_state <- function(data, hyper, seed = NULL) {
  set_seed(seed)
  assert_true(is.list(data), "data must be a list.")

  y_levels <- data$y_levels
  e1 <- data$e
  x1 <- data$x
  tree <- data$tree
  edges <- data$graph$edges

  assert_true(!is.null(y_levels) && is.list(y_levels), "data$y_levels must be a list.")
  L <- tree$L
  assert_true(length(y_levels) == L, "y_levels length must equal tree$L.")
  TT <- nrow(y_levels[[1]])
  n1 <- ncol(y_levels[[1]])
  p  <- dim(x1)[3]

  # Build H and basis
  H <- build_H_from_edges(edges, n1 = n1)
  Bb <- build_B_basis(n1)
  B <- Bb$B

  beta0 <- 0
  beta  <- rep(0, p)

  phi <- rep(0, n1)
  tau_phi <- 1

  # initial lambda: constant near 1
  lambda_tilde_0T <- matrix(1, nrow = TT + 1, ncol = n1)

  # initial kappa: 1
  kappa <- matrix(1, nrow = TT, ncol = n1)

  # initial r
  r1 <- rep(10, n1)

  # initial discounts
  gamma <- rep(hyper$discount$gamma_init, n1)
  delta <- init_delta_splits(tree, delta0 = hyper$discount$delta_init)

  # initial omega: smoothed empirical proportions at each time
  omega <- vector("list", L - 1)
  for (l in seq_len(L - 1)) {
    n_parent <- tree$n[l]
    omega_l <- vector("list", n_parent)
    for (j in seq_len(n_parent)) {
      ch <- tree$children[[l]][[j]]
      K <- length(ch)
      if (K <= 1) { omega_l[[j]] <- NULL; next }
      y_child <- y_levels[[l + 1]][, ch, drop = FALSE] # T x K
      om <- matrix(NA_real_, nrow = TT + 1, ncol = K)
      for (t in seq_len(TT)) {
        v <- y_child[t, ] + 1
        om[t + 1, ] <- v / sum(v)
      }
      om[1, ] <- om[2, ] # t=0 same as t=1
      omega_l[[j]] <- om
    }
    omega[[l]] <- omega_l
  }

  list(
    beta0 = beta0,
    beta  = beta,
    phi   = phi,
    tau_phi = tau_phi,
    lambda_tilde_0T = lambda_tilde_0T,
    kappa = kappa,
    r1 = r1,
    gamma = gamma,
    delta = delta,
    omega = omega,
    # cached
    H = H,
    B = B
  )
}
