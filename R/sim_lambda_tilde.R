# sim_lambda_tilde.R
# Simulate dynamic residual risk trajectories at level 1
# Aktekin, Soyer, Xu (2013)
# Irie, Aktekin (2025)
# Chen et al. (2018), JASA
sim_lambda_tilde <- function(TT, n1, gamma = 0.98, a0 = 20, b0 = 20,
                            mode = c("discount", "log_rw"),
                            cp_time = 0L, cp_strength = 1.0,
                            rw_sd = NULL,
                            seed = NULL) {
  set_seed(seed)
  mode <- match.arg(mode)

  assert_true(TT >= 1, "T must be at least 1.")
  assert_true(n1 >= 1, "n1 must be at least 1.")
  assert_true(gamma > 0 && gamma <= 1, "gamma must be in (0,1].")
  assert_true(a0 > 0 && b0 > 0, "a0 and b0 must be positive.")
  assert_true(cp_strength > 0, "cp_strength must be positive.")

  lam <- matrix(NA_real_, nrow = TT, ncol = n1)
  # Prior mean near 1 is recommended in the manuscript
  lam[1, ] <- rgamma(n1, shape = a0, rate = b0)

  if (mode == "discount") {
    if (gamma == 1) {
      for (t in 2:TT) lam[t, ] <- lam[t - 1, ]
    } else {
      a1 <- gamma * a0
      a2 <- (1 - gamma) * a0
      assert_true(a1 > 0 && a2 > 0, "a0 must be large enough for discounting.")
      for (t in 2:TT) {
        eta <- rbeta(n1, shape1 = a1, shape2 = a2)
        lam[t, ] <- lam[t - 1, ] * (eta / gamma)
      }
    }
  }

  if (mode == "log_rw") {
    if (is.null(rw_sd)) rw_sd <- 0.25 * sqrt(max(1e-8, 1 - gamma))
    for (t in 2:TT) {
      eps <- rnorm(n1, mean = -0.5 * rw_sd^2, sd = rw_sd)
      lam[t, ] <- lam[t - 1, ] * exp(eps)
    }
  }

  if (!is.null(cp_time) && cp_time > 0 && cp_time < TT) {
    lam[(cp_time + 1):TT, ] <- lam[(cp_time + 1):TT, , drop = FALSE] * cp_strength
  }

  lam
}

