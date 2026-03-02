# loglik.R
# Coarsest Poisson log-likelihood (up to additive constants) used in ESS blocks.

loglik_coarsest <- function(beta0, beta, phi, lambda_tilde_0T, kappa, y1, e1, x1, Xmat = NULL) {
  TT <- nrow(y1)
  n1 <- ncol(y1)
  p  <- dim(x1)[3]

  assert_true(nrow(e1) == TT && ncol(e1) == n1, "e1 must be T x n1.")
  assert_true(all(dim(x1)[1:2] == c(TT, n1)), "x1 must be T x n1 x p.")
  assert_true(length(phi) == n1, "phi length mismatch.")
  assert_true(length(beta) == p, "beta length mismatch.")
  assert_true(all(dim(kappa) == c(TT, n1)), "kappa must be T x n1.")
  assert_true(all(dim(lambda_tilde_0T) == c(TT + 1, n1)), "lambda_tilde_0T must be (T+1) x n1.")

  lam_t <- lambda_tilde_0T[-1, , drop = FALSE]
  y_vec <- as.numeric(y1)
  log_cmult_vec <- as.numeric(
    log(pmax(e1, 1e-300)) +
      log(pmax(kappa, 1e-300)) +
      log(pmax(lam_t, 1e-300))
  )

  if (p > 0) {
    if (is.null(Xmat)) {
      Xmat <- matrix(x1, nrow = TT * n1, ncol = p)
    } else {
      assert_true(nrow(Xmat) == TT * n1 && ncol(Xmat) == p, "Xmat dim mismatch.")
    }
    xb_vec <- as.numeric(Xmat %*% beta)
  } else {
    xb_vec <- 0
  }

  phi_rep <- rep.int(as.numeric(phi), each = TT)
  eta_vec <- beta0 + xb_vec + phi_rep
  ce <- clamp(eta_vec)

  log_mu_vec <- log_cmult_vec + ce
  mu_vec <- exp(clamp(log_mu_vec))

  val <- sum(y_vec * ce - mu_vec)
  if (!is.finite(val)) return(-Inf)
  val
}
