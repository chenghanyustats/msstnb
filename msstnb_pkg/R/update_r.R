
# update_r.R
# Update dispersion r_{1j} for each coarsest region using 1D random-walk MH on log r.
# Targets Appendix A.10 kernel.

logpost_logr <- function(logr, kappa_col, a_r, b_r) {
    r <- exp(logr)
    if (!is.finite(r) || r <= 0) return(-Inf)
    # ensure kappa positive
    kappa_col <- pmax(kappa_col, 1e-300)
    Tt <- length(kappa_col)

    lp <- (a_r - 1) * logr - b_r * r
    lp <- lp + sum(r * logr - lgamma(r) + (r - 1) * log(kappa_col) -
                     r * kappa_col)
    # Jacobian for log-transform
    lp <- lp + logr
    lp
}

update_r_mh <- function(r1, kappa, a_r = 2, b_r = 2, prop_sd = 0.2) {
  TT <- nrow(kappa); n1 <- ncol(kappa)
  assert_true(length(r1) == n1, "r1 length mismatch.")
  assert_true(prop_sd > 0, "prop_sd must be positive.")
  out <- r1
  acc <- logical(n1)

  for (j in seq_len(n1)) {
      logr <- log(pmax(r1[j], 1e-8))
      lp  <- logpost_logr(logr, kappa[, j], a_r, b_r)
      logr_p <- logr + rnorm(1, 0, prop_sd)
      lp_p <- logpost_logr(logr_p, kappa[, j], a_r, b_r)
      if (is.finite(lp_p) && log(runif(1)) < (lp_p - lp)) {
          out[j] <- exp(logr_p)
          acc[j] <- TRUE
      } else {
          out[j] <- exp(logr)
          acc[j] <- FALSE
      }
  }
  list(r1 = out, accept = mean(acc), accept_vec = acc)
}
