
# ess.R
# Elliptical slice sampling (Murray et al. 2010) as described in Algorithm 1.

# out <- ess(u0, loglik_fn = loglik_fn, draw_nu_fn = draw_nu,
#            m = rep(0, n1 - 1))
ess <- function(z, loglik_fn, draw_nu_fn, m = NULL, max_steps = 10000L) {
    # z: current vector
    # loglik_fn: function(z) returning scalar log-likelihood (up to constant)
    # draw_nu_fn: function() returning a draw nu ~ N(0, Sigma) in same dimension as z
    # m: prior mean; if NULL, use 0
    z <- as.numeric(z)
    d <- length(z)
    if (is.null(m)) m <- rep(0, d) else m <- as.numeric(m)

    nu <- as.numeric(draw_nu_fn())
    assert_true(length(nu) == d, "draw_nu_fn returned wrong dimension.")

    loglik0 <- loglik_fn(z)
    assert_true(is.finite(loglik0), "Initial log-likelihood is not finite.")
    logy <- loglik0 + log(runif(1))

    twopi <- 2*pi
    theta <- runif(1, 0, twopi)
    thetamin <- theta - twopi
    thetamax <- theta

    for (step in seq_len(max_steps)) {
        z_prop <- m + (z - m) * cos(theta) + nu * sin(theta)
        lp <- loglik_fn(z_prop)
        if (is.finite(lp) && lp >= logy) {
            return(list(z = as.numeric(z_prop), n_steps = step))
        }
        # shrink bracket
        if (theta < 0) thetamin <- theta else thetamax <- theta
        theta <- runif(1, thetamin, thetamax)
    }
    stop("ESS exceeded max_steps without acceptance.
         Check log-likelihood for NaNs/Inf.", call. = FALSE)
}

# Draw nu ~ N(0, Sigma) given a triangular factor.
# If cholSigma is upper triangular from chol(Sigma), then Sigma = t(R) %*% R and we return t(R) %*% z.
# If cholSigma is lower triangular, then Sigma = L %*% t(L) and we return L %*% z.
draw_nu_cholSigma <- function(cholSigma, tol = 1e-12) {
    cholSigma <- as.matrix(cholSigma)
    n <- nrow(cholSigma)
    assert_true(ncol(cholSigma) == n, "cholSigma must be square.")

    z <- rnorm(n)

    upper_ok <- all(abs(cholSigma[lower.tri(cholSigma)]) < tol)
    lower_ok <- all(abs(cholSigma[upper.tri(cholSigma)]) < tol)

    if (upper_ok && !lower_ok) {
      return(as.numeric(t(cholSigma) %*% z))
    }
    if (lower_ok && !upper_ok) {
      return(as.numeric(cholSigma %*% z))
    }

    stop("cholSigma must be triangular (upper or lower).", call. = FALSE)
}
