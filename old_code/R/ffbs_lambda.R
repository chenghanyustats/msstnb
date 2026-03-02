
# ffbs_lambda.R
# Forward filtering (29) and backward sampling (34) for lambda_tilde trajectories.

filter_lambda <- function(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0) {
    TT <- nrow(y1); n1 <- ncol(y1); p <- dim(x1)[3]
    assert_true(length(gamma) == n1, "gamma must have length n1.")
    assert_true(all(gamma > 0 & gamma <= 1), "gamma must be in (0, 1].")
    assert_true(length(a0) %in% c(1, n1), "a0 must be scalar or length n1.")
    assert_true(length(b0) %in% c(1, n1), "b0 must be scalar or length n1.")
    if (length(a0) == 1) a0 <- rep(a0, n1)
    if (length(b0) == 1) b0 <- rep(b0, n1)

    # compute xi_t = e * kappa * exp(beta0 + xbeta + phi)
    eta <- matrix(beta0, nrow = TT, ncol = n1)

    if (p > 0) {
        for (j in seq_len(n1)) eta[, j] <- eta[, j] +
            as.numeric(x1[, j, , drop = TRUE] %*% beta)
    }
    eta <- eta + matrix(rep(phi, each = TT), nrow = TT, ncol = n1)
    xi <- e1 * kappa * exp(clamp(eta))

    a <- matrix(NA_real_, nrow = TT + 1, ncol = n1)
    b <- matrix(NA_real_, nrow = TT + 1, ncol = n1)
    a[1, ] <- a0
    b[1, ] <- b0
    for (t in seq_len(TT)) {
        a[t + 1, ] <- gamma * a[t, ] + y1[t, ]
        b[t + 1, ] <- gamma * b[t, ] + xi[t, ]
    }
    list(a = a, b = b, xi = xi)
}

backward_sample_lambda <- function(a, b, gamma) {
    TT <- nrow(a) - 1
    n1 <- ncol(a)
    lam <- matrix(NA_real_, nrow = TT + 1, ncol = n1)
    # sample lambda_T
    lam[TT + 1, ] <- rgamma(n1, shape = pmax(a[TT + 1, ], 1e-8),
                            rate = pmax(b[TT + 1, ], 1e-12))
    # backward
    for (t in TT:1) {
        for (j in seq_len(n1)) {
            if (gamma[j] >= 1) {
                lam[t, j] <- lam[t + 1, j]
        } else {
            sh <- (1 - gamma[j]) * a[t, j]
            sh <- pmax(sh, 1e-8)
            rate <- pmax(b[t, j], 1e-12)
            U <- rgamma(1, shape = sh, rate = rate)
            lam[t, j] <- U + gamma[j] * lam[t + 1, j]
        }
      }
    }
    lam
}

ffbs_lambda <- function(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0,
                        seed = NULL) {
    set_seed(seed)
    filt <- filter_lambda(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0)
    lam <- backward_sample_lambda(filt$a, filt$b, gamma)
    list(lambda_tilde_0T = lam, a = filt$a, b = filt$b, xi = filt$xi)
}
