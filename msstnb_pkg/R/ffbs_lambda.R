# ffbs_lambda.R
# Forward filtering and backward sampling for lambda_tilde trajectories.
#
# Optimizations:
# - optional direct xi input to avoid recomputing exp(eta)
# - vectorized backward sampling over regions (one rgamma call per time step)

filter_lambda_from_xi <- function(y1, xi, gamma, a0, b0) {
    TT <- nrow(y1)
    n1 <- ncol(y1)

    assert_true(all(dim(xi) == c(TT, n1)), "xi must be T x n1.")
    assert_true(length(gamma) == n1, "gamma must have length n1.")
    assert_true(all(gamma > 0 & gamma <= 1), "gamma must be in (0, 1].")

    if (length(a0) == 1L) {
        a0 <- rep(a0, n1)
    }
    if (length(b0) == 1L) {
        b0 <- rep(b0, n1)
    }
    assert_true(length(a0) == n1, "a0 must be scalar or length n1.")
    assert_true(length(b0) == n1, "b0 must be scalar or length n1.")

    a <- matrix(NA_real_, nrow = TT + 1L, ncol = n1)
    b <- matrix(NA_real_, nrow = TT + 1L, ncol = n1)
    a[1, ] <- a0
    b[1, ] <- b0

    for (t in seq_len(TT)) {
        a[t + 1L, ] <- gamma * a[t, ] + y1[t, ]
        b[t + 1L, ] <- gamma * b[t, ] + xi[t, ]
    }

    list(a = a, b = b)
}

backward_sample_lambda <- function(a, b, gamma) {
    TT <- nrow(a) - 1L
    n1 <- ncol(a)

    lam <- matrix(NA_real_, nrow = TT + 1L, ncol = n1)

    lam[TT + 1L, ] <- rgamma(
        n1,
        shape = pmax(a[TT + 1L, ], 1e-8),
        rate  = pmax(b[TT + 1L, ], 1e-12)
    )

    keep <- gamma >= 1

    for (t in TT:1) {
        if (any(keep)) {
            lam[t, keep] <- lam[t + 1L, keep]
        }

        idx <- !keep
        if (any(idx)) {
            sh <- (1 - gamma[idx]) * a[t, idx]
            sh <- pmax(sh, 1e-8)
            rt <- pmax(b[t, idx], 1e-12)
            U <- rgamma(sum(idx), shape = sh, rate = rt)
            lam[t, idx] <- U + gamma[idx] * lam[t + 1L, idx]
        }
    }

    lam
}

# Fast path: user supplies xi.
ffbs_lambda_from_xi <- function(y1, xi, gamma, a0, b0, seed = NULL) {
    set_seed(seed)
    filt <- filter_lambda_from_xi(y1 = y1, xi = xi, gamma = gamma,
                                  a0 = a0, b0 = b0)
    lam <- backward_sample_lambda(filt$a, filt$b, gamma = gamma)
    list(lambda_tilde_0T = lam, a = filt$a, b = filt$b)
}


# Backward-compatible wrapper (keeps previous signature).
filter_lambda <- function(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0, Xmat = NULL) {
    TT <- nrow(y1)
    n1 <- ncol(y1)
    p  <- dim(x1)[3]

    if (p > 0L && is.null(Xmat)) {
        Xmat <- prep_Xmat(x1)
    }

    eta <- compute_eta(beta0 = beta0, beta = beta, phi = phi, TT = TT, n1 = n1, Xmat = Xmat)
    xi <- compute_xi(e1 = e1, kappa = kappa, eta = eta)

    filt <- filter_lambda_from_xi(y1 = y1, xi = xi, gamma = gamma, a0 = a0, b0 = b0)
    filt$xi <- xi
    filt
}

ffbs_lambda <- function(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0,
                        seed = NULL, Xmat = NULL, xi = NULL) {
    set_seed(seed)

    if (!is.null(xi)) {
        return(ffbs_lambda_from_xi(y1 = y1, xi = xi, gamma = gamma, a0 = a0, b0 = b0))
    }

    filt <- filter_lambda(y1, e1, x1, beta0, beta, phi, kappa, gamma, a0, b0, Xmat = Xmat)
    lam <- backward_sample_lambda(filt$a, filt$b, gamma)

    list(lambda_tilde_0T = lam, a = filt$a, b = filt$b, xi = filt$xi)
}
