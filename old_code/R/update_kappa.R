
# update_kappa.R
# Update NB random effects kappa_{t,1j} via (27).

update_kappa <- function(y1, e1, x1, beta0, beta, phi, lambda_tilde_0T, r1) {
    TT <- nrow(y1)
    n1 <- ncol(y1)
    p  <- dim(x1)[3]
    TTn1 <- TT * n1

    assert_true(length(r1) == n1, "r1 must have length n1.")
    assert_true(nrow(e1) == TT && ncol(e1) == n1, "e1 dim mismatch.")
    assert_true(length(phi) == n1, "phi must have length n1.")
    assert_true(nrow(lambda_tilde_0T) == TT + 1 && ncol(lambda_tilde_0T) == n1,
                "lambda_tilde_0T must be (TT+1) by n1.")

    # t = 1..T
    lam_t <- lambda_tilde_0T[-1, , drop = FALSE]

    # Build eta = beta0 + X beta + phi
    # x1 is TT x n1 x p. Reshape to (TT*n1) x p, multiply once, reshape back to TT x n1.
    if (p > 0) {
        Xmat <- matrix(x1, nrow = TTn1, ncol = p)     # column k is vec(x1[,,k])
        xb   <- as.numeric(Xmat %*% beta)             # length TT*n1
        xb   <- matrix(xb, nrow = TT, ncol = n1)      # fills column major: time within region
    } else {
        xb <- matrix(0, nrow = TT, ncol = n1)
    }

    eta <- beta0 + xb + matrix(phi, nrow = TT, ncol = n1, byrow = TRUE)

    # Effective offset excluding kappa
    zeta <- e1 * exp(clamp(eta)) * lam_t

    # Full cond
    # kappa[t,j] ~ Gamma(shape = y1[t,j] + r1[j], rate = r1[j] + zeta[t,j])
    r1_mat <- matrix(r1, nrow = TT, ncol = n1, byrow = TRUE)
    shape <- y1 + r1_mat
    rate  <- zeta + r1_mat

    # Safeguards
    shape <- pmax(shape, 1e-8)
    rate  <- pmax(rate,  1e-12)

    kappa_vec <- rgamma(TTn1, shape = as.numeric(shape),
                        rate = as.numeric(rate))
    kappa <- matrix(kappa_vec, nrow = TT, ncol = n1)

    list(kappa = kappa, zeta = zeta)
}


# update_kappa <- function(y1, e1, x1, beta0, beta, phi, lambda_tilde_0T, r1) {
#     TT <- nrow(y1); n1 <- ncol(y1); p <- dim(x1)[3]
#     assert_true(length(r1) == n1, "r1 must have length n1.")
#     # compute eta
#     eta <- matrix(beta0, nrow = TT, ncol = n1)
#     if (p > 0) {
#         for (j in seq_len(n1)) {
#             eta[, j] <- eta[, j] + as.numeric(x1[, j, , drop = TRUE] %*% beta)
#         }
#     }
#     eta <- eta + matrix(rep(phi, each = TT), nrow = TT, ncol = n1)
#
#     lam_t <- lambda_tilde_0T[-1, , drop = FALSE] # t=1..T
#     zeta <- e1 * exp(clamp(eta)) * lam_t  # effective offset excluding kappa
#
#     # full conditional: kappa | ... ~ Ga(y+r, r+zeta)
#     kappa <- matrix(NA_real_, nrow = TT, ncol = n1)
#
#     for (j in seq_len(n1)) {
#       shape <- y1[, j] + r1[j]
#       rate  <- r1[j] + zeta[, j]
#       # safeguard
#       shape <- pmax(shape, 1e-8)
#       rate  <- pmax(rate,  1e-12)
#       kappa[, j] <- rgamma(TT, shape = shape, rate = rate)
#     }
#
#     list(kappa = kappa, zeta = zeta)
# }
