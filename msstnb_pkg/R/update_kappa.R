# update_kappa.R
# Update NB random effects kappa_{t,j} (coarsest level) via (27).
#
# Key optimizations:
# - optional cached Xmat (TT*n1 x p) avoids reshaping x1 each iteration
# - avoids constructing a large phi matrix via sweep()

update_kappa <- function(y1, e1, x1, beta0, beta, phi, lambda_tilde_0T, r1,
                         Xmat = NULL) {
    TT <- nrow(y1)
    n1 <- ncol(y1)
    p  <- dim(x1)[3]

    assert_true(length(r1) == n1, "r1 must have length n1.")
    assert_true(all(dim(e1) == c(TT, n1)), "e1 dim mismatch.")
    assert_true(length(phi) == n1, "phi length mismatch.")
    assert_true(all(dim(lambda_tilde_0T) == c(TT + 1, n1)),
                "lambda_tilde_0T must be (TT+1) x n1.")

    # t = 1..T
    lam_t <- lambda_tilde_0T[-1, , drop = FALSE]

    if (p > 0L) {
        if (is.null(Xmat)) {
            Xmat <- prep_Xmat(x1)
        }
        xb_vec <- as.numeric(Xmat %*% beta)
        xb_mat <- matrix(xb_vec, nrow = TT, ncol = n1)
    } else {
        xb_mat <- matrix(0, nrow = TT, ncol = n1)
    }

    eta <- beta0 + xb_mat
    eta <- sweep(eta, 2, phi, FUN = "+")

    # Effective offset excluding kappa
    zeta <- e1 * exp(clamp_eta(eta)) * lam_t

    # Full cond
    # kappa[t,j] ~ Gamma(shape = y1[t,j] + r1[j], rate = r1[j] + zeta[t,j])
    r1_mat <- matrix(r1, nrow = TT, ncol = n1, byrow = TRUE)
    shape <- y1 + r1_mat
    rate  <- zeta + r1_mat

    # Safeguards
    shape <- pmax(shape, 1e-8)
    rate  <- pmax(rate,  1e-12)

    kappa_vec <- rgamma(TT * n1, shape = as.numeric(shape),
                        rate = as.numeric(rate))
    kappa <- matrix(kappa_vec, nrow = TT, ncol = n1)

    list(kappa = kappa, zeta = zeta)
}
