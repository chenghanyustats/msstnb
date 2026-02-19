# linear_predictor.R
# Fast, shape-stable computation of linear predictors and intensities.

prep_Xmat <- function(x1) {
    TT <- dim(x1)[1]
    n1 <- dim(x1)[2]
    p  <- dim(x1)[3]
    if (p == 0L) {
        return(NULL)
    }
    matrix(x1, nrow = TT * n1, ncol = p)
}

compute_eta <- function(beta0, beta, phi, TT, n1, Xmat = NULL) {
    # Returns TT x n1 matrix.
    p <- length(beta)

    if (p > 0L) {
        assert_true(!is.null(Xmat), "compute_eta: Xmat is required when p > 0.")
        xb_vec <- as.numeric(Xmat %*% beta)
        xb_mat <- matrix(xb_vec, nrow = TT, ncol = n1)
    } else {
        xb_mat <- matrix(0, nrow = TT, ncol = n1)
    }

    eta <- beta0 + xb_mat
    eta <- sweep(eta, 2, phi, FUN = "+")
    eta
}

compute_xi <- function(e1, kappa, eta) {
    # xi_tj = e_tj * kappa_tj * exp(eta_tj)
    assert_true(all(dim(e1) == dim(kappa)),
                "compute_xi: e1 and kappa dim mismatch.")
    assert_true(all(dim(e1) == dim(eta)),
                "compute_xi: eta dim mismatch.")

    xi <- e1 * kappa * exp(clamp_eta(eta))
    xi <- pmax(xi, 1e-300)
    xi
}
