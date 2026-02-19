# update_phi.R
# ICAR-aware ESS update for phi (via reduced coordinate u).
#
# Speed improvements:
# - Optional cached Xmat (TT*n1 x p) to avoid reshaping x1 repeatedly
# - Precompute beta0 + X beta term once (beta fixed within this block)

update_phi_ess <- function(phi, tau_phi, H, beta0, beta, lambda_tilde_0T, kappa,
                           y1, e1, x1, B = NULL, Xmat = NULL) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.", call. = FALSE)
    }
    Matrix <- getNamespace("Matrix")

    TT <- nrow(y1)
    n1 <- length(phi)
    p  <- dim(x1)[3]

    assert_true(Matrix::nrow(H) == n1 && Matrix::ncol(H) == n1,
                "H dimension mismatch.")
    assert_true(tau_phi > 0, "tau_phi must be positive.")

    # Build basis if needed
    if (is.null(B)) {
        bb <- build_B_basis(n1)
        B <- bb$B
    } else {
        assert_true(nrow(B) == n1 && ncol(B) == (n1 - 1), "B dimension mismatch.")
    }

    # Current reduced coordinate u (since B has orthonormal columns): u = B' phi
    u0 <- as.numeric(t(B) %*% phi)

    # Prior precision on u: Qsub = tau_phi * (B' H B)
    BH <- Matrix::t(Matrix::Matrix(B)) %*% H
    Qsub <- tau_phi * (BH %*% Matrix::Matrix(B))
    Qsub <- Matrix::forceSymmetric(Qsub)

    cholQ <- Matrix::Cholesky(Qsub, LDL = FALSE, perm = TRUE)
    draw_nu <- function() {
        z <- rnorm(n1 - 1)
        v <- as.numeric(Matrix::solve(cholQ, z, system = "Lt"))  # v = L^{-T} z
        as.numeric(Matrix::solve(cholQ, v, system = "Pt"))       # nu = P^T v
    }


    # Precompute components that do not depend on phi proposal
    ## Maybe no need to recompute (check it)
    lam_t <- lambda_tilde_0T[-1, , drop = FALSE]
    y_vec <- as.numeric(y1)
    log_cmult_vec <- as.numeric(
        log(pmax(e1, 1e-300)) + log(pmax(kappa, 1e-300)) +
          log(pmax(lam_t, 1e-300))
    )


    ## Maybe no need to recompute (check it)
    if (p > 0) {
        if (is.null(Xmat)) {
            Xmat <- matrix(x1, nrow = TT * n1, ncol = p)
        } else {
            assert_true(nrow(Xmat) == TT * n1 && ncol(Xmat) == p,
                        "Xmat dim mismatch.")
        }
            xb_vec <- as.numeric(Xmat %*% beta)
    } else {
        xb_vec <- 0
    }

    eta_base_vec <- beta0 + xb_vec

    ## out outside?
    loglik_fn <- function(u) {
        phi_u <- as.numeric(B %*% u)
        phi_u <- phi_u - mean(phi_u) ## centering (why)
        phi_rep <- rep.int(phi_u, each = TT)

        eta_vec <- eta_base_vec + phi_rep
        ce <- clamp(eta_vec)

        log_mu_vec <- log_cmult_vec + ce
        mu_vec <- exp(clamp(log_mu_vec))

        ## change to y_vec * phi_rep?
        val <- sum(y_vec * ce - mu_vec)
        if (!is.finite(val)) return(-Inf)
        val
    }

    out <- ess(u0, loglik_fn = loglik_fn, draw_nu_fn = draw_nu,
               m = rep(0, n1 - 1))
    u_new <- out$z
    phi_new <- as.numeric(B %*% u_new)
    phi_new <- phi_new - mean(phi_new)

    list(phi = phi_new, u = u_new, ess_steps = out$n_steps, B = B)
}
