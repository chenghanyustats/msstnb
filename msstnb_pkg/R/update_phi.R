# update_phi.R
# ICAR-aware ESS update for phi via reduced coordinates u = B' phi.
#
# Optimizations:
# - Cholesky of (B' H B) is cached once (chol_BHB)
# - no Matrix factorization inside the MCMC loop
# - optional rep_idx for fast phi replication to TT*n1 vector

update_phi_ess <- function(phi, tau_phi, H, beta0, beta, lambda_tilde_0T, kappa,
                           y1, e1, x1, B = NULL, Xmat = NULL,
                           chol_BHB = NULL, rep_idx = NULL) {
    TT <- nrow(y1)
    n1 <- length(phi)
    p  <- dim(x1)[3]

    assert_true(tau_phi > 0, "tau_phi must be positive.")

    if (is.null(B)) {
        B <- build_B_basis(n1)$B
    }

    assert_true(nrow(B) == n1 && ncol(B) == (n1 - 1L), "B dimension mismatch.")

    # Current reduced coordinate (B has orthonormal columns)
    u0 <- as.numeric(crossprod(B, phi))

    if (is.null(rep_idx)) {
        rep_idx <- rep(seq_len(n1), each = TT)
    }

    # Cache chol(B' H B) if not provided.
    # This is dense (n1-1)x(n1-1) and usually small at the coarsest level.
    if (is.null(chol_BHB)) {
        if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Package 'Matrix' is required when chol_BHB is not provided.",
                 call. = FALSE)
        }
        HB <- as.matrix(H %*% B)
        BHB <- crossprod(B, HB)
        chol_BHB <- chol_spd(BHB)
    }



    # Prior: u ~ N(0, (tau_phi * BHB)^{-1})
    # If chol_BHB is upper R where BHB = R'R,
    # then sample nu by solving (sqrt(tau)*R) nu = z.
    draw_nu <- function() {
        # z <- rnorm(n1 - 1L)
        backsolve(chol_BHB, rnorm(n1 - 1L) / sqrt(tau_phi))
    }

    # Likelihood pieces that do not depend on phi proposal
    lam_t <- lambda_tilde_0T[-1, , drop = FALSE]
    y_vec <- as.numeric(y1)

    log_cmult_vec <- as.numeric(
        log(pmax(e1, 1e-300)) + log(pmax(kappa, 1e-300)) +
          log(pmax(lam_t, 1e-300))
    )

    if (p > 0L) {
        if (is.null(Xmat)) {
            Xmat <- prep_Xmat(x1)
        }
        xb_vec <- as.numeric(Xmat %*% beta)
    } else {
        xb_vec <- 0
    }

    eta_base_vec <- beta0 + xb_vec

    loglik_fn <- function(u) {
        phi_u <- as.numeric(B %*% u)
        # Numerical safety: enforce sum-to-zero.
        phi_u <- phi_u - mean(phi_u)

        phi_rep <- phi_u[rep_idx]
        eta_vec <- eta_base_vec + phi_rep
        ce <- clamp_eta(eta_vec)

        log_mu_vec <- log_cmult_vec + ce
        mu_vec <- exp(clamp(log_mu_vec))

        val <- sum(y_vec * ce - mu_vec)
        if (!is.finite(val)) {
            return(-Inf)
        }
        val
    }

    out <- ess(u0, loglik_fn = loglik_fn, draw_nu_fn = draw_nu,
               m = rep(0, n1 - 1L))

    u_new <- out$z
    phi_new <- as.numeric(B %*% u_new)
    phi_new <- phi_new - mean(phi_new)

    list(phi = phi_new, u = u_new, ess_steps = out$n_steps, B = B)
}
