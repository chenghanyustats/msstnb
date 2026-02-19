# update_beta.R
# Elliptical slice sampling (ESS) update for (beta0, beta).
#
# Optimizations:
# - vectorized likelihood evaluation (no loop over regions)
# - optional cached Xmat (TT*n1 x p)
# - optional cached prior Cholesky for (beta0, beta)

prep_beta_prior <- function(p, m0 = 0, s0 = 10, m_beta = NULL,
                            V_beta = NULL) {
    if (is.null(m_beta)) {
        m_beta <- rep(0, p)
    }
    if (is.null(V_beta)) {
        V_beta <- diag(p) * 100
    }

    assert_true(s0 > 0, "s0 must be positive.")
    assert_true(length(m_beta) == p, "m_beta length mismatch.")
    assert_true(all(dim(V_beta) == c(p, p)), "V_beta must be p x p.")

    k <- 1L + p
    m <- c(m0, m_beta)

    Sigma <- matrix(0, nrow = k, ncol = k)
    Sigma[1, 1] <- s0^2
    if (p > 0L) {
        Sigma[-1, -1] <- V_beta
    }

    cholSigma <- chol_spd(Sigma)
    list(k = k, m = m, cholSigma = cholSigma)
}

prep_beta_ess <- function(y1, e1, kappa, lambda_tilde_0T, phi, Xmat = NULL,
                          x1 = NULL) {
    TT <- nrow(y1)
    n1 <- ncol(y1)

    assert_true(all(dim(e1) == c(TT, n1)), "e1 must be T x n1.")
    assert_true(all(dim(kappa) == c(TT, n1)), "kappa must be T x n1.")
    assert_true(all(dim(lambda_tilde_0T) == c(TT + 1, n1)),
                "lambda_tilde_0T must be (T+1) x n1.")
    assert_true(length(phi) == n1, "phi length mismatch.")

    lam_t <- lambda_tilde_0T[-1, , drop = FALSE]

    y_vec <- as.numeric(y1)

    log_cmult_vec <- as.numeric(
        log(pmax(e1, 1e-300)) +
            log(pmax(kappa, 1e-300)) +
            log(pmax(lam_t, 1e-300))
    )

    phi_rep <- rep.int(as.numeric(phi), each = TT)

    if (is.null(Xmat)) {
        assert_true(!is.null(x1), "Provide x1 when Xmat is NULL.")
        p <- dim(x1)[3]
        if (p > 0L) {
            Xmat <- prep_Xmat(x1)
        } else {
            Xmat <- NULL
        }
    } else {
        p <- ncol(Xmat)
    }

    list(
        TT = TT,
        n1 = n1,
        p = p,
        y_vec = y_vec,
        log_cmult_vec = log_cmult_vec,
        phi_rep = phi_rep,
        Xmat = Xmat
    )
}

update_beta_ess_fast <- function(beta0, beta, prep, prior) {
    p <- prep$p
    assert_true(length(beta) == p, "beta length mismatch.")
    assert_true(prior$k == (1L + p), "prior dimension mismatch.")

    z0 <- c(beta0, beta)

    draw_nu <- function() {
        draw_nu_cholSigma(prior$cholSigma)
    }

    loglik_fn <- function(z) {
        b0 <- z[1]
        b  <- z[-1]

        if (p > 0L) {
            xb_vec <- as.numeric(prep$Xmat %*% b)
        } else {
            xb_vec <- 0
        }

        eta_vec <- b0 + xb_vec + prep$phi_rep
        ce <- clamp_eta(eta_vec)

        log_mu_vec <- prep$log_cmult_vec + ce
        mu_vec <- exp(clamp(log_mu_vec))

        val <- sum(prep$y_vec * ce - mu_vec)
        if (!is.finite(val)) {
            return(-Inf)
        }
        val
    }

    out <- ess(z0, loglik_fn = loglik_fn, draw_nu_fn = draw_nu, m = prior$m)
    z_new <- out$z

    list(beta0 = z_new[1], beta = z_new[-1], ess_steps = out$n_steps)
}

# Backward-compatible wrapper.
update_beta_ess <- function(beta0, beta, phi, lambda_tilde_0T, kappa, y1, e1, x1,
                            m0 = 0, s0 = 10, m_beta = NULL, V_beta = NULL,
                            Xmat = NULL, prior_prep = NULL) {
    p <- dim(x1)[3]

    if (is.null(prior_prep)) {
        prior_prep <- prep_beta_prior(
            p = p,
            m0 = m0,
            s0 = s0,
            m_beta = m_beta,
            V_beta = V_beta
        )
    }

    prep <- prep_beta_ess(
        y1 = y1,
        e1 = e1,
        kappa = kappa,
        lambda_tilde_0T = lambda_tilde_0T,
        phi = phi,
        Xmat = Xmat,
        x1 = x1
    )

    update_beta_ess_fast(beta0 = beta0, beta = beta, prep = prep,
                         prior = prior_prep)
}
