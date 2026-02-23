
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
    TT <- nrow(kappa)
    n1 <- ncol(kappa)
    assert_true(length(r1) == n1, "r1 length mismatch.")
    assert_true(is.finite(prop_sd) && prop_sd > 0,
                "prop_sd must be positive.")

    out <- r1
    acc <- logical(n1)

    for (j in seq_len(n1)) {
        logr <- log(pmax(r1[j], 1e-8))
        lp <- logpost_logr(logr, kappa[, j], a_r, b_r)

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


# Adaptive tuning for the random-walk MH proposal sd on log(r).
#
# We adapt on the log scale using a Robbins-Monro style update:
#   log(sd_new) = log(sd_old) + gamma_t * (acc_bar - target)
# where gamma_t decays with iteration: gamma_t = gamma0 / (iter + t0)^power.
#
# Inputs
#   - prop_sd: current proposal sd
#   - acc_bar: recent acceptance rate summary (e.g., mean over last M iters)
#   - target: target acceptance rate (e.g., 0.30)
#   - iter: current iteration index (1-based)
#   - gamma0, t0, power: step-size parameters
#   - min_sd, max_sd: clamp bounds
#
# Output
#   - new proposal sd
adapt_prop_sd_r <- function(prop_sd, acc_bar, target = 0.30, iter,
                            gamma0 = 0.05, t0 = 10, power = 0.5,
                            min_sd = 0.01, max_sd = 2.0) {
    assert_true(is.finite(prop_sd) && prop_sd > 0,
                "prop_sd must be positive.")
    assert_true(is.finite(acc_bar),
                "acc_bar must be finite.")
    assert_true(is.finite(target) && target > 0 && target < 1,
                "target must be in (0, 1).")
    assert_true(is.finite(iter) && iter >= 1,
                "iter must be >= 1.")
    assert_true(is.finite(gamma0) && gamma0 > 0,
                "gamma0 must be positive.")
    assert_true(is.finite(t0) && t0 >= 0,
                "t0 must be >= 0.")
    assert_true(is.finite(power) && power > 0,
                "power must be positive.")
    assert_true(is.finite(min_sd) && min_sd > 0,
                "min_sd must be positive.")
    assert_true(is.finite(max_sd) && max_sd > min_sd,
                "max_sd must be > min_sd.")

    gamma_t <- gamma0 / ((as.numeric(iter) + as.numeric(t0)) ^ power)
    log_sd_new <- log(prop_sd) + gamma_t * (acc_bar - target)
    sd_new <- exp(log_sd_new)
    sd_new <- min(max(sd_new, min_sd), max_sd)
    sd_new
}
