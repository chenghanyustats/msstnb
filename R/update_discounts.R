
# update_discounts.R
# Grid-based sampling-importance-resampling for discount factors gamma_j and delta_lj (Section 8.6).
# Uses one-step predictive factors: NegBin for coarsest counts, DirMult for splits.

log_nb_predictive <- function(y, a_prev, b_prev, xi, gamma) {
    # a_prev, b_prev are scalars; xi >= 0
    alpha <- gamma * a_prev
    beta  <- gamma * b_prev
    alpha <- pmax(alpha, 1e-12)
    beta  <- pmax(beta,  1e-12)
    xi    <- pmax(xi,    1e-300)
    y     <- as.numeric(y)

    # log pmf
    lgamma(y + alpha) - lgamma(alpha) - lgamma(y + 1) +
      alpha * (log(beta) - log(beta + xi)) +
      y * (log(xi) - log(beta + xi))
}

update_gamma_grid <- function(y1, xi, a0, b0, a_prior = 20, b_prior = 20,
                              grid = seq(0.05, 0.999, length.out = 200L),
                              time_idx = NULL) {
    TT <- nrow(y1); n1 <- ncol(y1)
    assert_true(all(dim(xi) == c(TT, n1)), "xi must be T x n1.")
    if (is.null(time_idx)) time_idx <- seq_len(TT)

    if (length(a0) == 1) a0 <- rep(a0, n1)
    if (length(b0) == 1) b0 <- rep(b0, n1)

    # shrink grid adaptively, or update
    # 𝛾
    # γ less frequently, or vectorize parts of the computation
    gamma_new <- numeric(n1)
    for (j in seq_len(n1)) {
        logw <- rep(-Inf, length(grid))
        for (m in seq_along(grid)) {
            g <- grid[m]
            # Beta prior
            lp <- dbeta(g, shape1 = a_prior, shape2 = b_prior, log = TRUE)
            a_prev <- a0[j]; b_prev <- b0[j]
            for (t in time_idx) {
              lp <- lp + log_nb_predictive(y1[t, j], a_prev, b_prev, xi[t, j], g)
              # update filter state under this g
              a_prev <- g * a_prev + y1[t, j]
              b_prev <- g * b_prev + xi[t, j]
            }
            logw[m] <- lp
        }
        # resample from discrete weights
        logw <- logw - max(logw)
        w <- exp(logw); w <- w / sum(w)
        gamma_new[j] <- sample(grid, size = 1, prob = w)
    }
    gamma_new
}

log_dirmult_predictive <- function(yvec, alpha) {
    # yvec length K nonnegative integers, alpha length K positive
    yvec <- as.numeric(yvec)
    alpha <- pmax(alpha, 1e-12)
    N <- sum(yvec)
    A <- sum(alpha)
    logcoef <- lgamma(N + 1) - sum(lgamma(yvec + 1))
    logcoef + lgamma(A) - lgamma(N + A) + sum(lgamma(yvec + alpha) - lgamma(alpha))
}

update_delta_grid <- function(tree, y_levels, alpha_dir = 50, q0 = NULL,
                              delta, a_prior = 20, b_prior = 20,
                              grid = seq(0.05, 0.999, length.out = 200L),
                              time_idx = NULL) {
    TT <- nrow(y_levels[[1]])
    if (is.null(time_idx)) time_idx <- seq_len(TT)

    # need prior c0 to run recursion inside each candidate
    c0 <- init_c0_splits(tree, alpha_dir = alpha_dir, q0 = q0)

    L <- tree$L
    delta_new <- delta

    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K <= 1) next

            y_child <- y_levels[[l + 1]][, ch, drop = FALSE]  # T x K
            logw <- rep(-Inf, length(grid))
            for (m in seq_along(grid)) {
                d <- grid[m]
                lp <- dbeta(d, shape1 = a_prior, shape2 = b_prior, log = TRUE)
                c_prev <- c0[[l]][[j]]
                for (t in time_idx) {
                    alpha <- d * c_prev
                    lp <- lp + log_dirmult_predictive(y_child[t, ], alpha)
                    # c_prev <- d * c_prev + y_child[t, ]
                    c_prev <- alpha + y_child[t, ]
                }
                logw[m] <- lp
            }
            logw <- logw - max(logw)
            w <- exp(logw); w <- w / sum(w)
            delta_new[[l]][[j]] <- sample(grid, size = 1, prob = w)
        }
    }
    delta_new
}
