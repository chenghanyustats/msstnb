# update_discounts.R
# Grid-based sampling-importance-resampling for discount factors gamma_j and delta_{l,j}.
#
# Major optimization: vectorized recursion over grid points.
# - For each region/node, we maintain (a_prev, b_prev) or c_prev for all grid values.
# - Then we loop only over time (often much smaller than grid length).

log_nb_predictive_vec <- function(y, alpha, beta, xi) {
    # Vectorized version of log NB predictive.
    # y is scalar, alpha/beta are vectors (same length), xi is scalar.
    alpha <- pmax(alpha, 1e-12)
    beta  <- pmax(beta,  1e-12)
    xi    <- pmax(xi,    1e-300)

    y <- as.numeric(y)

    lgamma(y + alpha) - lgamma(alpha) - lgamma(y + 1) +
        alpha * (log(beta) - log(beta + xi)) +
        y * (log(xi) - log(beta + xi))
}

update_gamma_grid <- function(y1, xi, a0, b0,
                              a_prior = 20, b_prior = 20,
                              grid = seq(0.05, 0.999, length.out = 200L),
                              time_idx = NULL) {
    TT <- nrow(y1)
    n1 <- ncol(y1)

    assert_true(all(dim(xi) == c(TT, n1)), "xi must be T x n1.")

    if (is.null(time_idx)) {
        time_idx <- seq_len(TT)
    }

    grid <- as.numeric(grid)
    M <- length(grid)

    if (length(a0) == 1L) {
        a0 <- rep(a0, n1)
    }
    if (length(b0) == 1L) {
        b0 <- rep(b0, n1)
    }

    assert_true(length(a0) == n1, "a0 must be scalar or length n1.")
    assert_true(length(b0) == n1, "b0 must be scalar or length n1.")

    log_prior <- dbeta(grid, shape1 = a_prior, shape2 = b_prior, log = TRUE)

    gamma_new <- numeric(n1)

    # shrink grid adaptively, or update
    # 𝛾
    # γ less frequently, or vectorize parts of the computation
    for (j in seq_len(n1)) {
        logw <- log_prior

        a_prev <- rep(a0[j], M)
        b_prev <- rep(b0[j], M)

        for (t in time_idx) {
            y <- y1[t, j]
            xit <- xi[t, j]

            alpha <- grid * a_prev
            beta  <- grid * b_prev

            logw <- logw + log_nb_predictive_vec(y = y, alpha = alpha,
                                                 beta = beta, xi = xit)

            a_prev <- alpha + y
            b_prev <- beta + xit
        }

        logw <- logw - max(logw)
        w <- exp(logw)
        w <- w / sum(w)

        gamma_new[j] <- sample(grid, size = 1, prob = w)
    }

    gamma_new
}

log_dirmult_predictive_vec <- function(yvec, alpha_mat) {
    # yvec: length K
    # alpha_mat: M x K
    yvec <- as.numeric(yvec)
    alpha_mat <- pmax(alpha_mat, 1e-12)

    N <- sum(yvec)
    logcoef <- lgamma(N + 1) - sum(lgamma(yvec + 1))

    A <- rowSums(alpha_mat)
    term1 <- lgamma(A) - lgamma(N + A)

    y_mat <- matrix(yvec, nrow = nrow(alpha_mat), ncol = length(yvec),
                    byrow = TRUE)
    term2 <- rowSums(lgamma(y_mat + alpha_mat) - lgamma(alpha_mat))

    logcoef + term1 + term2
}

update_delta_grid <- function(tree, y_levels,
                              alpha_dir = 50, q0 = NULL,
                              delta,
                              a_prior = 20, b_prior = 20,
                              grid = seq(0.05, 0.999, length.out = 200L),
                              time_idx = NULL) {
    TT <- nrow(y_levels[[1]])
    if (is.null(time_idx)) {
        time_idx <- seq_len(TT)
    }

    grid <- as.numeric(grid)
    M <- length(grid)

    log_prior <- dbeta(grid, shape1 = a_prior, shape2 = b_prior, log = TRUE)

    c0 <- init_c0_splits(tree, alpha_dir = alpha_dir, q0 = q0)

    L <- tree$L
    delta_new <- delta

    for (l in seq_len(L - 1L)) {
        n_parent <- tree$n[l]

        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K <= 1L) next

            y_child <- y_levels[[l + 1L]][, ch, drop = FALSE]  # T x K

            # State for all grid values: M x K
            c_prev <- matrix(rep(c0[[l]][[j]], each = M), nrow = M, ncol = K)
            logw <- log_prior

            for (t in time_idx) {
                yvec <- y_child[t, ]

                # alpha_mat[m, ] = grid[m] * c_prev[m, ]
                alpha_mat <- c_prev * grid

                logw <- logw + log_dirmult_predictive_vec(yvec = yvec,
                                                          alpha_mat = alpha_mat)

                # update state: c_prev <- alpha + y
                c_prev <- alpha_mat + matrix(yvec, nrow = M, ncol = K,
                                             byrow = TRUE)
            }

            logw <- logw - max(logw)
            w <- exp(logw)
            w <- w / sum(w)

            delta_new[[l]][[j]] <- sample(grid, size = 1, prob = w)
        }
    }

    delta_new
}
