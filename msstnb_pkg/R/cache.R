# cache.R
# One-time caches to avoid recomputation inside the MCMC loop.

build_msstnb_cache <- function(data, hyper, mcmc) {
    y1 <- data$y_levels[[1]]
    x1 <- data$x
    tree <- data$tree

    TT <- nrow(y1)
    n1 <- ncol(y1)
    p  <- dim(x1)[3]

    # Design matrix: (TT*n1) x p
    Xmat <- if (p > 0L) prep_Xmat(x1) else NULL

    # Replication index: maps (TT*n1) vector positions to region index
    rep_idx <- rep(seq_len(n1), each = TT)

    # Delta mapping
    delta_map <- make_delta_map(tree)

    # Beta prior prep for ESS
    beta_prior_prep <- prep_beta_prior(
        p = p,
        m0 = hyper$beta$m0,
        s0 = hyper$beta$s0,
        m_beta = hyper$beta$m_beta,
        V_beta = hyper$beta$V_beta
    )

    # H and B for ICAR
    if (!is.null(data$graph$H)) {
        H <- data$graph$H
    } else {
        H <- build_H_from_edges(data$graph$edges, n1 = n1)
    }
    B <- build_B_basis(n1)$B

    # Cache chol(B' H B)
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.", call. = FALSE)
    }
    HB <- as.matrix(H %*% B)
    BHB <- crossprod(B, HB)
    chol_BHB <- chol_spd(BHB)

    # Discount grids and time indices
    gamma_grid <- as.numeric(mcmc$gamma_grid %||%
                               seq(0.05, 0.999, length.out = 200L))
    delta_grid <- as.numeric(mcmc$delta_grid %||%
                               seq(0.05, 0.999, length.out = 200L))
    discount_time_idx <- mcmc$discount_time_idx %||% NULL

    list(
        TT = TT,
        n1 = n1,
        p = p,
        Xmat = Xmat,
        rep_idx = rep_idx,
        tree = tree,
        delta_map = delta_map,
        beta_prior_prep = beta_prior_prep,
        H = H,
        B = B,
        BHB = BHB,
        chol_BHB = chol_BHB,
        gamma_grid = gamma_grid,
        delta_grid = delta_grid,
        discount_time_idx = discount_time_idx
    )
}
