# utils.R
# Core utilities used throughout the MSSTNB Algorithm 2 sampler.
#
# Design goals:
# - predictable shapes
# - numerical safeguards
# - minimal dependencies

set_seed <- function(seed = NULL) {
    if (!is.null(seed)) {
        set.seed(as.integer(seed))
    }
    invisible(NULL)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

assert_true <- function(cond, msg = "Assertion failed.") {
    if (!isTRUE(cond)) {
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

# Stable log-sum-exp
logsumexp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) {
        return(m)
    }
    m + log(sum(exp(x - m)))
}

# Clamp helper (generic)
clamp <- function(x, lo = -700, hi = 700) {
    pmin(pmax(x, lo), hi)
}

# Clamp specifically for linear predictors used inside exp(eta)
# (tighter bounds reduce Inf/NaN risk in intensities).
clamp_eta <- function(x, lo = -30, hi = 30) {
    pmin(pmax(x, lo), hi)
}

# Dirichlet draw with vector alpha (>0)
rdirichlet1 <- function(alpha) {
    alpha <- as.numeric(alpha)
    assert_true(all(alpha > 0), "Dirichlet concentration must be positive.")
    z <- rgamma(length(alpha), shape = alpha, rate = 1)
    z / sum(z)
}

# Dirichlet draws for a matrix of alpha: one draw per row.
# Vectorized: uses a single rgamma call and row normalization.
rdirichlet_mat <- function(alpha_mat) {
    alpha_mat <- as.matrix(alpha_mat)
    assert_true(all(alpha_mat > 0), "Dirichlet concentration must be positive.")

    n <- nrow(alpha_mat)
    k <- ncol(alpha_mat)

    z <- matrix(
        rgamma(n * k, shape = as.numeric(alpha_mat), rate = 1),
        nrow = n,
        ncol = k
    )

    rs <- rowSums(z)
    rs <- pmax(rs, 1e-300)
    z / rs
}

# Build ICAR intrinsic precision H from an undirected edge list (i,j), 1-based.
# Optional weights g_ij (defaults to 1).
build_H_from_edges <- function(edges, n1) {
    assert_true(n1 >= 1, "n1 must be >= 1.")
    assert_true(all(c("i", "j") %in% names(edges)), "edges must have columns i and j.")

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required.", call. = FALSE)
    }

    if (!("g" %in% names(edges))) {
        edges$g <- 1
    }

    edges$i <- as.integer(edges$i)
    edges$j <- as.integer(edges$j)
    edges$g <- as.numeric(edges$g)

    ii <- c(edges$i, edges$j)
    jj <- c(edges$j, edges$i)
    vv <- c(edges$g, edges$g)

    W <- Matrix::sparseMatrix(i = ii, j = jj, x = vv, dims = c(n1, n1))
    diag(W) <- 0

    d <- Matrix::rowSums(W)
    H <- Matrix::Diagonal(n = n1, x = d) - W
    Matrix::forceSymmetric(H)
}

# Orthonormal basis B spanning {phi: 1'phi = 0}.
# Uses QR of the 1-vector. Returns B (n1 x (n1-1)) and q1 (first column).
build_B_basis <- function(n1) {
    one <- matrix(1, nrow = n1, ncol = 1)
    Q <- qr(one)
    Bfull <- qr.Q(Q, complete = TRUE)

    list(
        q1 = Bfull[, 1, drop = FALSE],
        B  = Bfull[, -1, drop = FALSE]
    )
}

# Robust Cholesky for symmetric positive definite matrices.
# Adds diagonal jitter if needed.
chol_spd <- function(A, jitter = 1e-10, max_tries = 6L) {
    A <- as.matrix(A)
    d <- nrow(A)
    assert_true(ncol(A) == d, "A must be square.")

    for (k in 0:max_tries) {
        add <- if (k == 0) 0 else jitter * 10^(k - 1)
        out <- try(chol(A + diag(add, d)), silent = TRUE)
        if (!inherits(out, "try-error")) {
            return(out)
        }
    }

    stop("chol_spd failed: matrix not SPD even after jitter.",
         call. = FALSE)
}
