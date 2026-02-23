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
# build_H_from_edges <- function(edges, n1) {
#     assert_true(n1 >= 1, "n1 must be >= 1.")
#     assert_true(all(c("i", "j") %in% names(edges)), "edges must have columns i and j.")
#
#     if (!requireNamespace("Matrix", quietly = TRUE)) {
#         stop("Package 'Matrix' is required.", call. = FALSE)
#     }
#
#     if (!("g" %in% names(edges))) {
#         edges$g <- 1
#     }
#
#     edges$i <- as.integer(edges$i)
#     edges$j <- as.integer(edges$j)
#     edges$g <- as.numeric(edges$g)
#
#     ii <- c(edges$i, edges$j)
#     jj <- c(edges$j, edges$i)
#     vv <- c(edges$g, edges$g)
#
#     W <- Matrix::sparseMatrix(i = ii, j = jj, x = vv, dims = c(n1, n1))
#     diag(W) <- 0
#
#     d <- Matrix::rowSums(W)
#     H <- Matrix::Diagonal(n = n1, x = d) - W
#     Matrix::forceSymmetric(H)
# }


# Build ICAR intrinsic precision H from an undirected edge list.
# Accepts edges as:
#   1) data.frame with columns i, j (optional g)
#   2) data.frame with columns from, to (optional g) as in igraph
#   3) matrix with at least 2 columns (treated as i, j)
# Indices are assumed 1 based. If 0 based indices are detected, they are shifted to 1 based.
build_H_from_edges <- function(edges, n1) {
  assert_true(n1 >= 1L, "n1 must be >= 1.")
  assert_true(!is.null(edges), "edges is NULL.")

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }

  to_int <- function(z) {
    if (is.factor(z)) z <- as.character(z)
    z <- suppressWarnings(as.integer(z))
    z
  }

  # Normalize edges to data.frame
  if (is.matrix(edges)) {
    assert_true(ncol(edges) >= 2L, "edges matrix must have at least 2 columns.")
    edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  } else if (!is.data.frame(edges)) {
    edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  }

  nm <- names(edges)
  if (is.null(nm)) nm <- character(0)

  # Map common column names to i and j
  if (all(c("from", "to") %in% nm)) {
    names(edges)[match("from", names(edges))] <- "i"
    names(edges)[match("to", names(edges))] <- "j"
  } else if (all(c("v1", "v2") %in% nm)) {
    names(edges)[match("v1", names(edges))] <- "i"
    names(edges)[match("v2", names(edges))] <- "j"
  } else if (!all(c("i", "j") %in% nm)) {
    # No usable names, treat first two columns as i and j
    assert_true(ncol(edges) >= 2L, "edges must have at least 2 columns.")
    edges$i <- edges[[1]]
    edges$j <- edges[[2]]
  }

  if (!("g" %in% names(edges))) {
    edges$g <- 1
  }

  edges <- edges[, c("i", "j", "g"), drop = FALSE]

  edges$i <- to_int(edges$i)
  edges$j <- to_int(edges$j)
  edges$g <- as.numeric(edges$g)

  assert_true(!anyNA(edges$i) && !anyNA(edges$j), "edges i or j contain non numeric values.")
  assert_true(all(is.finite(edges$g)) && all(edges$g > 0), "edge weights g must be positive.")

  # Detect 0 based indexing and shift to 1 based
  min_idx <- min(c(edges$i, edges$j))
  if (is.finite(min_idx) && min_idx == 0L) {
    edges$i <- edges$i + 1L
    edges$j <- edges$j + 1L
  }

  assert_true(all(edges$i >= 1L & edges$i <= n1), "edges$i out of range.")
  assert_true(all(edges$j >= 1L & edges$j <= n1), "edges$j out of range.")

  # Remove self loops
  edges <- edges[edges$i != edges$j, , drop = FALSE]

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
