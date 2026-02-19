
# utils.R
# Core utilities used throughout the MSSTNB Algorithm 2 sampler.

set_seed <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  invisible(NULL)
}

assert_true <- function(cond, msg = "Assertion failed.") {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
  invisible(NULL)
}

# Stable log-sum-exp
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Clamp to avoid overflow in exp()
clamp <- function(x, lo = -700, hi = 700) pmin(pmax(x, lo), hi)

# Dirichlet draw with vector alpha (>0)
rdirichlet1 <- function(alpha) {
  assert_true(all(alpha > 0), "Dirichlet concentration must be positive.")
  z <- rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

# Dirichlet draw for matrix of alpha: one draw per row
rdirichlet_mat <- function(alpha_mat) {
  n <- nrow(alpha_mat)
  out <- matrix(NA_real_, nrow = n, ncol = ncol(alpha_mat))
  for (i in seq_len(n)) out[i, ] <- rdirichlet1(alpha_mat[i, ])
  out
}

# Build ICAR intrinsic precision H from an undirected edge list (i,j), 1-based
# and optional weights g_ij. If weights missing, set to 1.
build_H_from_edges <- function(edges, n1) {
  # edges: data.frame with columns i, j, optionally g
  assert_true(n1 >= 1, "n1 must be >= 1.")
  assert_true(all(c("i","j") %in% names(edges)), "edges must have columns i and j.")
  if (!("g" %in% names(edges))) edges$g <- 1
  edges$i <- as.integer(edges$i)
  edges$j <- as.integer(edges$j)
  edges$g <- as.numeric(edges$g)

  # Use sparse matrix
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  Matrix <- getNamespace("Matrix")

  ii <- c(edges$i, edges$j)
  jj <- c(edges$j, edges$i)
  vv <- c(edges$g, edges$g)
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = vv, dims = c(n1, n1))
  # ensure zero diagonal
  diag(W) <- 0

  d <- Matrix::rowSums(W)
  H <- Matrix::Diagonal(n = n1, x = d) - W
  Matrix::forceSymmetric(H)
}

# Orthonormal basis B spanning {phi: 1'phi = 0}.
# Uses QR of 1-vector. Returns B (n1 x (n1-1)) and q1 (first column).
build_B_basis <- function(n1) {
    one <- matrix(1, nrow = n1, ncol = 1)
    Q <- qr(one)
    Bfull <- qr.Q(Q, complete = TRUE)
    list(q1 = Bfull[, 1, drop = FALSE],
         B  = Bfull[, -1, drop = FALSE])
}
