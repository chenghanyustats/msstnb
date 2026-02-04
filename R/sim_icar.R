\
# sim_icar.R
# Simulate constrained ICAR spatial effects with sum to zero constraint

suppressWarnings({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package Matrix is required.", call. = FALSE)
  }
})

source_local <- function(path) {
  source(path, local = TRUE)
}

sim_phi_icar <- function(H, tau_phi, seed = NULL) {
  # H is intrinsic precision matrix (n1 by n1), symmetric
  # Simulate phi with constraint 1' phi = 0
  # Using reduced coordinates phi = B u, u ~ N(0, (tau_phi B' H B)^{-1})
  set_seed(seed)

  n1 <- nrow(H)
  assert_true(ncol(H) == n1, "H must be square.")
  assert_true(tau_phi > 0, "tau_phi must be positive.")

  one <- matrix(1, nrow = n1, ncol = 1)
  Q <- qr(one)
  Bfull <- qr.Q(Q, complete = TRUE)
  B <- Bfull[, -1, drop = FALSE]  # n1 by (n1 - 1)

  # Reduced precision
  # Using sparse multiplication
  BH <- Matrix::t(Matrix::Matrix(B)) %*% H
  Qsub <- tau_phi * (BH %*% B)
  Qsub <- Matrix::forceSymmetric(Qsub)

  # Cholesky of precision
  cholQ <- Matrix::Cholesky(Qsub, LDL = FALSE, perm = TRUE)
  z <- rnorm(n1 - 1)
  u <- Matrix::solve(cholQ, z)
  phi <- as.numeric(B %*% u)

  # Numerically enforce centering
  phi <- phi - mean(phi)
  phi
}

