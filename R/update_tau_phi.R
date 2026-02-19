
# update_tau_phi.R
# Update ICAR precision tau_phi via (31).

update_tau_phi <- function(phi, H, a_phi, b_phi) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required.", call. = FALSE)
    }
    Matrix <- getNamespace("Matrix")

    n1 <- length(phi)
    assert_true(Matrix::nrow(H) == n1 && Matrix::ncol(H) == n1, "H dimension mismatch.")
    assert_true(a_phi > 0 && b_phi > 0, "a_phi and b_phi must be positive.")

    quad <- as.numeric(Matrix::t(phi) %*% (H %*% phi))
    shape <- a_phi + (n1 - 1)/2
    rate  <- b_phi + 0.5 * quad
    rgamma(1, shape = shape, rate = rate)
}
