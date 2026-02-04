\
# utils_io.R
# Small helpers for simulation code

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

parse_num_vec <- function(x) {
  # Parse "0.2, -0.1" into numeric vector
  if (is.numeric(x)) return(as.numeric(x))
  if (is.na(x) || nchar(trimws(x)) == 0) return(numeric(0))
  parts <- strsplit(gsub("\"", "", x), ",")[[1]]
  as.numeric(trimws(parts))
}

parse_int_vec <- function(x) {
  as.integer(parse_num_vec(x))
}

rdirichlet1 <- function(alpha) {
  # single draw from Dirichlet(alpha)
  assert_true(all(alpha > 0), "Dirichlet alpha must be positive.")
  z <- rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

rdirichlet <- function(n, alpha) {
  # n draws, return n by K
  K <- length(alpha)
  out <- matrix(NA_real_, nrow = n, ncol = K)
  for (i in seq_len(n)) out[i, ] <- rdirichlet1(alpha)
  out
}

rmultinom1 <- function(size, prob) {
  # single multinomial draw, returns integer vector
  as.integer(rmultinom(1, size = size, prob = prob)[, 1])
}

set_seed <- function(seed) {
  if (!is.null(seed) && !is.na(seed)) set.seed(as.integer(seed))
}

