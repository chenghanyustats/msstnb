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

# rmultinom1 <- function(size, prob, max_chunk = .Machine$integer.max) {
#   prob <- as.numeric(prob)
#   K <- length(prob)
#   assert_true(K >= 1, "prob must have positive length.")
#   assert_true(all(is.finite(prob)), "prob must be finite.")
#   prob[prob < 0] <- 0
#   ps <- sum(prob)
#   assert_true(is.finite(ps) && ps > 0, "prob must sum to a positive finite value.")
#   prob <- prob / ps
#
#   assert_true(length(size) == 1, "size must be a scalar.")
#   assert_true(is.finite(size), "size must be finite.")
#   assert_true(size >= 0, "size must be nonnegative.")
#
#   size_round <- round(size)
#   if (abs(size - size_round) > 1e-8) {
#     warning("rmultinom1: size is not integer like; rounding to nearest integer.")
#   }
#   n <- as.double(size_round)
#
#   if (n == 0) return(rep(0, K))
#
#   if (n <= max_chunk) {
#     return(as.numeric(rmultinom(1, size = as.integer(n), prob = prob)[, 1]))
#   }
#
#   n_full <- floor(n / max_chunk)
#   rem <- n - n_full * max_chunk
#
#   out <- numeric(K)
#   if (n_full > 0) {
#     for (i in seq_len(n_full)) {
#       out <- out + as.numeric(rmultinom(1, size = as.integer(max_chunk), prob = prob)[, 1])
#     }
#   }
#   if (rem > 0) {
#     out <- out + as.numeric(rmultinom(1, size = as.integer(rem), prob = prob)[, 1])
#   }
#   out
# }


set_seed <- function(seed) {
  if (!is.null(seed) && !is.na(seed)) set.seed(as.integer(seed))
}

