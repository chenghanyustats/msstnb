
# ffbs_omega.R
# Forward filtering (30) and backward sampling (35)-(37) for omega split trajectories.

# Build prior c0 for each internal node (l,j): alpha_dir * q0
init_c0_splits <- function(tree, alpha_dir = 50, q0 = NULL) {
    L <- tree$L
    if (L <= 1) return(list())
    c0 <- vector("list", L - 1)
    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        c0_l <- vector("list", n_parent)
      for (j in seq_len(n_parent)) {
          ch <- tree$children[[l]][[j]]
          K <- length(ch)
          if (K <= 1) {
              c0_l[[j]] <- NULL
          } else {
              if (is.null(q0)) {
                  q <- rep(1 / K, K)
          } else {
              # Expect q0[[l]][[j]] length K
              q <- q0[[l]][[j]]
              assert_true(length(q) == K, "q0 structure mismatch.")
              q <- pmax(q, 1e-12); q <- q / sum(q)
          }
              c0_l[[j]] <- alpha_dir * q
          }
      }
          c0[[l]] <- c0_l
    }
    c0
}

# delta structure: list over levels, list over parents
init_delta_splits <- function(tree, delta0 = 0.98) {
    L <- tree$L
    if (L <= 1) return(list())
    delta <- vector("list", L - 1)
    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        dl <- vector("list", n_parent)
        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K <= 1) dl[[j]] <- NULL else dl[[j]] <- as.numeric(delta0)
        }
        delta[[l]] <- dl
    }
    delta
}

filter_omega <- function(tree, y_levels, delta, c0) {
    TT <- nrow(y_levels[[1]])
    L <- tree$L
    assert_true(length(y_levels) == L, "y_levels length must equal tree$L.")
    assert_true(length(delta) == L - 1, "delta must be list of length L-1.")
    assert_true(length(c0) == L - 1, "c0 must be list of length L-1.")

    c_time <- vector("list", L - 1)
    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        c_l <- vector("list", n_parent)
        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K <= 1) { c_l[[j]] <- NULL; next }
            del <- delta[[l]][[j]]
            assert_true(!is.null(del) && del > 0 && del <= 1,
                        "delta must be in (0,1].")

            # child counts y_{t, D(l,j)}
            y_child <- y_levels[[l + 1]][, ch, drop = FALSE]  # T x K
            # initialize
            cmat <- matrix(NA_real_, nrow = TT + 1, ncol = K)
            cmat[1, ] <- c0[[l]][[j]]
            for (t in seq_len(TT)) {
                cmat[t + 1, ] <- del * cmat[t, ] + y_child[t, ]
            }
            c_l[[j]] <- cmat
        }
    c_time[[l]] <- c_l
    }
    c_time
}

backward_sample_omega_one <- function(cmat, del) {
    # cmat: (T+1) x K
    TT <- nrow(cmat) - 1
    K <- ncol(cmat)
    omega <- matrix(NA_real_, nrow = TT + 1, ncol = K)
    # sample omega_T
    omega[TT + 1, ] <- rdirichlet1(pmax(cmat[TT + 1, ], 1e-12))
    for (t in TT:1) {
        if (del >= 1) {
          omega[t, ] <- omega[t + 1, ]
          next
        }
        c_prev <- pmax(cmat[t, ], 1e-12)
        C_prev <- sum(c_prev)
        aS <- del * C_prev
        bS <- (1 - del) * C_prev
        # Beta params must be > 0; if bS extremely small, force omega[t] = omega[t+1]
        if (bS < 1e-12) {
          omega[t, ] <- omega[t + 1, ]
          next
        }
        S <- rbeta(1, shape1 = pmax(aS, 1e-8), shape2 = pmax(bS, 1e-8))
        omega_tilde <- rdirichlet1((1 - del) * c_prev)
        omega[t, ] <- (1 - S) * omega_tilde + S * omega[t + 1, ]
        # renormalize for safety
        omega[t, ] <- pmax(omega[t, ], 0)
        omega[t, ] <- omega[t, ] / sum(omega[t, ])
    }
    omega
}

backward_sample_omega <- function(tree, c_time, delta) {
    L <- tree$L
    if (L <= 1) return(list())
    omega_time <- vector("list", L - 1)
    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        omega_l <- vector("list", n_parent)
        for (j in seq_len(n_parent)) {
            cmat <- c_time[[l]][[j]]
            if (is.null(cmat)) { omega_l[[j]] <- NULL; next }
            del <- delta[[l]][[j]]
            omega_l[[j]] <- backward_sample_omega_one(cmat, del)
        }
        omega_time[[l]] <- omega_l
    }
    omega_time
}

ffbs_omega <- function(tree, y_levels, alpha_dir = 50, delta, q0 = NULL,
                       seed = NULL) {
    set_seed(seed)
    c0 <- init_c0_splits(tree, alpha_dir = alpha_dir, q0 = q0)
    c_time <- filter_omega(tree, y_levels, delta = delta, c0 = c0)
    omega_time <- backward_sample_omega(tree, c_time, delta = delta)
    list(omega_0T = omega_time, c = c_time, c0 = c0)
}
