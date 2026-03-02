# sampler_algo2.R
# Blocked Gibbs sampler implementing Algorithm 2.

flatten_internal_nodes <- function(tree) {
    # returns a data.frame with columns l, j, K
    L <- tree$L
    rows <- list()
    idx <- 1
    for (l in seq_len(L - 1)) {
        n_parent <- tree$n[l]
        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K > 1) {
                rows[[idx]] <- data.frame(l = l, j = j, K = K,
                                          stringsAsFactors = FALSE)
                idx <- idx + 1
            }
        }
    }
    if (length(rows) == 0) return(data.frame(l = integer(), j = integer(),
                                             K = integer()))
    do.call(rbind, rows)
}

delta_to_vec <- function(tree, delta) {
    tab <- flatten_internal_nodes(tree)
    v <- numeric(nrow(tab))
    for (i in seq_len(nrow(tab))) {
        l <- tab$l[i]; j <- tab$j[i]
        v[i] <- delta[[l]][[j]]
    }
    v
}

vec_to_delta <- function(tree, v, delta_template) {
    tab <- flatten_internal_nodes(tree)
    out <- delta_template
    for (i in seq_len(nrow(tab))) {
        l <- tab$l[i]; j <- tab$j[i]
        out[[l]][[j]] <- v[i]
    }
    out
}

msstnb_algo2_mcmc <- function(data, hyper, mcmc,
                              init = NULL,
                              seed = NULL,
                              verbose = TRUE) {
    set_seed(seed)

    # unpack data
    y_levels <- data$y_levels
    y1 <- y_levels[[1]]
    e1 <- data$e
    x1 <- data$x
    tree <- data$tree
    TT <- nrow(y1); n1 <- ncol(y1); p <- dim(x1)[3]

    # Cached design matrix for ESS blocks (column-major: time within region)
    Xmat <- if (p > 0) matrix(x1, nrow = TT * n1, ncol = p) else NULL

    # Cached prior for (beta0, beta) ESS
    beta_prior_prep <- prep_beta_prior(
        p = p,
        m0 = hyper$beta$m0, s0 = hyper$beta$s0,
        m_beta = hyper$beta$m_beta, V_beta = hyper$beta$V_beta
    )


    # defaults
    n_iter <- mcmc$n_iter %||% 2000L
    burn   <- mcmc$burn   %||% floor(n_iter/2)
    thin   <- mcmc$thin   %||% 1L
    prop_sd_r <- mcmc$prop_sd_r %||% 0.2
    update_discounts_every <- mcmc$update_discounts_every %||% 1L
    store_phi <- isTRUE(mcmc$store_phi)
    store_lambda <- isTRUE(mcmc$store_lambda)
    store_gamma <- isTRUE(mcmc$store_gamma)
    store_delta <- isTRUE(mcmc$store_delta)
    store_r <- isTRUE(mcmc$store_r)

    assert_true(burn < n_iter, "burn must be < n_iter.")
    assert_true(thin >= 1, "thin must be >= 1.")
    n_save <- floor((n_iter - burn)/thin)

    # init
    if (is.null(init)) {
      state <- init_state(data, hyper, seed = seed)
    } else {
      state <- init
    }

    # allocate draws
    beta0_draw <- numeric(n_save)
    beta_draw  <- matrix(NA_real_, nrow = n_save, ncol = p)
    tau_phi_draw <- numeric(n_save)
    if (store_phi) {
        phi_draw <- matrix(NA_real_, nrow = n_save, ncol = n1)
    } else {
        phi_draw <- NULL
    }
    if (store_lambda) {
        lambda_draw <- array(NA_real_, dim = c(n_save, TT + 1, n1))
    } else {
        lambda_draw <- NULL
    }
    if (store_gamma) {
        gamma_draw <- matrix(NA_real_, nrow = n_save, ncol = n1)
    } else {
        gamma_draw <- NULL
    }
    if (store_r) {
        r_draw <- matrix(NA_real_, nrow = n_save, ncol = n1)
    } else {
        r_draw <- NULL
    }

    tab_nodes <- flatten_internal_nodes(tree)
    if (store_delta) {
        delta_draw <- matrix(NA_real_, nrow = n_save, ncol = nrow(tab_nodes))
    } else {
        delta_draw <- NULL
    }

    # diagnostics
    acc_r <- numeric(n_iter)
    ess_beta_steps <- integer(n_iter)
    ess_phi_steps  <- integer(n_iter)

    save_idx <- 0L

    for (iter in seq_len(n_iter)) {

        # 2) Update kappa via (27)
        kap_out <- update_kappa(y1, e1, x1,
                                beta0 = state$beta0, beta = state$beta,
                                phi = state$phi,
                                lambda_tilde_0T = state$lambda_tilde_0T,
                                r1 = state$r1)
        state$kappa <- kap_out$kappa

        # 3) Update (beta0, beta) via ESS targeting (32)
        bet_out <- update_beta_ess(beta0 = state$beta0, beta = state$beta,
                                   phi = state$phi,
                                   lambda_tilde_0T = state$lambda_tilde_0T,
                                   kappa = state$kappa,
                                   y1 = y1, e1 = e1, x1 = x1,
                                   m0 = hyper$beta$m0, s0 = hyper$beta$s0,
                                   m_beta = hyper$beta$m_beta,
                                   V_beta = hyper$beta$V_beta,
                                   Xmat = Xmat, prior_prep = beta_prior_prep)
        state$beta0 <- bet_out$beta0
        state$beta  <- bet_out$beta
        ess_beta_steps[iter] <- bet_out$ess_steps

        # 4) Update phi via ICAR-aware ESS targeting (33)
        phi_out <- update_phi_ess(phi = state$phi, tau_phi = state$tau_phi,
                                  H = state$H,
                                  beta0 = state$beta0, beta = state$beta,
                                  lambda_tilde_0T = state$lambda_tilde_0T,
                                  kappa = state$kappa, y1 = y1, e1 = e1, x1 = x1,
                                  B = state$B, Xmat = Xmat)
        state$phi <- phi_out$phi
        state$B   <- phi_out$B
        ess_phi_steps[iter] <- phi_out$ess_steps

        # 5) Update tau_phi via (31)
        state$tau_phi <- update_tau_phi(state$phi, state$H,
                                        a_phi = hyper$tau_phi$a_phi,
                                        b_phi = hyper$tau_phi$b_phi)

        # 6-7) FFBS for lambda_tilde via (29) and (34)
        lam_out <- ffbs_lambda(y1, e1, x1,
                               beta0 = state$beta0, beta = state$beta,
                               phi = state$phi,
                               kappa = state$kappa, gamma = state$gamma,
                               a0 = hyper$lambda$a0, b0 = hyper$lambda$b0)
        state$lambda_tilde_0T <- lam_out$lambda_tilde_0T

        # 6-7) FFBS for omega via (30) and (35)-(37)
        om_out <- ffbs_omega(tree = tree, y_levels = y_levels,
                             alpha_dir = hyper$splits$alpha_dir,
                             delta = state$delta, q0 = mcmc$q0_splits %||% NULL)
        state$omega <- om_out$omega_0T

        # 8) Update dispersion r via MH on log r (Appendix A.10)
        r_out <- update_r_mh(state$r1, state$kappa,
                             a_r = hyper$r$a_r, b_r = hyper$r$b_r,
                             prop_sd = prop_sd_r)
        state$r1 <- r_out$r1
        acc_r[iter] <- r_out$accept

        # 9) Update discounts gamma and delta (Section 8.6)
        if (update_discounts_every > 0 &&
            (iter %% update_discounts_every == 0)) {
          # gamma: use xi from lambda filtering
          # (depends on current beta, phi, kappa)
          state$gamma <- update_gamma_grid(
              y1 = y1, xi = lam_out$xi,
              a0 = hyper$lambda$a0, b0 = hyper$lambda$b0,
              a_prior = hyper$discount$a_gamma,
              b_prior = hyper$discount$b_gamma,
              grid = mcmc$gamma_grid %||% seq(0.05, 0.999, length.out = 200L),
              time_idx = mcmc$discount_time_idx %||% NULL)

          state$delta <- update_delta_grid(
              tree = tree, y_levels = y_levels,
              alpha_dir = hyper$splits$alpha_dir,
              q0 = mcmc$q0_splits %||% NULL,
              delta = state$delta,
              a_prior = hyper$discount$a_delta,
              b_prior = hyper$discount$b_delta,
              grid = mcmc$delta_grid %||% seq(0.05, 0.999, length.out = 200L),
              time_idx = mcmc$discount_time_idx %||% NULL)
        }

        # store
        if (iter > burn && ((iter - burn) %% thin == 0)) {
            save_idx <- save_idx + 1L
            beta0_draw[save_idx] <- state$beta0
            beta_draw[save_idx, ] <- state$beta
            tau_phi_draw[save_idx] <- state$tau_phi
            if (store_phi) phi_draw[save_idx, ] <- state$phi
            if (store_lambda) lambda_draw[save_idx, , ] <- state$lambda_tilde_0T
            if (store_gamma) gamma_draw[save_idx, ] <- state$gamma
            if (store_r) r_draw[save_idx, ] <- state$r1
            if (store_delta) {
                delta_draw[save_idx, ] <- delta_to_vec(tree, state$delta)
            }
        }

        if (verbose && (iter %% max(1L, floor(n_iter/10)) == 0)) {
            cat(sprintf("Iter %d/%d | acc_r=%.2f | ESS(beta)=%.0f | ESS(phi)=%.0f\n",
                        iter, n_iter, mean(acc_r[pmax(1, iter-99):iter]),
                        mean(ess_beta_steps[pmax(1, iter-99):iter]),
                        mean(ess_phi_steps[pmax(1, iter-99):iter])))
        }
    }

    draws <- list(
        beta0 = beta0_draw,
        beta  = beta_draw,
        tau_phi = tau_phi_draw,
        phi = phi_draw,
        lambda_tilde_0T = lambda_draw,
        gamma = gamma_draw,
        r1 = r_draw,
        delta = delta_draw,
        delta_nodes = tab_nodes,
        diagnostics = list(acc_r = acc_r,
                           ess_beta_steps = ess_beta_steps,
                           ess_phi_steps = ess_phi_steps),
        mcmc = list(n_iter = n_iter, burn = burn, thin = thin)
    )
    list(draws = draws, last_state = state)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
