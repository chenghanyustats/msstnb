# sampler_algo2_reordered.R
# Blocked Gibbs sampler for MSSTNB (Algorithm 2), with discounts updated before FFBS.
#
# Refactor goals:
# - package-ready organization
# - cache heavy objects once (Xmat, chol(B'HB), delta mapping)
# - reduce inner-loop overhead

msstnb_algo2_mcmc <- function(data, hyper, mcmc, init = NULL, seed = NULL,
                              verbose = TRUE) {
    set_seed(seed)

    assert_true(is.list(data), "data must be a list.")
    assert_true(is.list(hyper), "hyper must be a list.")
    assert_true(is.list(mcmc), "mcmc must be a list.")

    # Required MCMC controls
    n_iter <- as.integer(mcmc$n_iter %||% 2000L)
    burn   <- as.integer(mcmc$burn %||% floor(n_iter / 2))
    thin   <- as.integer(mcmc$thin %||% 1L)

    assert_true(n_iter >= 1L, "n_iter must be >= 1.")
    assert_true(burn >= 0L && burn < n_iter, "burn must be in [0, n_iter).");
    assert_true(thin >= 1L, "thin must be >= 1.")

    # Storage toggles
    store_phi    <- isTRUE(mcmc$store_phi %||% FALSE)
    store_lambda <- isTRUE(mcmc$store_lambda %||% FALSE)
    store_r      <- isTRUE(mcmc$store_r %||% FALSE)
    store_gamma  <- isTRUE(mcmc$store_gamma %||% TRUE)
    store_delta  <- isTRUE(mcmc$store_delta %||% TRUE)

    # Discount update schedule
    update_discounts_every <- as.integer(mcmc$update_discounts_every %||% 1L)
    assert_true(update_discounts_every >= 1L,
                "update_discounts_every must be >= 1.")

    # r proposal
    prop_sd_r_init <- as.numeric(mcmc$prop_sd_r %||% 0.2)
    assert_true(is.finite(prop_sd_r_init) && prop_sd_r_init > 0,
                "prop_sd_r must be positive.")

    # Adaptive MH for r
    adapt_r <- isTRUE(mcmc$adapt_r %||% FALSE)
    adapt_r_window <- as.integer(mcmc$adapt_r_window %||% 100L)
    adapt_r_target <- as.numeric(mcmc$adapt_r_target %||% 0.30)
    adapt_r_gamma0 <- as.numeric(mcmc$adapt_r_gamma0 %||% 0.05)
    adapt_r_t0 <- as.numeric(mcmc$adapt_r_t0 %||% 10)
    adapt_r_power <- as.numeric(mcmc$adapt_r_power %||% 0.5)
    adapt_r_min_sd <- as.numeric(mcmc$adapt_r_min_sd %||% 0.01)
    adapt_r_max_sd <- as.numeric(mcmc$adapt_r_max_sd %||% 2.0)
    adapt_r_until <- as.integer(mcmc$adapt_r_until %||% burn)

    if (adapt_r) {
        assert_true(adapt_r_window >= 5L, "adapt_r_window must be >= 5.")
        assert_true(adapt_r_target > 0 && adapt_r_target < 1,
                    "adapt_r_target must be in (0, 1).")
        assert_true(adapt_r_until >= 1L, "adapt_r_until must be >= 1.")
    }

    # Data
    y_levels <- data$y_levels
    y1 <- y_levels[[1]]
    e1 <- data$e
    x1 <- data$x
    tree <- data$tree

    TT <- nrow(y1)
    n1 <- ncol(y1)
    p  <- dim(x1)[3]

    # Build caches
    cache <- build_msstnb_cache(data = data, hyper = hyper, mcmc = mcmc)

    # Initialize state
    if (is.null(init)) {
        state <- init_state(data = data, hyper = hyper, seed = seed,
                            cache = cache)
    } else {
        state <- init
    }

    # Ensure required components exist
    needed <- c("beta0", "beta", "phi", "tau_phi", "lambda_tilde_0T",
                "kappa", "r1", "gamma", "delta", "omega")
    miss <- setdiff(needed, names(state))
    assert_true(length(miss) == 0L, paste("init is missing:",
                                          paste(miss, collapse = ", ")))

    # Preallocate storage
    n_save <- as.integer(floor((n_iter - burn) / thin))
    assert_true(n_save >= 1L, "No samples would be saved. Check burn/thin.")

    draws_beta0 <- numeric(n_save)
    draws_beta  <- matrix(NA_real_, nrow = n_save, ncol = p)

    draws_phi <- if (store_phi) matrix(NA_real_, nrow = n_save,
                                       ncol = n1) else NULL
    draws_lambda <- if (store_lambda) array(NA_real_,
                                            dim = c(n_save, TT + 1L, n1)) else NULL
    draws_r <- if (store_r) matrix(NA_real_, nrow = n_save,
                                   ncol = n1) else NULL
    draws_gamma <- if (store_gamma) matrix(NA_real_, nrow = n_save,
                                           ncol = n1) else NULL

    delta_map <- cache$delta_map
    n_delta_nodes <- nrow(delta_map$tab)
    draws_delta <- if (store_delta) matrix(NA_real_, nrow = n_save,
                                           ncol = n_delta_nodes) else NULL

    # Diagnostics
    ess_steps_beta <- integer(n_iter)
    ess_steps_phi  <- integer(n_iter)
    accept_r <- numeric(n_iter)

    accept_r_window <- rep(NA_real_, n_iter)
    prop_sd_r_trace <- rep(NA_real_, n_iter)
    acc_r_count_postburn <- integer(n1)
    n_postburn <- 0L

    prop_sd_r_curr <- prop_sd_r_init

    save_idx <- 0L

    if (verbose) {
        message("MSSTNB Algorithm 2 MCMC")
        message("  TT = ", TT, ", n1 = ", n1, ", p = ", p)
        message("  n_iter = ", n_iter, ", burn = ", burn, ", thin = ",
                thin, ", n_save = ", n_save)
    }

    for (iter in seq_len(n_iter)) {
        # 1) kappa | ...
        kappa_out <- update_kappa(
            y1 = y1,
            e1 = e1,
            x1 = x1,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            lambda_tilde_0T = state$lambda_tilde_0T,
            r1 = state$r1,
            Xmat = cache$Xmat
        )
        state$kappa <- kappa_out$kappa

        # 2) beta0, beta | ... via ESS
        beta_out <- update_beta_ess(
            beta0 = state$beta0,
            beta  = state$beta,
            phi = state$phi,
            lambda_tilde_0T = state$lambda_tilde_0T,
            kappa = state$kappa,
            y1 = y1,
            e1 = e1,
            x1 = x1,
            m0 = hyper$beta$m0,
            s0 = hyper$beta$s0,
            m_beta = hyper$beta$m_beta,
            V_beta = hyper$beta$V_beta,
            Xmat = cache$Xmat,
            prior_prep = cache$beta_prior_prep
        )
        state$beta0 <- beta_out$beta0
        state$beta  <- beta_out$beta
        ess_steps_beta[iter] <- beta_out$ess_steps

        # 3) phi | ... via ESS (ICAR)
        phi_out <- update_phi_ess(
            phi = state$phi,
            tau_phi = state$tau_phi,
            H = cache$H,
            beta0 = state$beta0,
            beta = state$beta,
            lambda_tilde_0T = state$lambda_tilde_0T,
            kappa = state$kappa,
            y1 = y1,
            e1 = e1,
            x1 = x1,
            B = cache$B,
            Xmat = cache$Xmat,
            chol_BHB = cache$chol_BHB,
            rep_idx = cache$rep_idx
        )
        state$phi <- phi_out$phi
        ess_steps_phi[iter] <- phi_out$ess_steps

        # 4) tau_phi | phi
        state$tau_phi <- update_tau_phi(
            phi = state$phi,
            H = cache$H,
            a_phi = hyper$tau_phi$a_phi,
            b_phi = hyper$tau_phi$b_phi
        )

        # 5) r1 | kappa via MH
        r_out <- update_r_mh(
            r1 = state$r1,
            kappa = state$kappa,
            a_r = hyper$r$a_r,
            b_r = hyper$r$b_r,
            prop_sd = prop_sd_r_curr
        )
        state$r1 <- r_out$r1
        accept_r[iter] <- r_out$accept
        prop_sd_r_trace[iter] <- prop_sd_r_curr

        # post-burn acceptance by region
        if (iter > burn) {
            acc_r_count_postburn <- acc_r_count_postburn +
              as.integer(r_out$accept_vec)
            n_postburn <- n_postburn + 1L
        }

        # adapt proposal sd for r using recent M iterations acceptance rate
        if (adapt_r && iter <= adapt_r_until) {
            w0 <- max(1L, iter - adapt_r_window + 1L)
            acc_bar <- mean(accept_r[w0:iter])
            accept_r_window[iter] <- acc_bar
            prop_sd_r_curr <- adapt_prop_sd_r(
                prop_sd = prop_sd_r_curr,
                acc_bar = acc_bar,
                target = adapt_r_target,
                iter = iter,
                gamma0 = adapt_r_gamma0,
                t0 = adapt_r_t0,
                power = adapt_r_power,
                min_sd = adapt_r_min_sd,
                max_sd = adapt_r_max_sd
            )
        }

        # Precompute eta and xi once (used by discounts and lambda FFBS)
        eta <- compute_eta(
            beta0 = state$beta0,
            beta  = state$beta,
            phi   = state$phi,
            TT = TT,
            n1 = n1,
            Xmat = cache$Xmat
        )
        # xi = e * kappa * exp(eta)
        xi <- compute_xi(e1 = e1, kappa = state$kappa, eta = eta)

        # 6) discounts (gamma, delta) before FFBS
        if (iter %% update_discounts_every == 0L) {
            state$gamma <- update_gamma_grid(
                y1 = y1,
                xi = xi,
                a0 = hyper$lambda$a0,
                b0 = hyper$lambda$b0,
                a_prior = hyper$discount$a_gamma,
                b_prior = hyper$discount$b_gamma,
                grid = cache$gamma_grid,
                time_idx = cache$discount_time_idx
            )

            state$delta <- update_delta_grid(
                tree = tree,
                y_levels = y_levels,
                alpha_dir = hyper$splits$alpha_dir,
                q0 = NULL,
                delta = state$delta,
                a_prior = hyper$discount$a_delta,
                b_prior = hyper$discount$b_delta,
                grid = cache$delta_grid,
                time_idx = cache$discount_time_idx
            )
        }

        # 7) lambda FFBS
        lam_out <- ffbs_lambda_from_xi(
            y1 = y1,
            xi = xi,
            gamma = state$gamma,
            a0 = hyper$lambda$a0,
            b0 = hyper$lambda$b0
        )
        state$lambda_tilde_0T <- lam_out$lambda_tilde_0T

        # 8) omega FFBS
        omega_out <- ffbs_omega(
            tree = tree,
            y_levels = y_levels,
            alpha_dir = hyper$splits$alpha_dir,
            delta = state$delta,
            q0 = NULL
        )
        state$omega <- omega_out$omega_0T

        # Save
        if (iter > burn && ((iter - burn) %% thin == 0L)) {
            save_idx <- save_idx + 1L

            draws_beta0[save_idx] <- state$beta0
            if (p > 0L) {
                draws_beta[save_idx, ] <- state$beta
            }

            if (store_phi) {
                draws_phi[save_idx, ] <- state$phi
            }

            if (store_lambda) {
                draws_lambda[save_idx, , ] <- state$lambda_tilde_0T
            }

            if (store_r) {
                draws_r[save_idx, ] <- state$r1
            }

            if (store_gamma) {
                draws_gamma[save_idx, ] <- state$gamma
            }

            if (store_delta) {
                draws_delta[save_idx, ] <- delta_to_vec(state$delta,
                                                        delta_map = delta_map)
            }
        }

        if (verbose && (iter %% (mcmc$progress_every %||% 50L) == 0L)) {
            message(
                sprintf(
                    "iter %d/%d | ESS(beta) %d | ESS(phi) %d | acc(r) %.2f | sd(r) %.3f",
                    iter, n_iter, ess_steps_beta[iter], ess_steps_phi[iter],
                    accept_r[iter], prop_sd_r_trace[iter]
                )
            )
        }
    }

    accept_r_by_region_postburn <- acc_r_count_postburn / max(1L, n_postburn)

    out <- list(
        draws = list(
            beta0 = draws_beta0,
            beta  = draws_beta,
            phi   = draws_phi,
            lambda_tilde_0T = draws_lambda,
            r1 = draws_r,
            gamma = draws_gamma,
            delta = draws_delta,
            delta_map = delta_map$tab
        ),
        diag = list(
            ess_steps_beta = ess_steps_beta,
            ess_steps_phi  = ess_steps_phi,
            accept_r = accept_r,
            acc_r = accept_r,
            accept_r_window = accept_r_window,
            prop_sd_r = prop_sd_r_trace,
            accept_r_by_region_postburn = accept_r_by_region_postburn
        ),
        state_last = state,
        cache = list(
            TT = TT,
            n1 = n1,
            p = p,
            n_save = n_save
        ),
        call = match.call()
    )

    class(out) <- c("msstnb_fit", class(out))
    out
}
