# 00_write_scenario_grid.R
# Writes a starter scenario grid to config/scenario_grid.csv

dir.create("config", showWarnings = FALSE, recursive = TRUE)

create_scenario <- function(scenario_id = c("S0"),
                            TT = c(50),
                            n1 = c(10),
                            L = c(2),
                            branching = c("4"),
                            refine_prob = c(1),
                            graph_type = c("grid"),
                            tau_phi = c(2),
                            r = c(10),
                            alpha_dir = c(10),
                            gamma_lambda = c(0.95),
                            lambda_a0 = c(10),
                            lambda_b0 = c(10),
                            p = c(2),
                            beta0 = c(-1),
                            beta_vals = c("0.5, 1"),
                            x_sd = c(1),
                            loge_mean = c(5),
                            loge_sd = c(1),
                            zi_prob = c(0),
                            comp_mode = c("static"),
                            cp_time = c(0),
                            cp_strength = c(1),
                            q_drift_sd = c(0),
                            q_conc = c(10),
                            n_rep = c(30),
                            seed_base = c(1000),
                            stringsAsFactors = FALSE) {
  grid <- data.frame(
    scenario_id = scenario_id, TT = TT, n1 = n1, L = L, branching = branching,
    refine_prob = refine_prob, graph_type = graph_type, tau_phi = tau_phi,
    r = r, alpha_dir = alpha_dir, gamma_lambda = gamma_lambda,
    lambda_a0 = lambda_a0, lambda_b0 = lambda_b0, p = p,
    beta0 = beta0, beta_vals = beta_vals, x_sd = x_sd,
    loge_mean = loge_mean, loge_sd = loge_sd, zi_prob = zi_prob,
    comp_mode = comp_mode, cp_time = cp_time, cp_strength = cp_strength,
    q_drift_sd = q_drift_sd, q_conc = q_conc, n_rep = n_rep,
    seed_base = seed_base, stringsAsFactors = stringsAsFactors
  )
  grid
}

s0 <- create_scenario()
s1_t_small <- create_scenario(scenario_id = "S1_t_small", TT = 10)
s1_t_large <- create_scenario(scenario_id = "S1_t_large", TT = 100)
s2_n1_small <- create_scenario(scenario_id = "S2_n1_small", n1 = 5)
s2_n1_large <- create_scenario(scenario_id = "S2_n1_large", n1 = 30)
s3_overdisp_strong <- create_scenario(scenario_id = "S3_overdisp_strong", r = 1)
s3_overdisp_pois <- create_scenario(scenario_id = "S3_overdisp_pois", r = 100)
s4_alpha_spiky <- create_scenario(scenario_id = "S4_alpha_spiky", alpha_dir = 1)
s4_alpha_smooth <- create_scenario(scenario_id = "S4_alpha_smooth", alpha_dir = 100)
s5_icar_smooth <- create_scenario(scenario_id = "S5_icar_smooth", tau_phi = 100)
s5_icar_heter <- create_scenario(scenario_id = "S5_icar_heter", tau_phi = 0.5)
s6_zi <- create_scenario(scenario_id = "S6_zi", zi_prob = 0.2)
s6_cp <- create_scenario(scenario_id = "S6_cp", cp_time = 25, cp_strength = 2,
                         comp_mode = c("changepoint"))
s7_mix <- create_scenario(scenario_id = "S7_mix", refine_prob = 0.8)

grid <- rbind(s0, s1_t_small, s1_t_large, s2_n1_small, s2_n1_large,
              s3_overdisp_strong, s3_overdisp_pois,
              s4_alpha_spiky, s4_alpha_smooth,
              s5_icar_smooth, s5_icar_heter, s6_zi, s6_cp, s7_mix)
write.csv(grid, file = "config/scenario_grid.csv", row.names = FALSE, quote = TRUE)
message("Wrote config/scenario_grid.csv")
