# 00_write_scenario_grid.R
# Writes a starter scenario grid to config/scenario_grid.csv

dir.create("config", showWarnings = FALSE, recursive = TRUE)

grid <- data.frame(
    scenario_id = c("S0"),
    T = c(100),
    n1 = c(10),
    L = c(2),
    branching = c("4"),
    refine_prob = c(1.0),
    graph_type = c("grid"),
    tau_phi = c(4.0),
    r_size = c(20),
    alpha_dir = c(50),
    gamma_lambda = c(0.98),
    p = c(2),
    beta0 = c(-7.0),
    beta_vals = c("0.25, -0.15"),
    x_sd = c(1.0),
    loge_mean = c(9.5),
    loge_sd = c(0.3),
    zi_prob = c(0.00),
    comp_mode = c("static"),
    cp_time = c(0),
    cp_strength = c(1.0),
    q_drift_sd = c(0.00),
    n_rep = c(30),
    seed_base = c(1000),
    stringsAsFactors = FALSE
)

write.csv(grid, file = "config/scenario_grid_test.csv",
          row.names = FALSE, quote = TRUE)
message("Wrote config/scenario_grid.csv")
