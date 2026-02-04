\
# 01_generate_all_data.R
# Driver script: generate data for all scenarios and replicates and write to disk

message("Starting MSSTNB simulation data generation")

# Source modules
source("R/utils_io.R")
source("R/make_graph.R")
source("R/sim_icar.R")
source("R/make_tree.R")
source("R/sim_lambda_tilde.R")
source("R/sim_splits.R")
source("R/generate_dataset.R")

dir_create("data_sim")

grid_path <- file.path("config", "scenario_grid.csv")
assert_true(file.exists(grid_path), paste("Missing", grid_path))

scen_grid <- read.csv(grid_path, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("scenario_id","T","n1","L","branching","refine_prob","graph_type","tau_phi",
                   "r_size","alpha_dir","gamma_lambda","p","beta0","beta_vals","x_sd","loge_mean",
                   "loge_sd","zi_prob","comp_mode","cp_time","cp_strength","q_drift_sd","n_rep","seed_base")

missing_cols <- setdiff(required_cols, names(scen_grid))
assert_true(length(missing_cols) == 0, paste("Scenario grid missing columns:", paste(missing_cols, collapse = ", ")))

for (s in seq_len(nrow(scen_grid))) {
  scen <- scen_grid[s, , drop = FALSE]
  scen_id <- as.character(scen$scenario_id[1])
  n_rep <- as.integer(scen$n_rep[1])

  out_dir <- file.path("data_sim", scen_id)
  dir_create(out_dir)

  # Save scenario row
  write.csv(scen, file = file.path(out_dir, "scenario_row.csv"), row.names = FALSE, quote = TRUE)

  message("Scenario ", scen_id, " with ", n_rep, " replicates")

  for (r in seq_len(n_rep)) {
    ds <- generate_dataset(scen, rep_id = r)
    out_file <- file.path(out_dir, sprintf("rep%04d.rds", r))
    saveRDS(ds, out_file)
    if (r %% 10 == 0) message("  wrote ", r, " / ", n_rep)
  }

  # Save session info for reproducibility
  sink(file.path(out_dir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

message("Done. Outputs are in data_sim/")
