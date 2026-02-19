MSSTNB simulation data generation

This folder contains R code to generate coherent multiscale count data under the MSSTNB generative assumptions.

Core ingredients
1. Coarsest level counts follow a negative binomial model with NB2 mean size parameterization, implemented via a Poisson-Gamma augmentation.
2. Mean structure uses exposure offsets, covariates, and a constrained ICAR spatial effect at level 1.
3. Multiscale coherence is enforced by multinomial splits down a rooted tree.
4. Split probabilities can be static, drifting, or have a change point.
5. Optional misspecification includes zero inflation at the coarsest level and change points in split compositions.

Folder map
config
  scenario_grid.csv  Scenario definitions. Edit this first.
R
  generate_dataset.R  Main generator that returns a list with data and truth
  make_graph.R        Build adjacency and ICAR precision H
  sim_icar.R          Simulate constrained ICAR spatial effects
  make_tree.R         Build a balanced or mixed resolution tree
  sim_lambda_tilde.R  Simulate dynamic residual risk
  sim_splits.R        Simulate Dirichlet baseline compositions and split probabilities
  utils_io.R          Small helpers
scripts
  00_write_scenario_grid.R  Optional helper to write a starter scenario grid
  01_generate_all_data.R    Driver: loops over scenarios and replicates and writes RDS files

How to run
1. Open config/scenario_grid.csv and set scenario values
2. From the project root, run
   Rscript scripts/01_generate_all_data.R

Outputs
data_sim
  <scenario_id>
    rep0001.rds
    rep0002.rds
    ...

Each rep file contains
$list$data
  y_levels  list of count matrices by level
  e         exposure matrix at level 1
  x         covariate array at level 1
  tree      tree object including children mapping
  graph     adjacency and H
$list$truth
  beta0, beta, phi, tau_phi, r1, kappa, lambda_tilde
  omega, q, alpha_dir
  scenario  the scenario row used
