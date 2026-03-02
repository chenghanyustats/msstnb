
# hyper.R
# Default hyperparameter constructors for the MSSTNB Algorithm 2 sampler.

default_hyper <- function(p,
                          # regression prior (21)
                          m0 = 0, s0 = 10,
                          m_beta = NULL, V_beta = NULL,
                          # ICAR precision prior (31)
                          a_phi = 1, b_phi = 1,
                          # lambda initial prior (18)
                          a0 = 20, b0 = 20,
                          # Dirichlet split prior (20)
                          alpha_dir = 50,
                          # r prior (23)
                          a_r = 2, b_r = 2,
                          # discount priors (24)
                          a_gamma = 20, b_gamma = 20,
                          a_delta = 20, b_delta = 20,
                          # initial discounts
                          gamma_init = 0.98, delta_init = 0.98) {
  if (is.null(m_beta)) m_beta <- rep(0, p)
  if (is.null(V_beta)) V_beta <- diag(p) * 100
  list(
    beta = list(m0 = m0, s0 = s0, m_beta = m_beta, V_beta = V_beta),
    tau_phi = list(a_phi = a_phi, b_phi = b_phi),
    lambda = list(a0 = a0, b0 = b0),
    splits = list(alpha_dir = alpha_dir),
    r = list(a_r = a_r, b_r = b_r),
    discount = list(a_gamma = a_gamma, b_gamma = b_gamma,
                    a_delta = a_delta, b_delta = b_delta,
                    gamma_init = gamma_init, delta_init = delta_init)
  )
}
