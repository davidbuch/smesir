# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

solve_infections <- function(beta, gamma, T_1, IIP, N, vaccinations) {
    .Call(`_smesir_solve_infections`, beta, gamma, T_1, IIP, N, vaccinations)
}

solve_susceptible <- function(beta, gamma, T_1, IIP, N, vaccinations) {
    .Call(`_smesir_solve_susceptible`, beta, gamma, T_1, IIP, N, vaccinations)
}

solve_events <- function(nu, psi) {
    .Call(`_smesir_solve_events`, nu, psi)
}

smesir_mcmc <- function(Y, Design_Matrices, TildeOff, vaccinations, lambda, V0param, IGSR, gamma, T_1, expected_iip, expected_disp, N, psi, discount_period_length, discount_period_disp, ncycles, samps_per_cycle, nchain, iter, warmup, thin, sr_style, quiet) {
    .Call(`_smesir_smesir_mcmc`, Y, Design_Matrices, TildeOff, vaccinations, lambda, V0param, IGSR, gamma, T_1, expected_iip, expected_disp, N, psi, discount_period_length, discount_period_disp, ncycles, samps_per_cycle, nchain, iter, warmup, thin, sr_style, quiet)
}

