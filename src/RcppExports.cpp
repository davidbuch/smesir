// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// solve_infections
NumericVector solve_infections(const NumericVector beta, const double gamma, const int t_1, const double initial_impulse, const int N);
RcppExport SEXP _smesir_solve_infections(SEXP betaSEXP, SEXP gammaSEXP, SEXP t_1SEXP, SEXP initial_impulseSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int >::type t_1(t_1SEXP);
    Rcpp::traits::input_parameter< const double >::type initial_impulse(initial_impulseSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_infections(beta, gamma, t_1, initial_impulse, N));
    return rcpp_result_gen;
END_RCPP
}
// solve_events
NumericVector solve_events(const NumericVector nu, const NumericVector psi);
RcppExport SEXP _smesir_solve_events(SEXP nuSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_events(nu, psi));
    return rcpp_result_gen;
END_RCPP
}
// log_poisd
double log_poisd(const NumericVector y, const NumericVector lambda);
RcppExport SEXP _smesir_log_poisd(SEXP ySEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_poisd(y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// log_llh
double log_llh(const NumericMatrix B, const NumericMatrix Y, const double gamma, const IntegerVector T_1, const NumericVector Initial_Impulse, const IntegerVector N, const NumericVector psi, const NumericMatrix frailty);
RcppExport SEXP _smesir_log_llh(SEXP BSEXP, SEXP YSEXP, SEXP gammaSEXP, SEXP T_1SEXP, SEXP Initial_ImpulseSEXP, SEXP NSEXP, SEXP psiSEXP, SEXP frailtySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type B(BSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type T_1(T_1SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Initial_Impulse(Initial_ImpulseSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type frailty(frailtySEXP);
    rcpp_result_gen = Rcpp::wrap(log_llh(B, Y, gamma, T_1, Initial_Impulse, N, psi, frailty));
    return rcpp_result_gen;
END_RCPP
}
// smesir_mcmc
List smesir_mcmc(const NumericMatrix Y, const List Design_Matrices, const arma::vec& lambda, const arma::vec& V0param, const NumericMatrix IGSR, const double gamma, const IntegerVector T_1, const NumericVector Initial_Impulse, const double dispersion, const IntegerVector N, const NumericVector psi, const double tempering_ratio, const int ncycles, const int samps_per_cycle, const int nchain, const int iter, const int warmup, const int thin, const bool sr_style, const bool quiet);
RcppExport SEXP _smesir_smesir_mcmc(SEXP YSEXP, SEXP Design_MatricesSEXP, SEXP lambdaSEXP, SEXP V0paramSEXP, SEXP IGSRSEXP, SEXP gammaSEXP, SEXP T_1SEXP, SEXP Initial_ImpulseSEXP, SEXP dispersionSEXP, SEXP NSEXP, SEXP psiSEXP, SEXP tempering_ratioSEXP, SEXP ncyclesSEXP, SEXP samps_per_cycleSEXP, SEXP nchainSEXP, SEXP iterSEXP, SEXP warmupSEXP, SEXP thinSEXP, SEXP sr_styleSEXP, SEXP quietSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const List >::type Design_Matrices(Design_MatricesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type V0param(V0paramSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type IGSR(IGSRSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type T_1(T_1SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Initial_Impulse(Initial_ImpulseSEXP);
    Rcpp::traits::input_parameter< const double >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type tempering_ratio(tempering_ratioSEXP);
    Rcpp::traits::input_parameter< const int >::type ncycles(ncyclesSEXP);
    Rcpp::traits::input_parameter< const int >::type samps_per_cycle(samps_per_cycleSEXP);
    Rcpp::traits::input_parameter< const int >::type nchain(nchainSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type sr_style(sr_styleSEXP);
    Rcpp::traits::input_parameter< const bool >::type quiet(quietSEXP);
    rcpp_result_gen = Rcpp::wrap(smesir_mcmc(Y, Design_Matrices, lambda, V0param, IGSR, gamma, T_1, Initial_Impulse, dispersion, N, psi, tempering_ratio, ncycles, samps_per_cycle, nchain, iter, warmup, thin, sr_style, quiet));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_smesir_solve_infections", (DL_FUNC) &_smesir_solve_infections, 5},
    {"_smesir_solve_events", (DL_FUNC) &_smesir_solve_events, 2},
    {"_smesir_log_poisd", (DL_FUNC) &_smesir_log_poisd, 2},
    {"_smesir_log_llh", (DL_FUNC) &_smesir_log_llh, 8},
    {"_smesir_smesir_mcmc", (DL_FUNC) &_smesir_smesir_mcmc, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_smesir(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
