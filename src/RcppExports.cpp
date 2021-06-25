// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CumSum
NumericMatrix CumSum(NumericMatrix m);
RcppExport SEXP _iPath_CumSum(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(CumSum(m));
    return rcpp_result_gen;
END_RCPP
}
// caliES2
double caliES2(const std::vector<double>& ranks, const std::vector<int>& pos);
RcppExport SEXP _iPath_caliES2(SEXP ranksSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ranks(ranksSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(caliES2(ranks, pos));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iPath_CumSum", (DL_FUNC) &_iPath_CumSum, 1},
    {"_iPath_caliES2", (DL_FUNC) &_iPath_caliES2, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_iPath(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}