// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// make_trails_cpp
DataFrame make_trails_cpp(DataFrame field_df, DataFrame particles_df, int width, int height, int steps, double step_length, double min_dist, int seed);
RcppExport SEXP _flow_make_trails_cpp(SEXP field_dfSEXP, SEXP particles_dfSEXP, SEXP widthSEXP, SEXP heightSEXP, SEXP stepsSEXP, SEXP step_lengthSEXP, SEXP min_distSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type field_df(field_dfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type particles_df(particles_dfSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    Rcpp::traits::input_parameter< int >::type height(heightSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type step_length(step_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type min_dist(min_distSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(make_trails_cpp(field_df, particles_df, width, height, steps, step_length, min_dist, seed));
    return rcpp_result_gen;
END_RCPP
}
// make_trail
DataFrame make_trail(double x0, double y0, DataFrame field_df, DataFrame existing_points, double step_length, int max_steps, std::string direction, double dtest);
RcppExport SEXP _flow_make_trail(SEXP x0SEXP, SEXP y0SEXP, SEXP field_dfSEXP, SEXP existing_pointsSEXP, SEXP step_lengthSEXP, SEXP max_stepsSEXP, SEXP directionSEXP, SEXP dtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< DataFrame >::type field_df(field_dfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type existing_points(existing_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type step_length(step_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< std::string >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    rcpp_result_gen = Rcpp::wrap(make_trail(x0, y0, field_df, existing_points, step_length, max_steps, direction, dtest));
    return rcpp_result_gen;
END_RCPP
}
// make_trails_rcpp
DataFrame make_trails_rcpp(DataFrame particles, DataFrame field_df, double step_length, int max_steps, std::string direction, double dtest);
RcppExport SEXP _flow_make_trails_rcpp(SEXP particlesSEXP, SEXP field_dfSEXP, SEXP step_lengthSEXP, SEXP max_stepsSEXP, SEXP directionSEXP, SEXP dtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type field_df(field_dfSEXP);
    Rcpp::traits::input_parameter< double >::type step_length(step_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< std::string >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    rcpp_result_gen = Rcpp::wrap(make_trails_rcpp(particles, field_df, step_length, max_steps, direction, dtest));
    return rcpp_result_gen;
END_RCPP
}
// pack_circles_cpp
DataFrame pack_circles_cpp(int width, int height, int max_circles, double r, int max_attempts, int seed);
RcppExport SEXP _flow_pack_circles_cpp(SEXP widthSEXP, SEXP heightSEXP, SEXP max_circlesSEXP, SEXP rSEXP, SEXP max_attemptsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    Rcpp::traits::input_parameter< int >::type height(heightSEXP);
    Rcpp::traits::input_parameter< int >::type max_circles(max_circlesSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type max_attempts(max_attemptsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(pack_circles_cpp(width, height, max_circles, r, max_attempts, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flow_make_trails_cpp", (DL_FUNC) &_flow_make_trails_cpp, 8},
    {"_flow_make_trail", (DL_FUNC) &_flow_make_trail, 8},
    {"_flow_make_trails_rcpp", (DL_FUNC) &_flow_make_trails_rcpp, 6},
    {"_flow_pack_circles_cpp", (DL_FUNC) &_flow_pack_circles_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_flow(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
