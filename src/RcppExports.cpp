// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
DataFrame make_trails_rcpp(DataFrame particles, List flowfields, double step_length, int max_steps, std::string direction, double dtest, Nullable<DataFrame> existing_trails);
RcppExport SEXP _flow_make_trails_rcpp(SEXP particlesSEXP, SEXP flowfieldsSEXP, SEXP step_lengthSEXP, SEXP max_stepsSEXP, SEXP directionSEXP, SEXP dtestSEXP, SEXP existing_trailsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type particles(particlesSEXP);
    Rcpp::traits::input_parameter< List >::type flowfields(flowfieldsSEXP);
    Rcpp::traits::input_parameter< double >::type step_length(step_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< std::string >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    Rcpp::traits::input_parameter< Nullable<DataFrame> >::type existing_trails(existing_trailsSEXP);
    rcpp_result_gen = Rcpp::wrap(make_trails_rcpp(particles, flowfields, step_length, max_steps, direction, dtest, existing_trails));
    return rcpp_result_gen;
END_RCPP
}
// get_angle2
double get_angle2(double x, double y, int w);
RcppExport SEXP _flow_get_angle2(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(get_angle2(x, y, w));
    return rcpp_result_gen;
END_RCPP
}
// check_neighbors
bool check_neighbors(double x, double y, DataFrame existing_points, double dtest);
RcppExport SEXP _flow_check_neighbors(SEXP xSEXP, SEXP ySEXP, SEXP existing_pointsSEXP, SEXP dtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< DataFrame >::type existing_points(existing_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    rcpp_result_gen = Rcpp::wrap(check_neighbors(x, y, existing_points, dtest));
    return rcpp_result_gen;
END_RCPP
}
// update_current_seeds
DataFrame update_current_seeds(DataFrame line, DataFrame flowfield_df, double dsep);
RcppExport SEXP _flow_update_current_seeds(SEXP lineSEXP, SEXP flowfield_dfSEXP, SEXP dsepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type line(lineSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type flowfield_df(flowfield_dfSEXP);
    Rcpp::traits::input_parameter< double >::type dsep(dsepSEXP);
    rcpp_result_gen = Rcpp::wrap(update_current_seeds(line, flowfield_df, dsep));
    return rcpp_result_gen;
END_RCPP
}
// path_to_polygon
DataFrame path_to_polygon(DataFrame path, double thickness);
RcppExport SEXP _flow_path_to_polygon(SEXP pathSEXP, SEXP thicknessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type path(pathSEXP);
    Rcpp::traits::input_parameter< double >::type thickness(thicknessSEXP);
    rcpp_result_gen = Rcpp::wrap(path_to_polygon(path, thickness));
    return rcpp_result_gen;
END_RCPP
}
// check_valid
bool check_valid(double x, double y, int w, int h, DataFrame existing_points, double dtest);
RcppExport SEXP _flow_check_valid(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP hSEXP, SEXP existing_pointsSEXP, SEXP dtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type existing_points(existing_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    rcpp_result_gen = Rcpp::wrap(check_valid(x, y, w, h, existing_points, dtest));
    return rcpp_result_gen;
END_RCPP
}
// make_trail_oop
DataFrame make_trail_oop(double x0, double y0, int direction, double step_length, int max_steps, DataFrame flowfield, DataFrame existing_points, double dtest);
RcppExport SEXP _flow_make_trail_oop(SEXP x0SEXP, SEXP y0SEXP, SEXP directionSEXP, SEXP step_lengthSEXP, SEXP max_stepsSEXP, SEXP flowfieldSEXP, SEXP existing_pointsSEXP, SEXP dtestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< int >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< double >::type step_length(step_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type flowfield(flowfieldSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type existing_points(existing_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type dtest(dtestSEXP);
    rcpp_result_gen = Rcpp::wrap(make_trail_oop(x0, y0, direction, step_length, max_steps, flowfield, existing_points, dtest));
    return rcpp_result_gen;
END_RCPP
}
// test_class_Flowfield
void test_class_Flowfield(DataFrame df, int x, int y);
RcppExport SEXP _flow_test_class_Flowfield(SEXP dfSEXP, SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    test_class_Flowfield(df, x, y);
    return R_NilValue;
END_RCPP
}
// test_class_pointOnFlowfield
void test_class_pointOnFlowfield(double x, double y, DataFrame ff_df, int step_length);
RcppExport SEXP _flow_test_class_pointOnFlowfield(SEXP xSEXP, SEXP ySEXP, SEXP ff_dfSEXP, SEXP step_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< DataFrame >::type ff_df(ff_dfSEXP);
    Rcpp::traits::input_parameter< int >::type step_length(step_lengthSEXP);
    test_class_pointOnFlowfield(x, y, ff_df, step_length);
    return R_NilValue;
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
    {"_flow_make_trail", (DL_FUNC) &_flow_make_trail, 8},
    {"_flow_make_trails_rcpp", (DL_FUNC) &_flow_make_trails_rcpp, 7},
    {"_flow_get_angle2", (DL_FUNC) &_flow_get_angle2, 3},
    {"_flow_check_neighbors", (DL_FUNC) &_flow_check_neighbors, 4},
    {"_flow_update_current_seeds", (DL_FUNC) &_flow_update_current_seeds, 3},
    {"_flow_path_to_polygon", (DL_FUNC) &_flow_path_to_polygon, 2},
    {"_flow_check_valid", (DL_FUNC) &_flow_check_valid, 6},
    {"_flow_make_trail_oop", (DL_FUNC) &_flow_make_trail_oop, 8},
    {"_flow_test_class_Flowfield", (DL_FUNC) &_flow_test_class_Flowfield, 3},
    {"_flow_test_class_pointOnFlowfield", (DL_FUNC) &_flow_test_class_pointOnFlowfield, 4},
    {"_flow_pack_circles_cpp", (DL_FUNC) &_flow_pack_circles_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_flow(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
