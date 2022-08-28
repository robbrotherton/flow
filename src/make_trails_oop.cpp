#include <Rcpp.h>
// #include <algorithm>
using namespace Rcpp;
using namespace std;

class Flowfield {
  DataFrame ff;
  NumericVector x;
  NumericVector y;
  NumericVector a;
  int w;
  int h;

public:
  Flowfield() {}
  Flowfield(DataFrame flowfield_df) {
    ff = flowfield_df;
    x = ff["x"];
    y = ff["y"];
    a = ff["angle"];
    w = max(x);
    h = max(y);
  }

  int getWidth() {
    return w;
  }

  int getHeight() {
    return h;
  }

  double getAngle(double x, double y) {
    int index = floor(x)+(floor(y)-1) * w - 1;
    return a[index];
  }

};




class PointOnFlowfield {

public:
  double x;
  double y;
  double a;
  Flowfield ff;
  int step = 0;
  int w;
  int h;

  PointOnFlowfield(double x0, double y0, Flowfield ff_in) {
    ff = ff_in;
    x = x0;
    y = y0;
    a = ff.getAngle(x, y);

    w = ff.getWidth();
    h = ff.getHeight();
  }

  void takeStep(double step_length) {
    x = x + step_length * cos(a);
    y = y + step_length * sin(a);
    a = ff.getAngle(x, y);
  }

  bool isValid(DataFrame existing_points, double dtest) {
    if (pointIsOutOfBounds()) return false;

    return isSufficientDistanceFromNeighbors(existing_points, dtest);
  }

  bool pointIsOutOfBounds() {
    return
    (x < 1) |
      (x > w) |
      (y < 1) |
      (y > h);
  }

  bool isSufficientDistanceFromNeighbors(DataFrame existing_points, double dtest) {
    // code goes here
    return true;
  }
};



class Trail{
  vector<double> x_out;
  vector<double> y_out;
  vector<double> a_out;
  int step;
  int max_steps;

public:
  Trail(int _max_steps) {
    x_out.resize(_max_steps);
    y_out.resize(_max_steps);
    a_out.resize(_max_steps);
    max_steps = _max_steps;
    step = 0;
  }


  void addPointToTrail(PointOnFlowfield p) {
    x_out[step] = p.x;
    y_out[step] = p.y;
    a_out[step] = p.a;
    Rcout << step;
    step++;
  }

  bool maxStepsNotReached() {
    return step < max_steps;
  }

  DataFrame makeDataFrame() {
    x_out.resize(step);
    y_out.resize(step);
    a_out.resize(step);

    return DataFrame::create(_["x"]= x_out,
                             _["y"]= y_out,
                             _["a"]= a_out);
  }

};

// [[Rcpp::export]]
DataFrame make_trail_oop(double x0, double y0,
                         int direction,
                         double step_length,
                         int max_steps,
                         DataFrame flowfield,
                         DataFrame existing_points,
                         double dtest) {

  Flowfield ff (flowfield);
  PointOnFlowfield point (x0, y0, ff);
  Trail trail (max_steps);

  while (point.isValid(existing_points, dtest) & trail.maxStepsNotReached()) {

    point.takeStep(step_length);

    if (point.isValid(existing_points, dtest)) {
      trail.addPointToTrail(point);
    }
  }

  return trail.makeDataFrame();

}



// [[Rcpp::export]]
void test_class_Flowfield(DataFrame df, int x, int y) {
  Flowfield ff (df);
  Rcout << "max x value is: " << ff.getWidth() << endl;
  Rcout << "max y value is: " << ff.getHeight() << endl;
  Rcout << "angle: " << ff.getAngle(x, y) << endl;
}


// [[Rcpp::export]]
void test_class_pointOnFlowfield(double x, double y, DataFrame ff_df, int step_length) {
  Flowfield ff (ff_df);
  PointOnFlowfield p (x, y, ff);
  DataFrame existing_points = DataFrame::create();
  Rcout << "current x is: " << p.x << endl;
  Rcout << "current y is: " << p.y << endl;
  Rcout << "current a is: " << p.a << endl;
  Rcout << "point is valid: " << p.isValid(existing_points, 1) << endl;
  p.takeStep(step_length);
  Rcout << "next point is: " << p.x << "," << p.y << "," << p.a << ";" << p.isValid(existing_points, 1) <<endl;
}
