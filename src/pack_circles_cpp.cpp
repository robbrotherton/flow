#include <Rcpp.h>
#include <math.h>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame pack_circles_cpp(int width,
                           int height,
                           int max_circles,
                           double r,
                           int max_attempts = 2000,
                           int seed = 1) {

  int n;
  int attempt = 0;
  bool overlap;
  double new_x;
  double new_y;
  double dist;
  NumericVector x_out(max_circles+1);
  NumericVector y_out(max_circles+1);
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> randX(1,width);
  std::uniform_real_distribution<double> randY(1,height);

  // start with an initial point
  x_out[0] = randX(generator);
  y_out[0] = randY(generator);
  // x_out[0] = (rand() % (width*100) + 1);
  // y_out[0] = (rand() % (height*100) + 1);
  n = 1;

  while(attempt <= max_attempts) {

    ++attempt;
    // Rcout << "circle" << n << "//attempt" << attempt << "\n";
    overlap = FALSE;

    new_x = randX(generator);
    new_y = randY(generator);

    // now check the values against existing circles
    // returns TRUE if circles overlap
    for(int i = 0; i < n; ++i) {

      dist = sqrt(pow(x_out[i]-new_x, 2) + pow(y_out[i]-new_y, 2));

      overlap = dist < (2*r);
      // if there's overlap, stop checking
      if(overlap) break;

    }

    // if there's no overlap after checking every row...
    if(!overlap) {

      x_out[n] = new_x;
      y_out[n] = new_y;

      ++n;
      attempt = 0;

    }

    if(n >= max_circles) break; // all the circles fit

  }

  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());

  return DataFrame::create(_["x"]= x_out, _["y"]= y_out);

}
