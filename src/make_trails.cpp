#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame make_trails_cpp(DataFrame field_df,
                          DataFrame particles_df,
                          int width,
                          int height,
                          int steps,
                          double step_length,
                          double min_dist,
                          int seed) {

  NumericVector x = field_df["x"];
  NumericVector y = field_df["y"];
  NumericVector a = field_df["angle"];
  NumericVector px = particles_df["x"];
  NumericVector py = particles_df["y"];
  int pn = px.size();
  NumericVector x_out(pn * steps);
  NumericVector y_out(pn * steps);
  NumericVector group(pn * steps);
  // int x_start;
  // int y_start;
  double grid_angle;
  double x_step;
  double y_step;
  double prev_x;
  double prev_y;
  double new_x;
  double new_y;
  NumericVector collisions;
  int collisions_length;

  srand(seed);

  for(int p = 0; p < pn; ++p) {

    bool c = FALSE;
    // x_start = rand() % width + 1;
    // y_start = rand() % height + 1;

    for(int s = 0; s < steps; ++s) {

      if(s==0) {
        //if it's the first step of a trail, pick a random start point
        new_x = px[p];
        new_y = py[p];
      } else {

        prev_x = x_out[p * steps + (s - 1)];
        prev_y = y_out[p * steps + (s - 1)];

        if(prev_x > width | prev_x < 1) break;
        if(prev_y > height | prev_y < 1) break;

        grid_angle = a[floor(prev_x)+(floor(prev_y)-1)*width - 1];

        x_step = step_length * cos(grid_angle);
        y_step = step_length * sin(grid_angle);

        new_x = prev_x + x_step;
        new_y = prev_y + y_step;

      }

      // now check of new values collide with any existing values
      // break if an overlap is found
      if(p > 0 & min_dist > 0) {
        collisions = x_out[group < p &
          x_out > new_x - min_dist &
          x_out < new_x + min_dist &
          y_out > new_y - min_dist &
          y_out < new_y + min_dist];

        // for(int yc = 0; yc < xc.size() + 1; ++yc) {
        //
        // }

        collisions_length = collisions.size();

        if(collisions_length != 0) {
          s = steps;
          c = TRUE;
        }
      }

      if(c) break;

      x_out[p * steps + s] = new_x;
      y_out[p * steps + s] = new_y;
      group[p * steps + s] = p+1;

    }
  }

  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  group.erase(std::remove(group.begin(), group.end(), 0), group.end());

  return DataFrame::create(_["x"]= x_out, _["y"]= y_out, _["group"]= group);

}


// [[Rcpp::export]]
DataFrame make_trail(double x0,
                     double y0,
                     DataFrame field_df,
                     DataFrame existing_points,
                     double step_length,
                     double dtest) {

  NumericVector ff_x = field_df["x"];
  NumericVector ff_y = field_df["y"];
  NumericVector ff_a = field_df["angle"];

  NumericVector existing_x = existing_points["x"];
  NumericVector existing_y = existing_points["y"];
  int n_existing_points = existing_x.size();

  int min_x = min(ff_x);
  int max_x = max(ff_x);
  int min_y = min(ff_y);
  int max_y = max(ff_y);

  NumericVector x_out(1000);
  NumericVector y_out(1000);

  double grid_angle;
  double x_step;
  double y_step;
  double prev_x;
  double prev_y;
  double new_x;
  double new_y;

  // int step = 500;

  x_out[500] = x0;
  y_out[500] = y0;

  bool valid = TRUE;

  // grow forward =
  for(int step = 501; step < 1000; ++step) {

    valid = TRUE;
    prev_x = x_out[step - 1];
    prev_y = y_out[step - 1];

    grid_angle = ff_a[floor(prev_x)+(floor(prev_y)-1)*max_x - 1];

    x_step = step_length * cos(grid_angle);
    y_step = step_length * sin(grid_angle);

    new_x = prev_x + x_step;
    new_y = prev_y + y_step;

    // check if it's outside the bounds of the flowfield
    if(new_x < min_x |
       new_x > max_x |
       new_y < min_y |
       new_y > max_y) {

      // valid = FALSE;
      break;
    }

    // now check if it's too close to any existing points
    for(int i = 0; i < n_existing_points; ++i) {

      double dx = new_x - existing_x[i];
      double dy = new_y - existing_y[i];

      double dist = dx * dx + dy * dy;

      if(dist < (dtest * dtest)) {
        valid = FALSE;
        break;
      }

    }

    if(valid) {
      x_out[step] = new_x;
      y_out[step] = new_y;
    }
  }


  // grow backward
  for(int step = 499; step > 0; --step) {

    valid = TRUE;
    prev_x = x_out[step + 1];
    prev_y = y_out[step + 1];

    grid_angle = ff_a[floor(prev_x)+(floor(prev_y)-1)*max_x - 1];

    x_step = step_length * cos(grid_angle + M_PI);
    y_step = step_length * sin(grid_angle + M_PI);

    new_x = prev_x + x_step;
    new_y = prev_y + y_step;

    // check if it's outside the bounds of the flowfield
    if(new_x < min_x |
       new_x > max_x |
       new_y < min_y |
       new_y > max_y) {
      break;
    }

    // now check if it's too close to any existing points
    for(int i = 0; i < n_existing_points; ++i) {

      double dx = new_x - existing_x[i];
      double dy = new_y - existing_y[i];

      double dist = dx * dx + dy * dy;

      if(dist < (dtest * dtest)) {
        valid = FALSE;
        break;
      }

    }

    if(valid) {
      x_out[step] = new_x;
      y_out[step] = new_y;
    }
  }

  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  // group.erase(std::remove(group.begin(), group.end(), 0), group.end());

  return(DataFrame::create(_["x"]= x_out, _["y"]= y_out));

}
