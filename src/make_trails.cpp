#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

bool check_neighbors(double x, double y, DataFrame existing_points, double dtest);

// [[Rcpp::export]]
DataFrame make_trail(double x0,
                     double y0,
                     DataFrame field_df,
                     DataFrame existing_points,
                     double step_length = 1,
                     int max_steps = 1000,
                     std::string direction = "both",
                     double dtest = 1) {

  NumericVector ff_x = field_df["x"];
  NumericVector ff_y = field_df["y"];
  NumericVector ff_a = field_df["angle"];

  // NumericVector existing_x = existing_points["x"];
  // NumericVector existing_y = existing_points["y"];
  // int n_existing_points = existing_x.size();

  int min_x = min(ff_x);
  int max_x = max(ff_x);
  int min_y = min(ff_y);
  int max_y = max(ff_y);

  NumericVector x_out(max_steps);
  NumericVector y_out(max_steps);
  int initial_step;

  double grid_angle;
  double x_step;
  double y_step;
  double prev_x;
  double prev_y;
  double new_x;
  double new_y;

  // int step = 500;
  if(direction == "forward") {
    initial_step = 0;
  } else if(direction == "backward") {
    initial_step = max_steps - 1;
  } else if(direction == "both") {
    initial_step = floor(max_steps / 2);
    max_steps = initial_step;
  } else {
    stop("'direction' must be one of 'forward', 'backward', or 'both'.");
  }

  x_out[initial_step] = x0;
  y_out[initial_step] = y0;

  bool valid = TRUE;

  // grow forward
  if(direction == "forward" | direction == "both") {
    for(int step = (initial_step + 1); step < (initial_step + max_steps); ++step) {

      valid = TRUE;
      prev_x = x_out[step - 1];
      prev_y = y_out[step - 1];

      grid_angle = ff_a[floor(prev_x)+(floor(prev_y)-1)*max_x - 1];

      x_step = step_length * cos(grid_angle);
      y_step = step_length * sin(grid_angle);

      new_x = prev_x + x_step;
      new_y = prev_y + y_step;

      // Check if it's outside the bounds of the flowfield
      if((new_x < min_x) |
         (new_x > max_x) |
         (new_y < min_y) |
         (new_y > max_y)) {

        break;
        // If it's out of bounds, the point won't be recorded
      }

      // Now check if it's too close to any existing points
      if(dtest > 0) {
        valid = check_neighbors(new_x, new_y, existing_points, dtest);
      }

      if(valid) {
        x_out[step] = new_x;
        y_out[step] = new_y;
      } else {
        break;
      }
    }
  }


  // grow backward
  if(direction == "backward" | direction == "both") {
    for(int step = (initial_step - 1); step > (initial_step - max_steps); --step) {

      valid = TRUE;
      prev_x = x_out[step + 1];
      prev_y = y_out[step + 1];

      grid_angle = ff_a[floor(prev_x)+(floor(prev_y)-1)*max_x - 1];

      x_step = step_length * cos(grid_angle + M_PI);
      y_step = step_length * sin(grid_angle + M_PI);

      new_x = prev_x + x_step;
      new_y = prev_y + y_step;

      // check if it's outside the bounds of the flowfield
      if((new_x < min_x) |
         (new_x > max_x) |
         (new_y < min_y) |
         (new_y > max_y)) {
        break;
      }

      // now check if it's too close to any existing points
      if(dtest > 0) {
        valid = check_neighbors(new_x, new_y, existing_points, dtest);
      }

      if(valid) {
        x_out[step] = new_x;
        y_out[step] = new_y;
      } else {
        break;
      }
    }
  }

  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  // group.erase(std::remove(group.begin(), group.end(), 0), group.end());

  return(DataFrame::create(_["x"]= x_out, _["y"]= y_out));

}


// [[Rcpp::export]]
DataFrame make_trails_rcpp(DataFrame particles,
                           DataFrame field_df,
                           double step_length = 1,
                           int max_steps = 1000,
                           std::string direction = "both",
                           double dtest = 1) {

  NumericVector particles_x = particles["x"];
  NumericVector particles_y = particles["y"];
  int n_particles = particles_x.size();

  NumericVector all_points_x(n_particles * max_steps);
  NumericVector all_points_y(n_particles * max_steps);
  NumericVector all_points_g(n_particles * max_steps);

  DataFrame init_df = DataFrame::create(_["x"]= 0, _["y"]= 0);
  int n_rows = 0;

  for(int i = 0; i < n_particles; ++i) {

    // Rcout << i;

    DataFrame new_line;

    double x = particles_x[i];
    double y = particles_y[i];

    if(i == 0) {

      new_line = make_trail(x, y,
                            field_df,
                            init_df,
                            step_length, max_steps, direction, dtest);

    } else {

      new_line = make_trail(x, y,
                            field_df,
                            DataFrame::create(_["x"]= all_points_x[Rcpp::Range(0, n_rows - 1)],
                                              _["y"]= all_points_y[Rcpp::Range(0, n_rows - 1)]),
                            step_length, max_steps, direction, dtest);

    }

    NumericVector new_x = new_line["x"];
    NumericVector new_y = new_line["y"];
    int n_new_rows = new_x.size();
    // NumericVector group(n_new_rows);

    for(int j = 0; j < n_new_rows; ++j) {
      all_points_x[n_rows + j] = new_x[j];
      all_points_y[n_rows + j] = new_y[j];
      all_points_g[n_rows + j] = i;
    }

    n_rows = n_rows + n_new_rows;

  }

  return(DataFrame::create(_["x"]= all_points_x[Rcpp::Range(0, n_rows - 1)],
                           _["y"]= all_points_y[Rcpp::Range(0, n_rows - 1)],
                           _["group"]= all_points_g[Rcpp::Range(0, n_rows - 1)]));

}


bool check_neighbors(double x, double y, DataFrame existing_points, double dtest) {

  NumericVector existing_x = existing_points["x"];
  NumericVector existing_y = existing_points["y"];
  int n_existing = existing_x.size();
  // Rcout << n_existing;

  for(int i = 0; i < n_existing; ++i) {

    // Rcout << i;
    double dx = x - existing_x[i];
    double dy = y - existing_y[i];

    double dist = dx * dx + dy * dy;
    //
    if(dist < (dtest * dtest)) {
      return(FALSE);
    }
  }
  return(TRUE);
}


