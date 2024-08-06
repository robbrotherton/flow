#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

bool check_valid(double x, double y, int w, int h, DataFrame existing_points, double dtest);
bool check_neighbors(double x, double y, DataFrame existing_points, double dtest);
double get_angle(double x, double y, DataFrame flowfield_df);
double get_angle2(double x, double y, int w);

// [[Rcpp::export]]
DataFrame make_long_trail(double x0,
                          double y0,
                          DataFrame flowfield,
                          DataFrame existing_points,
                          double step_length,
                          double dtest,
                          int max_steps = 1000) {

  NumericVector ff_x = flowfield["x"];
  NumericVector ff_y = flowfield["y"];
  NumericVector ff_a = flowfield["angle"];

  int w = max(ff_x);
  int h = max(ff_y);

  int step;

  double x;
  double y;
  double a;
  double a0 = ff_a[get_angle2(x0, y0, w)];

  NumericVector x_back(1000);
  NumericVector y_back(1000);
  NumericVector a_back(1000);

  NumericVector x_out(1000);
  NumericVector y_out(1000);
  NumericVector a_out(1000);

  bool valid;


  // ------------------------------------- GROW LINE BACKWARD
  x = x0;
  y = y0;
  a = a0;
  step = 0;
  x_back[step] = x;
  y_back[step] = y;
  a_back[step] = a;

  valid = TRUE;

  while(valid) {

    x = x - step_length * cos(a);
    y = y - step_length * sin(a);

    valid = check_valid(x, y, w, h, existing_points, dtest);

    if(valid) {
      step++;
      a = ff_a[get_angle2(x, y, w)];

      x_back[step] = x;
      y_back[step] = y;
      a_back[step] = a;

    }
  }

  // Now reverse the points so far

  x_back = rev(x_back[Range(0, step)]);
  y_back = rev(y_back[Range(0, step)]);
  a_back = rev(a_back[Range(0, step)]);

  for(int i = 0; i < step + 1; i++) {
    x_out[i] = x_back[i];
    y_out[i] = y_back[i];
    a_out[i] = a_back[i];
  }

  // for(int i = 0; i < step + 1; i++) {
  //   x_out[i] = x_back[step - i];
  //   y_out[i] = y_back[step - i];
  //   a_out[i] = a_back[step - i];
  // }


  // GROW LINE FORWARD --------------------------------------

  x = x0;
  y = y0;
  a = a0;
  valid = TRUE;

  while(valid) {

    x = x + step_length * cos(a);
    y = y + step_length * sin(a);

    valid = check_valid(x, y, w, h, existing_points, dtest);

    if(valid) {
      step++;
      a = ff_a[get_angle2(x, y, w)];

      x_out[step] = x;
      y_out[step] = y;
      a_out[step] = a;
    }
  }

  return(DataFrame::create(_["x"]= x_out[Range(0, step)],
                           _["y"]= y_out[Range(0, step)],
                           _["a"]= a_out[Range(0, step)]));

}


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
  NumericVector a_out(max_steps);
  int initial_step;

  double grid_angle;
  double prev_a;
  double x_step;
  double y_step;
  double prev_x;
  double prev_y;
  double new_x;
  double new_y;

  // need to keep track of what indices are actually used, so that only these
  // can be returned
  int start_index;
  int end_index;

  // int step = 500;
  if(direction == "forward") {
    initial_step = 0;
    start_index = initial_step;
  } else if(direction == "backward") {
    initial_step = max_steps - 1;
    end_index = initial_step;
  } else if(direction == "both") {
    initial_step = floor(max_steps / 2);
    max_steps = initial_step;
    start_index = initial_step;
    end_index = initial_step;
  } else {
    stop("'direction' must be one of 'forward', 'backward', or 'both'.");
  }

  grid_angle = get_angle(x0, y0, field_df);

  x_out[initial_step] = x0;
  y_out[initial_step] = y0;
  a_out[initial_step] = grid_angle;

  bool valid = TRUE;

  // grow forward
  if(direction == "forward" | direction == "both") {
    for(int step = (initial_step + 1); step < (initial_step + max_steps); ++step) {

      valid = TRUE;
      prev_x = x_out[step - 1];
      prev_y = y_out[step - 1];
      prev_a = a_out[step - 1];

      // grid_angle = get_angle(prev_x, prev_y, field_df);

      x_step = step_length * cos(prev_a);
      y_step = step_length * sin(prev_a);

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
        a_out[step] = get_angle(new_x, new_y, field_df);

        end_index = step;
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
      prev_a = a_out[step + 1];

      // grid_angle = get_angle(prev_x, prev_y, field_df);
      // grid_angle = ff_a[floor(prev_x)+(floor(prev_y)-1)*max_x - 1];

      x_step = step_length * cos(prev_a + M_PI);
      y_step = step_length * sin(prev_a + M_PI);

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
        a_out[step] = get_angle(new_x, new_y, field_df);

        start_index = step;
      } else {
        break;
      }
    }
  }

  // x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  // y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  // a_out.erase(std::remove(a_out.begin(), a_out.end(), 0), a_out.end());
  // group.erase(std::remove(group.begin(), group.end(), 0), group.end());

  // return(DataFrame::create(_["x"]= x_out,
  //                          _["y"]= y_out,
  //                          _["a"]= a_out));

  return(DataFrame::create(_["x"]= x_out[Rcpp::Range(start_index, end_index)],
                           _["y"]= y_out[Rcpp::Range(start_index, end_index)],
                           _["a"]= a_out[Rcpp::Range(start_index, end_index)]));

}


// [[Rcpp::export]]
DataFrame make_trails_rcpp(DataFrame particles,
                           DataFrame field_df,
                           double step_length = 1,
                           int max_steps = 1000,
                           std::string direction = "both",
                           double dtest = 1,
                           Nullable<DataFrame> existing_trails = R_NilValue) {

  NumericVector particles_x = particles["x"];
  NumericVector particles_y = particles["y"];
  int n_particles = particles_x.size();

  NumericVector all_points_x(n_particles * max_steps);
  NumericVector all_points_y(n_particles * max_steps);
  NumericVector all_points_g(n_particles * max_steps);

  // DataFrame init_df = DataFrame::create(_["x"]= 0, _["y"]= 0);

  // use existing_trails if provided, otherwise initialize empty dataframe
  NumericVector existing_x, existing_y;
  bool has_existing_trails = !existing_trails.isNull();
  Rcout << "existing?" << has_existing_trails;
  if (has_existing_trails) {
    DataFrame existing_df = as<DataFrame>(existing_trails);
    existing_x = existing_df["x"];
    existing_y = existing_df["y"];
  }

  int n_rows = 0;

  for(int i = 0; i < n_particles; ++i) {

    // Rcout << i;

    DataFrame new_line;

    double x = particles_x[i];
    double y = particles_y[i];

    if (has_existing_trails) {
      // combine existing and newly generated points
      NumericVector combined_x(existing_x.size() + n_rows);
      NumericVector combined_y(existing_y.size() + n_rows);

      // copy existing trails
      std::copy(existing_x.begin(), existing_x.end(), combined_x.begin());
      std::copy(existing_y.begin(), existing_y.end(), combined_y.begin());

      // copy newly generated points
      std::copy(all_points_x.begin(), all_points_x.begin() + n_rows, combined_x.begin() + existing_x.size());
      std::copy(all_points_y.begin(), all_points_y.begin() + n_rows, combined_y.begin() + existing_y.size());

      new_line = make_trail(x, y,
                            field_df,
                            DataFrame::create(_["x"] = combined_x,
                                              _["y"] = combined_y),
                                              step_length, max_steps, direction, dtest);

    } else {
      if(i == 0) {

        new_line = make_trail(x, y,
                              field_df,
                              DataFrame::create(_["x"]= 0, _["y"]= 0), // formerly init_df
                              step_length, max_steps, direction, dtest);

      } else {

        new_line = make_trail(x, y,
                              field_df,
                              DataFrame::create(_["x"]= all_points_x[Rcpp::Range(0, n_rows - 1)],
                                                _["y"]= all_points_y[Rcpp::Range(0, n_rows - 1)]),
                                                step_length, max_steps, direction, dtest);

      }
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

double get_angle(double x, double y, DataFrame flowfield_df) {

  NumericVector angles = flowfield_df["angle"];
  NumericVector x_col = flowfield_df["x"];
  int max_x = max(x_col);

  return angles[floor(x)+(floor(y)-1) * max_x - 1];

}

// [[Rcpp::export]]
double get_angle2(double x, double y, int w) {

  return floor(x)+(floor(y)-1) * w - 1;

}

// [[Rcpp::export]]
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


// [[Rcpp::export]]
DataFrame update_current_seeds(DataFrame line, DataFrame flowfield_df, double dsep) {
  // take the first line in the queue (the oldest line) as the new current line

  NumericVector x = line["x"];
  NumericVector y = line["y"];
  NumericVector a = line["a"];

  int n = x.length();
  NumericVector x_out(n * 2 - 1);
  NumericVector y_out(n * 2 - 1);

  // now compute all possible seed points; for each row of the df
  for(int i = 0; i < n; i++) {

    x_out[i]   = x[i] + dsep * cos(a[i] + M_PI * 0.5);
    x_out[i+n] = x[i] + dsep * cos(a[i] - M_PI * 0.5);
    y_out[i]   = y[i] + dsep * sin(a[i] + M_PI * 0.5);
    y_out[i+n] = y[i] + dsep * sin(a[i] - M_PI * 0.5);

    // check if seed is inbounds?

  }

  // randomize order before returning?

  // all_possible_seeds <- data.frame(x = x, y = y) |>
  //   dplyr::filter(inbounds(x, y, w, h)) |>
  //   dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
  //   dplyr::arrange(rand)

  return(DataFrame::create(_["x"] = x_out,
                           _["y"] = y_out));

}


// [[Rcpp::export]]
DataFrame path_to_polygon(DataFrame path, double thickness = 1) {
  // take the first line in the queue (the oldest line) as the new current line

  NumericVector x = path["x"];
  NumericVector y = path["y"];
  NumericVector a = path["a"];

  int n = x.length();
  NumericVector x_out(n * 2);
  NumericVector y_out(n * 2);

  for(int i = 0; i < n; i++) {

    x_out[i]   = x[i] + thickness * cos(a[i] + M_PI * 0.5);
    y_out[i]   = y[i] + thickness * sin(a[i] + M_PI * 0.5);

    x_out[i + n] = x[n - i - 1] + thickness * cos(a[n - i - 1] - M_PI * 0.5);
    y_out[i + n] = y[n - i - 1] + thickness * sin(a[n - i - 1] - M_PI * 0.5);

  }

  return(DataFrame::create(_["x"] = x_out,
                           _["y"] = y_out));

}



// [[Rcpp::export]]
bool check_valid(double x, double y, int w, int h, DataFrame existing_points, double dtest) {

  // Check if point is out of bounds
  if((x < 1) |
     (x > w) |
     (y < 1) |
     (y > h)) {
    return false;
  }

  // Now check if it's too close to any existing points
  if(dtest > 0) {
    return check_neighbors(x, y, existing_points, dtest);
  } else {
    return true;
  }

}
