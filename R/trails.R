#' Generate flowfield trails
#'
#' @param flowfield_df
#' @param particles
#' @param start_positions
#' @param max_steps
#' @param step_length
#' @param direction
#' @param dtest
#'
#' @return
#' @export
#'
#' @examples
make_trails <- function(flowfield_df,
                        particles = 100,
                        start_positions = c("grid", "runif", "poisson"),
                        max_steps = 1,
                        step_length = .05,
                        direction = c("both", "forward", "backwards"),
                        dtest = 0) {

  start_positions <- match.arg(start_positions,
                               choices = c("grid",
                                           "runif",
                                           "poisson"))

  direction <- match.arg(direction,
                         choices = c("both",
                                     "forward",
                                     "backward"))

  ff_width <- max(flowfield_df$x) # - min(flowfield_df$x)
  ff_height <- max(flowfield_df$y) # - min(flowfield_df$y)
  max_steps <- as.integer(ff_width * max_steps)
  step_length <- ff_width * step_length

  if(start_positions == "poisson") {
    particles <- pack_circles_cpp(ff_width, ff_height,
                                  max_circles = particles,
                                  r = ff_width / sqrt(particles) / 2)
  }

  if(start_positions == "grid") {
    particles <- ambient::long_grid(x = seq(1, ff_width, length.out=sqrt(particles)),
                                    y = seq(1, ff_height, length.out=sqrt(particles)))
  }

  if(start_positions == "runif") {
    particles <- data.frame(x = runif(particles, min = 1, max = ff_width),
                            y = runif(particles, min = 1, max = ff_height))
  }

  make_trails_rcpp(field_df = flowfield_df,
                   particles = particles,
                   max_steps = max_steps,
                   step_length = step_length,
                   direction = direction,
                   dtest = dtest)

}

# make_flowfield(w = 200, h = 100, f = .015) |>
#   # dplyr::mutate(angle = round(angle, 0)) |>
#   make_trails(particles = 100,
#               start_positions = "p",
#               max_steps = .1,
#               direction = "forward",
#               step_length = .01,
#               dtest = 1) |>
#   draw_trails()


# ggplot2::ggplot() +
#   ggplot2::geom_path(data = t, ggplot2::aes(x, y, group = group), alpha = .5) +
#   ggplot2::coord_fixed() +
#   ggplot2::theme_void()

# some way of figuring max circles and radius
# seq(1, 20, length.out = 10)
# seq(1, 10, length.out = 10)

# pack_circles_cpp(40, 40, max_circles = 100, r = 40/sqrt(100)/2, seed = runif(1, 1, 10000)) |>
#   ggplot2::ggplot(ggplot2::aes(x, y)) +
#   ggplot2::geom_point() +
#   ggplot2::coord_fixed()
#
# make_flowfield(w = 150, h = 150, seed = 2, angle = 2) |>
#   make_trails(start = "p", particles = 1500, steps = .2, step_length = .01, min_dist = 1) |>
#   # dplyr::bind_rows(dplyr::mutate(draw::polygon(r = 70), group = -1) + 50) |>
#   ggplot2::ggplot(ggplot2::aes(x, y, group = group)) +
#   ggplot2::geom_path() +
#   ggplot2::coord_fixed()

# ff <- make_flowfield(w = 200, h = 200, plot = TRUE)
#
# make_trail(ff)|>
#     ggplot2::ggplot(ggplot2::aes(x, y, group = group)) +
#     ggplot2::geom_path() +
#     ggplot2::coord_fixed()
#
# ff |>
#   make_trails(particles = 100) |>
#   ggplot2::ggplot(ggplot2::aes(x, y, group = group)) +
#   ggplot2::geom_path() +
#   ggplot2::coord_fixed()
