#' Generate flowfield trails
#'
#' @param flowfield A data frame containing three columns, x, y, and a: the
#'   coordinates of the flowfield grid and the angle associated with each grid
#'   reference.
#' @param particles A data frame containing two columns, x and y: the
#'   coordinates of the starting points (particles) from which trails will be
#'   generated on the flowfield.
#' @param max_steps The maximum number of steps to take from each trail seed
#'   point. Trails will be terminated prior to taking this many steps if the go
#'   out of bounds of the flowfield, or, if \code{dtest} is specified, when a
#'   trail becomes too close to an existing point.
#' @param step_length The length of each 'step' taken as a trail is generated.
#'   Smaller values result in smoother lines.
#' @param direction Which direction should a trail be drawn from a seedpoint?
#'   One of 'forward', 'backward' or 'both'.
#' @param dtest The test distance for terminating a trail when it becomes too
#'   close to an existing point. Larger values result in sparser trails. A value
#'   of zero means lines will not be terminated.
#'
#' @return
#' @export
#'
#' @examples make_flowfield() |> make_trails() |> draw_trails()
make_trails <- function(flowfield,
                        particles = particles_poisson(100, lims(flowfield)),
                        max_steps = 1,
                        step_length = .01,
                        direction = c("both", "forward", "backward"),
                        dtest = 0) {


  direction <- match.arg(direction)

  ff_width <- max(flowfield$x) # - min(flowfield_df$x)
  ff_height <- max(flowfield$y) # - min(flowfield_df$y)
  max_steps <- as.integer(ff_width * max_steps)
  step_length <- ff_width * step_length


  make_trails_rcpp(field_df = flowfield,
                   particles = particles,
                   max_steps = max_steps,
                   step_length = step_length,
                   direction = direction,
                   dtest = dtest)

}


particles_poisson <- function(n, limits) {
  pack_circles_cpp(limits[1], limits[2],
                   max_circles = n,
                   r = limits[1] / sqrt(n) / 2)
}

particles_grid <- function(n, limits) {
  ambient::long_grid(x = seq(1, limits[1], length.out=sqrt(n)),
                     y = seq(1, limits[2], length.out=sqrt(n)))
}

particles_unif <- function(n, limits) {
  data.frame(x = runif(n, min = 1, max = limits[1]),
             y = runif(n, min = 1, max = limits[2]))
}

lims <- function(df) {
  c(max(df$x), max(df$y))
}

