#' Generate flowfield trails
#'
#' @param flowfield A data frame, or a list of data frames. Each flowfield data
#'   frame must contain three columns, x, y, and a: the coordinates of the
#'   flowfield grid and the angle associated with each grid reference. If a list
#'   of data frames is supplied, trails will be drawn in alternating order.
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
#' @param existing_trails_df A data frame containing previously-drawn trails.
#'   New trails will avoid collisions with these according to dtest.
#'
#' @return A data frame.
#' @export
#'
#' @examples make_flowfield() |> make_trails() |> draw_trails()
make_trails <- function(flowfield,
                        particles = trail_seeds(100, ff_limits(flowfield)),
                        step_length = 1,
                        direction = c("both", "forward", "backward"),
                        existing_trails_df = NULL) {


  if ("data.frame" %in% class(flowfield)) {
    flowfield <- list(flowfield)
  }

  direction <- match.arg(direction)

  # ff_width <- max(flowfield$x) # - min(flowfield_df$x)
  # ff_height <- max(flowfield$y) # - min(flowfield_df$y)
  # max_steps <- as.integer(ff_width * max_steps)
  # step_length <- ff_width * step_length


  make_trails_rcpp(flowfields = flowfield,
                   particles = particles,
                   step_length = step_length,
                   direction = direction,
                   existing_trails = existing_trails_df)

}

trail_seeds <- function(n,
                        limits,
                        distribution = c("poisson", "grid", "uniform"),
                        size = 0,
                        max_length = 100) {

  distribution <- match.arg(distribution, choices = c("poisson", "grid", "uniform"))
  distribution <- switch(distribution,
                         "poisson"  = particles_poisson(n, limits),
                         "grid" = particles_grid(n, limits),
                         "uniform"   = particles_unif(n, limits))

  distribution |>
    dplyr::mutate(size = rep(size, length.out = dplyr::n()),
                  max_length = rep(max_length, length.out = dplyr::n()))
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

ff_limits <- function(df) {
  if (!is.data.frame(df)) {
    df <- df[[1]]
  }
  c(max(df$x), max(df$y))
}
