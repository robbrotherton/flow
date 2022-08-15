#' Generate flowfield trails
#'
#' @param flowfield
#' @param particles
#' @param max_steps
#' @param step_length
#' @param direction
#' @param dtest
#'
#' @return
#' @export
#'
#' @examples
make_trails <- function(flowfield,
                        particles = particles_poisson(10, lims(flowfield)),
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

