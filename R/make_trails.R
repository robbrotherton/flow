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
                        start_positions,
                        particles = 100,
                        max_steps = 1,
                        step_length = .01,
                        direction = c("both", "forward", "backward"),
                        dtest = 0) {

  # start_positions <- match.arg(start_positions,
  #                              choices = c("grid",
  #                                          "runif",
  #                                          "poisson"))

  direction <- match.arg(direction)

  ff_width <- max(flowfield_df$x) # - min(flowfield_df$x)
  ff_height <- max(flowfield_df$y) # - min(flowfield_df$y)
  max_steps <- as.integer(ff_width * max_steps)
  step_length <- ff_width * step_length

  if (is.data.frame(start_positions)) {
    particles <- start_positions
  } else if (start_positions == "poisson")  {
    particles <- pack_circles_cpp(ff_width, ff_height,
                                  max_circles = particles,
                                  r = ff_width / sqrt(particles) / 2)

  } else if (start_positions == "grid") {
    particles <- ambient::long_grid(x = seq(1, ff_width, length.out=sqrt(particles)),
                                    y = seq(1, ff_height, length.out=sqrt(particles)))
  } else if (start_positions == "runif") {
    particles <- data.frame(x = runif(particles, min = 1, max = ff_width),
                            y = runif(particles, min = 1, max = ff_height))
  } else {
    stop("start positions not specified correctly")
  }

  make_trails_rcpp(field_df = flowfield_df,
                   particles = particles,
                   max_steps = max_steps,
                   step_length = step_length,
                   direction = direction,
                   dtest = dtest)

}
