#' Create a flowfied
#'
#' @param width
#' @param height
#' @param noise
#' @param angle
#' @param octaves
#' @param frequency
#' @param fractal
#' @param seed
#' @param plot
#'
#' @return A data.frame with three columns: x and y coordinates, and angle, the
#'   angle of flow at that grid coordinate
#' @export
#'
#' @examples
make_flowfield <- function(width = 100, height = 100,
                           noise = c("perlin", "simplex", "cubic", "value", "waves", "spheres", "white", "checkerboard"),
                           angle = 1,
                           octaves = 1,
                           frequency = .01,
                           fractal = c("fbm", "billow", "clamped", "ridged"),
                           seed = runif(1, max = 10000),
                           plot = FALSE) {

  fractal <- match.arg(fractal)
  # fractal <- switch(fractal,
  #                   "none" = ambient::none,
  #                   "fbm" = ambient::fbm,
  #                   "billow" = ambient::billow,
  #                   "clamped" = ambient::clamped,
  #                   "ridged" = ambient::ridged)

  noise <- match.arg(noise)
  noise <- switch(noise,
                  "perlin"  = ambient::gen_perlin,
                  "simplex" = ambient::gen_simplex,
                  "cubic"   = ambient::gen_cubic,
                  "value"   = ambient::gen_value,
                  "waves"   = ambient::gen_waves,
                  "spheres" = ambient::gen_spheres,
                  "white"   = ambient::gen_white,
                  "checkerboard" = ambient::gen_checkerboard)

  ff <- ambient::long_grid(x = 1:height, y = 1:width) |>
    dplyr::rename(x = y, y = x) |>
    dplyr::mutate(n = ambient::fracture(noise = noise,
                                        x = x*frequency,
                                        y = y*frequency,
                                        octaves = octaves,
                                        fractal = fractal,
                                        seed = seed),
                  angle = n*angle*pi)

  if(plot) {
    pic <- draw_field(ff)
    plot(pic)
  }

  return(ff)

}


draw_field <- function(df) {
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = df,
                          ggplot2::aes(x = x, xend = x + cos(angle)*.9, y = y, yend = y + sin(angle)*.9),
                          # arrow = ggplot2::arrow(angle = df$angle, length = ggplot2::unit(.1, "inches"))
                          ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}

draw_trails <- function(df, size = 1, alpha = 1) {
  ggplot2::ggplot(data = df, ggplot2::aes(x, y, group = group)) +
    ggplot2::geom_path(alpha = alpha, size = size) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    NULL
}
