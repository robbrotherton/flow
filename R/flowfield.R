#' Title
#'
#' @param w
#' @param h
#' @param angle
#' @param octaves
#' @param frequency
#' @param seed
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
make_flowfield <- function(w = 100, h = 100,
                           angle = 1,
                           octaves = 1,
                           frequency = .01,
                           seed = runif(1, max = 10000),
                           plot = FALSE) {

  ff <- ambient::long_grid(x = 1:w, y = 1:h) |>
    dplyr::rename(x = y, y = x) |>
    dplyr::mutate(n = ambient::fracture(ambient::gen_perlin,
                                        x = x*frequency,
                                        y = y*frequency,
                                        octaves = octaves,
                                        fractal = ambient::fbm,
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
                          ggplot2::aes(x = x, xend = x + cos(angle), y = y, yend = y + sin(angle))) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}

draw_trails <- function(df) {
  ggplot2::ggplot(data = df, ggplot2::aes(x, y, group = group)) +
    ggplot2::geom_path() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}

# make_flowfield(f = .01, angle = 1) |>
#   make_trails(start_positions = "p", particles = 500, step_length = .02, steps = 1) |>
#   draw_trails()
