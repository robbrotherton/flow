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
