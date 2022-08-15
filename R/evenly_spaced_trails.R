# https://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf

#' Title
#'
#' @param flowfield_df A data.frame containing rows with x and y coordinates and
#'   a corresponding angle.
#' @param dsep The separating distance of newly-generated lines, as a proportion
#'   of the width of the flowfield.
#' @param dtest The minimal distance that lines are allowed to converge. A
#'   proportion of \code{dsep}.
#' @param max_steps The maximum number of steps (i.e. increments of \code{dsep})
#'   that a line is allowed to grow. (Of course, a line will end earlier if it
#'   hits the boundary of the flowfield or gets too close to another line
#'   first.)
#'
#' @return
#' @export
#'
#' @examples
make_trails_even <- function(flowfield_df, max_steps = 100, dsep = .05, dtest = .5) {

  # need some properties of the flowfield and parameters as variables
  width <- max(flowfield_df$x)
  height <- max(flowfield_df$y)
  dsep <- width * dsep
  dtest <- dsep * dtest
  step_length <- 1

  # start by making the grid (that'll make looking for collisions more efficient)
  grid <- list(
    xvals = vector(mode = "list", length = width*height*(1/dsep^2)),
    yvals = vector(mode = "list", length = width*height*(1/dsep^2))
  )

  # for the first line, just pick a random starting point
  line_number <- 1
  starting_x <- runif(1, min = 1, max = width-1)
  starting_y <- runif(1, min = 1, max = height-1)

  # Now generate a trail from that starting point, working backwards and forwards
  new_line <- make_long_trail(starting_x,
                              starting_y,
                              flowfield_df,
                              existing_points = data.frame(x = -100, y = -100),
                              step_length = step_length,
                              dtest = dtest)

  grid <- add_to_grid(grid, new_line$x, new_line$y, width, dsep)

  # return(grid)
  # Initialize a list to collect each line as it it made. This will become the
  # output
  lines <- list(new_line)

  queue_position <- 1


  while(queue_position <= line_number) {

    # print(paste("queue position: ", queue_position))
    seeds <- get_seeds(lines[[queue_position]], dsep)

    for(i in 1:length(seeds$x)) {

      starting_x <- seeds$x[i]
      starting_y <- seeds$y[i]

      # check if the point is valid

      # valid <- check_neighbors(starting_x, starting_y,
      #                          existing_points = dplyr::bind_rows(lines),
      #                          dtest = dtest) &
      #          inbounds(starting_x, starting_y, width, height)

      valid <- check_valid(starting_x, starting_y,
                           width, height,
                           existing_points = dplyr::bind_rows(lines),
                           dtest = dtest)

      # valid <- !query_grid(grid, starting_x, starting_y, width, height, dsep, dtest)

      if(valid) {

        line_number <- line_number + 1
        # print(line_number)

        new_line <- make_trail(starting_x,
                               starting_y,
                               flowfield_df,
                               existing_points = dplyr::bind_rows(lines),
                               step_length = step_length,
                               dtest = dtest)

        # grid <- add_to_grid(grid, new_line$x, new_line$y, width, dsep)
        lines[[line_number]] <- new_line

      }
    }

    # once we tried all seed points for the current queue line, increment the
    # queue position
    queue_position <- queue_position + 1

  }

  dplyr::bind_rows(lines, .id = "group") |>
    dplyr::mutate(group = as.numeric(group))

}



inbounds <- function(x, y, w, h) {
  x > 1 & x < w & y > 1 & y < h
}

# returns TRUE if the nearest neighbor is too close
# check_neighbor <- function(data, qtree, x, y, test_distance) {
#   index <- SearchTrees::knnLookup(qtree, x, y, k = 1)
#   dist <- sqrt((x-data$x[index])^2+(y-data$y[index])^2)
#   dist < test_distance
# }

get_seeds <- function(line, dsep) {
  # take the first line in the queue (the oldest line) as the new current line

  n <- nrow(line)
  x <- numeric(n*2)
  y <- numeric(n*2)
  # now compute all possible seed points; for each row of the df
  for(i in 1:n) {

    x[i]   <- line$x[i] + dsep * cos(line$a[i] + pi/2)
    x[i+n] <- line$x[i] + dsep * cos(line$a[i] - pi/2)
    y[i]   <- line$y[i] + dsep * sin(line$a[i] + pi/2)
    y[i+n] <- line$y[i] + dsep * sin(line$a[i] - pi/2)
  }

  # all_possible_seeds <- data.frame(x = x, y = y) |>
  #   dplyr::filter(inbounds(x, y, w, h)) |>
  #   dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
  #   dplyr::arrange(rand)
  #
  # return(all_possible_seeds)

  data.frame(x = x, y = y) |>
    dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
    dplyr::arrange(rand)
}

# create new seedpoints
# update_current_line <- function(queue, dsep, w, h) {
#   # take the first line in the queue (the oldest line) as the new current line
#   line <- queue[[1]]
#   n <- nrow(line)
#   x <- numeric()
#   y <- numeric()
#   # now compute all possible seed points; for each row of the df
#   for(i in 1:n) {
#
#     x[i]   <- line$x[i] + dsep * cos(line$a[i] + pi/2)
#     x[i+n] <- line$x[i] + dsep * cos(line$a[i] - pi/2)
#     y[i]   <- line$y[i] + dsep * sin(line$a[i] + pi/2)
#     y[i+n] <- line$y[i] + dsep * sin(line$a[i] - pi/2)
#   }
#
#   all_possible_seeds <- data.frame(x = x, y = y) |>
#     dplyr::filter(inbounds(x, y, w, h)) |>
#     dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
#     dplyr::arrange(rand)
#
#   return(all_possible_seeds)
# }


# w <- 50
# h <- 50
# dsep <- w*.05 #w*.1 # can be .1, .05
# dtest <- dsep*.9
#
# df <- make_flowfield(w, h)
# draw_field(df)
#
# dft <- make_trails(df, start_positions = "grid")
# dfx <- make_trail_x(df)
#
# draw_trails(dfx)
#
#
# ffdf <- make_flowfield(50, 50, angle = 1, octaves = 1, frequency = .01, plot = FALSE)
# #
# t <- make_trail(ffdf) |>
#   tibble::rownames_to_column("point")
#
# ggplot2::ggplot(t, ggplot2::aes(x, y, group = line)) +
#   ggplot2::geom_path() +
#   ggplot2::coord_fixed()

# pic <- t |>
#   ggplot2::ggplot() +
#   ggforce::geom_circle(mapping = ggplot2::aes(x0 = x, y0 = y, r = dtest/2)) +
#   # geom_point(mapping = aes(x, y)) +
#   # geom_text(data = t, mapping = aes(x, y, label = group), size = 8, color = 'red') +
#   ggplot2::geom_path(mapping = ggplot2::aes(x, y, group = line, color = ifelse(line==2, 'red', 'black'))) +
#   ggplot2::coord_fixed() +
#   # scale_x_continuous(breaks = 1:10) +
#   # scale_y_continuous(breaks = 1:10) +
#   ggplot2::scale_color_identity() +
#   ggplot2::theme_void() +
#   NULL
# pic


grid <- function(x, y, width, dsep) {
  x_ind <- ceiling(x*1/dsep)
  y_ind <- ceiling(y*1/dsep)

  grid_ind <- x_ind + (y_ind-1) * width*1/dsep

}


get_index <- function(x, dsep) {
  ceiling(x*1/dsep)
}

add_to_grid <- function(grid, x, y, width, dsep) {

  for(i in 1:length(x)) {
    x0 <- get_index(x[i], dsep) #round(x, 0)
    y0 <- get_index(y[i], dsep) #round(y,0)
    ind <- x0 + (y0-1) * width*1/dsep
    l <- length(grid$xvals[[ind]])

    grid$xvals[[ind]][l+1] <- x[i]
    grid$yvals[[ind]][l+1] <- y[i]

  }

  return(grid)

}

grid_indices <- function(x, y, width, height, dsep) {
  cols <- width*1/dsep
  rows <- height*1/dsep
  # edges_l <- 1 + 0:(rows-1) * rows
  # edges_r <- rows + 0:(rows-1) * rows
  # edges_b <- 1:cols + 0 * cols
  # edges_t <- 1:cols + (cols-1)*cols

  xout <- ifelse(x == 1, -1,
                 ifelse(x == cols, 1, NA))

  yout <- ifelse(y == 1, -1,
                 ifelse(y == rows, 1, NA))

  # return(c(xout, yout))

  z <- tibble::tibble(xn = c(0,-1,0,1,1, 1, 0,-1,-1),
              yn = c(0, 1,1,1,0,-1,-1,-1, 0),
              xi = xn + x,
              yi = yn + y,
              grid_index = xi + (yi-1) * width*1/dsep)

  z <- dplyr::filter(z, !xn %in% xout & !yn %in% yout)

  z$grid_index
  # return(out)

}

# grid_indices(x = 2, y = 19)
# edges <- grid_edges(w, h, dsep)
#
# z <- tibble(x = c(0,-1,0,1,1, 1, 0,-1,-1),
#             y = c(0, 1,1,1,0,-1,-1,-1, 0))
# z[!z$x %in% -1 & !z$y %in% NA,]

# takes tentative new x, y coords, looks up close neighbors from the grid
# calculates distance to each existing point. Returns T if there's a point < min dist/F
# could go on to query surrounding cells if no collision is detected...
query_grid <- function(grid, x, y, width, height, dsep, dtest) {

  x0 <- get_index(x, dsep)
  y0 <- get_index(y, dsep)
  # ind <- x0 + (y0-1) * w*1/dsep

  # compute all 9 indices up front, starting with own cell?
  # need to check if out of bounds. If x0 is on an edge
  ind <- grid_indices(x0, y0, width, height, dsep)

  # then make this loop through the neighboring cells
  for(i in 1:length(ind)) {
    index <- ind[i]
    x_vals <- grid$xvals[[index]]
    y_vals <- grid$yvals[[index]]

    if(length(x_vals)!=0) { #return(FALSE) # or check surrounding cells?

      points <- data.frame(x_existing = x_vals, y_existing = y_vals) |>
        dplyr::mutate(dist = sqrt((x - x_existing)^2+(y-y_existing)^2))

      if(min(points$dist) < dtest) return(TRUE)
    }
  }
  return(FALSE)
}
