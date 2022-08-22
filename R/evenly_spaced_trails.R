# https://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf

#' Make evenly-spaced trails on a flowfield
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
#' @return A data frame with a 'group' column identifying each individual trail,
#'   x and y coordinates, and the grid angle associated with each coordinate
#' @export
#'
#' @examples make_flowfield() |> make_trails_even() |> draw_trails()
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

      valid <- check_valid(starting_x, starting_y,
                           width, height,
                           existing_points = dplyr::bind_rows(lines),
                           dtest = dtest)

      if(valid) {

        line_number <- line_number + 1

        new_line <- make_trail(starting_x,
                               starting_y,
                               flowfield_df,
                               existing_points = dplyr::bind_rows(lines),
                               step_length = step_length,
                               dtest = dtest)

        lines[[line_number]] <- new_line

      }
    }

    queue_position <- queue_position + 1

  }

  dplyr::bind_rows(lines, .id = "group") |>
    dplyr::mutate(group = as.numeric(group))

}



inbounds <- function(x, y, w, h) {
  x > 1 & x < w & y > 1 & y < h
}

get_seeds <- function(line, dsep) {

  n <- nrow(line)
  x <- numeric(n*2)
  y <- numeric(n*2)

  for(i in 1:n) {
    x[i]   <- line$x[i] + dsep * cos(line$a[i] + pi/2)
    x[i+n] <- line$x[i] + dsep * cos(line$a[i] - pi/2)
    y[i]   <- line$y[i] + dsep * sin(line$a[i] + pi/2)
    y[i+n] <- line$y[i] + dsep * sin(line$a[i] - pi/2)
  }

  data.frame(x = x, y = y) |>
    dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
    dplyr::arrange(rand)

}


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
