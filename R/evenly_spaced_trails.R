# https://web.cs.ucdavis.edu/~ma/SIGGRAPH02/course23/notes/papers/Jobard.pdf

make_trail <- function(w = 50, h = 50, lines, steps = 10, step_length = .5, d = 1) {

  dsep <- w*.05

  out <- data.frame(x = 0, y = 0, a = 0, line = 0)
  qtree <- update_tree(out)
  # for the first line, just pic a random starting point
  starting_x <- w/2 #runif(1, min = 1, max = w-1)
  starting_y <- h/2 #runif(1, min = 1, max = h-1)
  line <- 1
  new_line <- grow_line(starting_x, starting_y, out, qtree, ff, steps, step_length, d, line = line, w, h)
  # grid <- add_to_grid(grid, new_line$x, new_line$y)
  out <- new_line
  qtree <- update_tree(out)
  queue <- list(new_line)
  seeds <- update_current_line(queue, dsep, w, h)
  queue <- queue[-1]
  finished <- FALSE

  while(!finished) {

    # use seedpoints from line from queue
    for(i in 1:length(seeds$x)) {

      starting_x <- seeds$x[i]
      starting_y <- seeds$y[i]

      # check if the point is valid
      valid <- !check_neighbor(out, qtree, starting_x, starting_y, test_distance = d)

      if(valid) {
        # message("valid")
        line <- line + 1
        new_line <- grow_line(starting_x, starting_y, out, qtree, ff, steps, step_length, d, line = line, w, h)
        # grid <- add_to_grid(grid, new_line$x, new_line$y)
        out <- dplyr::bind_rows(out, new_line)
        qtree <- update_tree(out)
        queue[[length(queue)+1]] <- new_line
      }
    }

    # finished = T
    # message(length(queue))
    # once we tried all seeds, need to update the seeds
    if(length(queue)==0) break
    seeds <- update_current_line(queue, dsep, w, h)
    queue <- queue[-1]

  }

  return(out)

}



inbounds <- function(x, y, w, h) {
  x > 0 & x < w & y > 0 & y < h
}

# returns TRUE if the nearest neighbor is too close
check_neighbor <- function(data, qtree, x, y, test_distance) {
  index <- SearchTrees::knnLookup(qtree, x, y, k = 1)
  dist <- sqrt((x-data$x[index])^2+(y-data$y[index])^2)
  dist < test_distance
}

update_tree <- function(new_data) {
  SearchTrees::createTree(new_data)
}


grow_line <- function(starting_x, starting_y, data, qtree, ff, steps, step_length, d, line, w, h) {

  x <- starting_x
  x2 <- starting_x
  y <- starting_y
  y2 <- starting_y
  a <- ff$angle[ceiling(x[1]) + (ceiling(y[1])-1) * w]
  a2 <- a

  for(i in 2:steps) {
    prev_x <- x[i-1]
    prev_y <- y[i-1]
    angle <- ff$angle[ceiling(x[i-1]) + (ceiling(y[i-1])-1) * w]
    new_x <- prev_x + step_length * cos(angle)
    new_y <- prev_y + step_length * sin(angle)

    # if it reaches the edges, we're done
    if(!inbounds(new_x, new_y, w, h)) break

    # if it's too close to an existing line, we're done
    if(line > 1){
      if(check_neighbor(data, qtree, new_x, new_y, test_distance = d)) break
    }
    x[i] <- new_x
    y[i] <- new_y
    a[i] <- angle
  }

  # moving the other way
  for(h in 2:steps) {
    prev_x <- x2[h-1]
    prev_y <- y2[h-1]
    angle <- ff$angle[ceiling(x2[h-1]) + (ceiling(y2[h-1])-1) * w]
    new_x <- prev_x + step_length * cos(angle+pi)
    new_y <- prev_y + step_length * sin(angle+pi)
    # if it reaches the edges, we're done
    if(!inbounds(new_x, new_y, w, h)) break

    # if it's too close to an existing line, we're done
    if(line > 1) {
      if(check_neighbor(data, qtree, new_x, new_y, test_distance = d)) break
    }
    x2[h] <- new_x
    y2[h] <- new_y
    a2[h] <- angle

  }

  x <- c(rev(x2), x)
  y <- c(rev(y2), y)
  a <- c(rev(a2), a)

  data.frame(x, y, a, line)

}




# create new seedpoints
update_current_line <- function(queue, dsep, w, h) {
  # take the first line in the queue (the oldest line) as the new current line
  line <- queue[[1]]
  n <- nrow(line)
  x <- numeric()
  y <- numeric()
  # now compute all possible seed points; for each row of the df
  for(i in 1:n) {

    x[i]   <- line$x[i] + dsep * cos(line$a[i] + pi/2)
    x[i+n] <- line$x[i] + dsep * cos(line$a[i] - pi/2)
    y[i]   <- line$y[i] + dsep * sin(line$a[i] + pi/2)
    y[i+n] <- line$y[i] + dsep * sin(line$a[i] - pi/2)
  }

  all_possible_seeds <- data.frame(x = x, y = y) |>
    dplyr::filter(inbounds(x, y, w, h)) |>
    dplyr::mutate(rand = sample(1:dplyr::n(), dplyr::n())) |>
    dplyr::arrange(rand)

  return(all_possible_seeds)
}


# w <- 50
# h <- 50
# dsep <- w*.05 #w*.1 # can be .1, .05
# dtest <- dsep*.9
#
#
ff <- make_flowfield(50, 50, angle = 1, octaves = 1, frequency = .01, plot = FALSE)
#
t <- make_trail(lines = 1, steps = 20, step_length = 1, d = .5) |>
  tibble::rownames_to_column("point")

ggplot2::ggplot(t, ggplot2::aes(x, y, group = line)) +
  ggplot2::geom_path() +
  ggplot2::coord_fixed()

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
