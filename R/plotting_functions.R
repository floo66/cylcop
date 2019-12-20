
#' Make scatterplot of turning angles ans steplengths
#'
#' @param traj A data.framecontaining the trajectory. Must contain collumns \code{traj$angle} and \code{traj$steplength}.
#' @param periodic A logical value denoting whether the plot should be periodically extended past -pi and pi.
#'
#' @return The scatterplot
#' @export
#'
scat_plot <- function(traj, periodic = FALSE) {
  plot_theme <- list(
    geom_hline(
      yintercept = 0,
      colour = "grey60",
      size = 0.2
    ),
    geom_point(shape="."),
    theme(legend.position = "none"),
    theme_bw(),
    xlab("X"),
    ylab(bquote(Theta)),
    theme(
      axis.title.x = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 17, colour = "black"),
      axis.text = element_text(size = 10, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    )
  )

#circular marginal density estimated with vonMises kernel density and rough bandwidth estimate
  marginal_angle_dens <- traj$angle %>%
    half2full_circ() %>%
    circular::circular(zero = 0, rotation = "counter") %>%
    circular::density.circular(
      #circular marginal density estimated with vonMises kernel density and rough bandwidth estimate
      bw = suppressWarnings(opt_circ_bw(traj$angle, loss = "adhoc", kappa.est="trigmoments")),
      kernel = "vonmises",
      na.rm = TRUE,
      control.circular = list(rotation = "counter", zero =
                                0)
    )
  marginal_angle_dens$x <- full2half_circ(marginal_angle_dens$x)
  marginal_angle_dens$y <- as.double(marginal_angle_dens$y)
  marginal_angle_dens <-
    cbind(x = marginal_angle_dens$x, y = marginal_angle_dens$y) %>%
    as.data.frame()



  if (periodic == TRUE) {

#calculate periodic images
    periodic_image <-
      traj %>%
      dplyr::filter((.data$angle <= -0.5 * pi) | (.data$angle >= 0.5 * pi)) %>%
      mutate(angle = .data$angle - (sign(.data$angle) * 2 * pi)) %>%
      mutate(period = TRUE)
    traj <- mutate(traj, period = FALSE) %>%
      rbind(periodic_image)

# main plot
    p <-
      ggplot(traj, aes(
        x = .data$steplength,
        y = as.numeric(.data$angle),
        colour = .data$period
      )) +
      geom_hline(yintercept = c(pi, -pi),
                 colour = "black",
                 size = 0.2) +
      scale_y_continuous(breaks = seq(-1.5 * pi, 1.5 * pi, 0.25 * pi),
                         labels = c(
                           expression(
                             0.5 * pi,
                             0.75 * pi,
                             -pi / pi,
                             -0.75 * pi,-0.5 * pi,
                             -0.25 * pi,
                             0,
                             0.25 * pi,
                             0.5 * pi,
                             0.75 * pi,
                             pi / -pi,
                             -0.75 * pi,
                             -0.5 * pi
                           )
                         )) +
      plot_theme +
      scale_color_manual(values = c("black", "grey"))

# add densities on the margins of the main plot
    xmarg <- cowplot::axis_canvas(p, axis = "x") +
      geom_density(
        data = traj[traj$period == FALSE,],
        aes(x = .data$steplength),
        fill = "red",
        alpha = 0.1,
        size = 0.2
      )
    ymarg <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_ribbon(
        data = marginal_angle_dens,
        aes(x = .data$x, ymin = 0, ymax = .data$y),
        alpha = 0.1,
        fill = "blue",
        color = "black"
      ) +
      coord_flip()

    p1 <-
      cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top")
    p2 <-
      cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right")
  }

#if no periodic image should be plotted
  else{
#main plot
    p <-
      ggplot(traj, aes(x = .data$steplength, y = as.numeric(.data$angle))) +
      scale_y_continuous(breaks = seq(-pi, pi, 0.25 * pi),
                         labels = c(
                           expression(
                             -pi,
                             -0.75 * pi,
                             -0.5 * pi,
                             -0.25 * pi,
                             0,
                             0.25 * pi,
                             0.5 * pi,
                             0.75 * pi,
                             pi
                           )
                         )) +

      plot_theme

# add densities on the margins of the main plot
    xmarg <- cowplot::axis_canvas(p, axis = "x") +
      geom_density(
        data = traj,
        aes(x = .data$steplength),
        fill = "red",
        alpha = 0.1,
        size = 0.2
      )
    ymarg <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_ribbon(
        data = marginal_angle_dens,
        aes(x = .data$x, ymin = 0, ymax = .data$y),
        alpha = 0.1,
        fill = "blue",
        color = "black"
      ) +
      geom_hline(yintercept = 0) +
      coord_flip()

    p1 <-
      cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top")
    p2 <-
      cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right")
  }

  return(cowplot::ggdraw(p2))
}


#' Plot the trajectories locations in x-y space
#'
#' @param traj A data.framecontaining the trajectory. Must contain collumns \code{traj$angle} and \code{traj$steplength}.
#'
#' @return The plot
#' @export
#'
traj_plot <- function(traj) {
  ggplot(traj, aes(x = .data$pos_x, y = .data$pos_y)) +
    geom_point(aes(colour = 1:length(traj$pos_x))) +
    geom_path(aes(colour = 1:length(traj$pos_x))) +
    geom_point(
      data = dplyr::slice(traj, 1, nrow(traj)),  #mark first and last point of trajectory
      aes(x = .data$pos_x, y = .data$pos_y),
      size = 4,
      color = "red"
    ) +
    scale_colour_gradientn(colours = inferno(1000)) +
    theme_bw() +
    xlab("X-position") +
    ylab("Y-position") +
    labs(color = 'Step number') +
    theme(
      axis.title = element_text(size = 12, colour = "black"),
      axis.text = element_text(size = 10, colour = "black")
    )
}



#' Make circular scatterplot of turning angles an steplengths
#'
#' @param traj A data.framecontaining the trajectory. Must contain collumns \code{traj$angle} and \code{traj$steplength}.
#'
#' @return the plot.
#' @export
#'
circ_plot <- function(traj) {

#circular marginal density estimated with vonMises kernel density and rough bandwidth estimate
  marginal_angle_dens <- traj$angle %>%
    half2full_circ() %>%
    circular::circular(zero = 0, rotation = "counter") %>%
    circular::density.circular(
      bw = suppressWarnings(opt_circ_bw(traj$angle, loss = "adhoc", kappa.est="trigmoments")),
      kernel = "vonmises",
      na.rm = TRUE,
      control.circular = list(rotation = "counter", zero =0)
    )
  marginal_angle_dens$x <- full2half_circ(marginal_angle_dens$x)
  marginal_angle_dens$y <- as.double(marginal_angle_dens$y)

#set the breaks for the gridlines of the plot
  y_breaks <- pretty(c(0, max(traj$steplength, na.rm = TRUE)), n = 5)
  y_breaks <- c(y_breaks, tail(y_breaks, 1) + y_breaks[2])
  x_breaks <- seq(-0.75 * pi, pi, 0.25 * pi)

#blow up the marginal density, so its zero is the the outermost gridline of the circular plot
  marginal_angle_dens <- cbind(x = marginal_angle_dens$x,
                               y = (0.2 * marginal_angle_dens$y + 1) * tail(y_breaks, 1)) %>%
    as.data.frame()


#set theme

  circ_plot_layers <- list(
    geom_hline( yintercept = y_breaks, colour = "grey90", size = 0.2),
    #geom_vline(xintercept = x_breaks, colour = "grey90", size = 0.2),
    coord_polar(start = pi, clip = "off"),

    scale_x_continuous(limits = c(-pi, pi),
                       breaks = x_breaks,
                       labels = c(expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
                                             0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
                                  )
      ),

    scale_y_continuous(breaks = y_breaks[1:(length(y_breaks) - 1)]),
    geom_hline(
      yintercept = y_breaks[(length(y_breaks) - 1)],
      colour = "grey60",
      size = 0.2
    ),

    geom_hline(
      yintercept = tail(y_breaks, 1),
      colour = "grey60",
      size = 0.2
    ),

    theme_bw(),

    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 10, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.border = element_blank()
    ),

    geom_segment(
      aes(
        x = x_breaks,
        y = 0,
        xend = x_breaks,
        yend = y_breaks[(length(y_breaks) - 1)]
      ),
      colour = "grey90",
      size = 0.2
    ),

    geom_segment(
      aes(
        x = x_breaks,
        y = tail(y_breaks, 1),
        xend = x_breaks,
        yend = 1.05 * max(marginal_angle_dens$y)
      ),
      colour = "grey90",
      size = 0.2
    ),

    #hjust= and vjust= don't seem to work with coord_polar()
    #neither does element_text(margin = margin())
    #The font size has different units for axes and labels, that's why I divide by 2.834646 below
    geom_label(
      data = cbind(x = -0.875 * pi,
                   y = y_breaks[2:(length(y_breaks) - 1)]) %>%
        as.data.frame(),
      aes(.data$x, .data$y, label = .data$y),
      size = 10 / 2.834646,
      label.size = 0
    )
  )


# plot

  ggplot() +
    circ_plot_layers +
    geom_point(data = traj, aes(y = .data$steplength, x = as.numeric(.data$angle))) +
    geom_line(data = marginal_angle_dens,
              aes(y = .data$y, x = .data$x),
              colour = "grey60",
              size = 0.2) +
    geom_ribbon(
      data = marginal_angle_dens,
      aes(
        x = .data$x,
        ymin = tail(y_breaks, 1),
        ymax = .data$y
      ),
      alpha = 0.3,
      fill = "blue"
    )
}


#' Make scatterplot of copula values
#'
#' @param traj #' @param traj A data.framecontaining the trajectory. Must contain collumns \code{traj$angle} and \code{traj$steplength}.
#'
#' @return The scatterplot.
#' @export
#'
cop_scat_plot <- function(traj) {
  plot_theme <- list(
    geom_point(),
    theme(legend.position = "none"),
    theme_bw(),
    xlab("v"),
    ylab("u"),
    theme(
      axis.title = element_text(size = 12, colour = "black"),
      axis.text = element_text(size = 10, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    )
  )

  p <-
    ggplot(traj, aes(x = .data$cop_v, y = .data$cop_u)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    plot_theme +
    coord_fixed()
  return(cowplot::ggdraw(p))
}



#' Make surface plot or heatmap of cdf or pdf of a copula
#'
#' @param copula A \code{cyl_copula} object.
#' @param type A string of characters, either "pdf or "cdf".
#' @param plot_type A string of characters. Available plottypes are:
#'   "rgl": surface plot,
#'   "plotly": interactive surface plot, or
#'   "ggplot": heatmap
#' @param resolution A numerical value. The density or distribution will be calculated at \code{resolution^2} points.
#' @param n_gridlines A numerical value giving the number of gridlines drawn in u and v direction.
#'
#' @return The plot.
#' @export
#'
cop_plot <- function(copula,
                     type = c("pdf", "cdf"),
                     plot_type = "rgl",
                     resolution = 50,
                     n_gridlines = 11) {

#check input

   if (type == "pdf")
    fun <- dCopula
  else if (type == "cdf")
    fun <- pCopula
  else
    stop(cylcop::error_sound(), "Type must be either cdf or pdf")

  if (!any(is(copula) == "cyl_copula") &&
      !any(is(copula) == "Copula")) {
    stop(
      cylcop::error_sound(),
      "copula must be a cyl_copula-object or a copula-object from the 'copkla'-package"
    )
  }

  if (plot_type != "rgl" &&
      plot_type != "plotly" && plot_type != "ggplot") {
    stop(cylcop::error_sound(), "plot_type must be rgl, plotly or ggplot")
  }

  printCop(copula)


# Generate matrix of values

  u <- seq(0, 1, length.out = resolution)
  v <- seq(0, 1, length.out = resolution)
  u[1] <- 0.00000001
  u[length(u)] <- 0.99999999
  v[1] <- 0.00000001
  v[length(v)] <- 0.99999999

  u_grid <- seq(0, 1, length.out = n_gridlines)
  v_grid <- seq(0, 1, length.out = n_gridlines)
  u_grid[1] <- 0.00000001
  u_grid[length(u_grid)] <- 0.99999999
  v_grid[1] <- 0.00000001
  v_grid[length(v_grid)] <- 0.99999999

  if (plot_type == "rgl" || plot_type == "plotly") {
    mat <- matrix(ncol = resolution, nrow = resolution)

    #set time for progress bar
    ptm <- proc.time()

    #calculate first and last slice of values to get an idea of the time it takes on average
    for (i in c(1, resolution)) {
      mat[i,] <- map2_dbl(u[i], v,  ~ fun(c(.x, .y), copula))
    }

    #generate progressbar
    time <- (proc.time() - ptm)[3] * resolution / 2 %>% as.integer()
    if (time > 10) {
      cat(
        "estimated time for calculation of surface: ",
        floor(time / 60),
        " minutes, ",
        time - 60 * floor(time / 60),
        " seconds\n"
      )
      pb <- utils::txtProgressBar(min = 2, max = resolution - 1)
    }

    #calculate rest of grid
    for (i in seq(2, resolution - 1)) {
      mat[i,] <- map2_dbl(u[i], v,  ~ fun(c(.x, .y), copula))
      if (time > 10){
        utils::setTxtProgressBar(pb, i)
        if(i==resolution/2) cylcop::waiting_sound()
      }
    }
    if (time > 10){
      close(pb)
      cylcop::done_sound()
    }


#------------Make surface plot with rgl------

    if (plot_type == "rgl") {

      #if gridlines are at already calculated gridvalues, no need to recalculate them
      if ((resolution - 1) %% (n_gridlines - 1) == 0) {
        mat_grid <-
          mat[seq(1, resolution, (resolution - 1) / (n_gridlines - 1)),
              seq(1, resolution, (resolution - 1) / (n_gridlines - 1))]
      }

      #calculate gridlines
      else{
        time <- time / (resolution ^ 2)

        #set progress bar
        if (time * (n_gridlines ^ 2) > 10) {
          cat(
            "gridlines are not a multiple of the surface points and have to be newly calculated.\nThis will take ",
            floor(time * (n_gridlines ^ 2) / 60),
            " minutes, ",
            (time * (n_gridlines ^ 2)) - 60 * floor(time * (n_gridlines ^ 2) / 60),
            " seconds\n"
          )
        }

        #calculate values for gridlines
        mat_grid <- matrix(ncol = n_gridlines, nrow = n_gridlines)
        for (i in seq(1, n_gridlines)) {
          mat_grid [i,] <- map2_dbl(u_grid[i], v_grid,  ~ fun(c(.x, .y), copula))
        }
      }


# Generate title

      if (any(is(copula) == "Copula")) {
        title <-
          paste(
            type,
            "of",
            copula@class[1],
            "with",
            copula@param.names,
            "=",
            copula@parameters
          )
      }
      else if (any(is(copula) == "cyl_copula")) {
        parameter_summary = paste(copula@param.names[1], " = ", copula@parameters[1])
        if (length(copula@parameters) > 1) {
          for (i in 2:length(copula@parameters)) {
            parameter_summary <-
              paste(parameter_summary,
                    ", ",
                    copula@param.names[i],
                    " = ",
                    copula@parameters[i])
          }
        }
        title <-
          paste(type, "of", copula@name, "with", parameter_summary)
      }


# Make plot

      persp3d(
        u,
        v,
        mat,
        col = inferno(50)[cut(mat, breaks = 50)],
        fit = "smooth",
        lit = TRUE,
        alpha = 0.9,
        xlab = "u",
        ylab = "v",
        zlab = if (type == "pdf")
          "c(u,v)"
        else
          "C(u,v)",
        expand = 0
      ) +
        surface3d(
          u_grid,
          v_grid,
          mat_grid,
          color = "black",
          lit = FALSE,
          front = "lines",
          back = "lines"
        ) +
        light3d(theta = 0, phi = 30) +
        title3d(main = title)
    }


#------------Make interactive surface plot with plotly------

    else if (plot_type == "plotly") {

      #set color gradient
      col <- vector("list", 50)
      for (i in 1:50) {
        col[[i]] <- c((i - 1) * (1 / 49), inferno(50)[i])
      }


# Calculate values for gridlines

      gridlines_x <- expand.grid(u = u_grid, v = v_grid)
      zg <- pmap_dbl(gridlines_x,  ~ fun(c(.x, .y), copula))
      gridlines_x <-
        mutate(gridlines_x, zg) %>% dplyr::rename(ug = u) %>% dplyr::rename(vg = v)

      #for aesthetics, add a small number so the gridlines are slightly above the surface
      gridlines_x$zg <- gridlines_x$zg + 0.01

      gridlines_y <- expand.grid(v = u_grid, u = v_grid)
      gridlines_y <- gridlines_y[c("u", "v")]
      zg <- pmap_dbl(gridlines_y,  ~ fun(c(.x, .y), copula))
      gridlines_y <-
        mutate(gridlines_y, zg) %>% dplyr::rename(ug = u) %>% dplyr::rename(vg = v)

      #for aesthetics, add a small number so the gridlines are slightly above the surface
      gridlines_y$zg <- gridlines_y$zg + 0.01


# Generate title

      if (any(is(copula) == "Copula")) {
        title <-
          paste(
            type,
            "of",
            copula@class[1],
            "with\n",
            copula@param.names,
            "=",
            copula@parameters
          )
      }
      else if (any(is(copula) == "cyl_copula")) {
        parameter_summary = paste(copula@param.names[1], " = ", copula@parameters[1])
        if (length(copula@parameters) > 1) {
          for (i in 2:length(copula@parameters)) {
            parameter_summary <-
              paste(parameter_summary,
                    ", ",
                    copula@param.names[i],
                    " = ",
                    copula@parameters[i])
          }
        }
        title <-
          paste(type, "of", copula@name, "with\n", parameter_summary)
      }


# Make plot

      p <-
        plot_ly(
          x = u,
          y = v,
          z = mat,
          showscale = F,
          colorscale = col,
          opacity = 0.9,
          lighting = list(ambient = 0.9, specular = 0)
        ) %>%
        plotly::layout(
          title = title,
          scene = list(
            xaxis = list(title = "u"),
            yaxis = list(title = "v"),
            zaxis = if (type == "pdf")
              list(title = "c(u,v)")
            else
              list(title = "C(u,v)"),
            camera = list(eye = list(
              x = -0.54,
              y = -1.8,
              z = 1.62
            ))
          )
        ) %>% add_surface()


# Add gridlines

      for (i in seq(n_gridlines + 1, ((n_gridlines - 1) * n_gridlines) + 1, n_gridlines)) {
        p <-
          add_trace(
            p,
            data = gridlines_x[i:(i + n_gridlines - 1), ],
            x = ~ vg,
            y = ~ ug,
            z = ~ zg,
            type = 'scatter3d',
            mode = 'lines',
            opacity = 0.95,
            line = list(
              width = 0.9,
              color = "white",
              reverscale = FALSE
            )
          )
      }
      for (i in seq(n_gridlines + 1, ((n_gridlines - 1) * n_gridlines) + 1, n_gridlines)) {
        p <-
          add_trace(
            p,
            data = gridlines_y[i:(i + n_gridlines - 1), ],
            x = ~ vg,
            y = ~ ug,
            z = ~ zg,
            type = 'scatter3d',
            mode = 'lines',
            opacity = 0.95,
            line = list(
              width = 0.9,
              color = "white",
              reverscale = FALSE
            )
          )
      }
      p %>% plotly::layout(showlegend = FALSE)
    }
  }

#------------Make heatmap with ggplot------
  else{

# Calculate values on grid. Different grid layout, than plotly and rgl

    outp <- expand.grid(u = u, v = v)
    outp$z <- pmap_dbl(outp,  ~ fun(c(.x, .y), copula))


# Make plot

    ggplot(outp, aes(v, u)) +
      geom_raster(aes(fill = .data$z),interpolate = T,hjust=0,vjust=0) +
      coord_fixed() +
      theme_bw() +
      labs(fill = if (type == "pdf")
        "c(u,v)"
        else
          "C(u,v)") +
      scale_fill_gradientn(colours = inferno(1000),limits=c(0, max(outp$z))) +
      theme(
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE,
        panel.grid = element_blank()
      ) +
      geom_raster(aes(fill = .data$z), alpha = 0.2) +
      geom_vline(
        xintercept = seq(0, 1, by = 0.125),
        colour = "grey60",
        size = 0.2,
        alpha = 0.25
      ) +
      geom_hline(
        yintercept = seq(0, 1, by = 0.125),
        colour = "grey60",
        size = 0.2,
        alpha = 0.25
      ) +
      scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
      scale_y_continuous(limits=c(0,1),expand=c(0,0))

  }
}
