#' Scatterplot of Turn Angles ans Step Lengths
#'
#' This function produces a scatterplot ('\code{\link[ggplot2]{ggplot}}' object) of
#' the turn angles and step lengths.
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{make_traj}()}. It must contain
#'  the columns \code{traj$angle} and \code{traj$steplength}.
#' @param periodic \link[base]{logical} value denoting whether the plot should
#' be periodically extended past -pi and pi.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object, the scatterplot.
#'
#' @examples set.seed(123)
#'
#' traj <- make_traj(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape=3)
#' )
#' plot1 <- scat_plot(traj)
#' plot2 <- scat_plot(traj, periodic = TRUE)
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cop_scat_plot}()}, \code{\link{traj_plot}()},
#' \code{\link{circ_plot}()}, \code{\link{cop_plot}()}.
#'
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
        color = "black"
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
      suppressWarnings(cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top"))
    p2 <-
      suppressWarnings(cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right"))
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
        color = "black"
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
      suppressWarnings(cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top"))
    p2 <-
      suppressWarnings(cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right"))
  }

  return(suppressWarnings(cowplot::ggdraw(p2)))
}


#' Plot a Trajectory in Euclidean Space
#'
#' This function plots the locations of a trajectory. The first and last point are
#' marked in red.
#'
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{make_traj}()}. It must contain
#'  the columns \code{traj$pos_x} and \code{traj$pos_y}.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object.
#'
#' @examples set.seed(123)
#'
#' traj <- make_traj(50,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape=3)
#' )
#' plot1 <- traj_plot(traj)
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cop_scat_plot}()},
#' \code{\link{circ_plot}()}, \code{\link{cop_plot}()}, \code{\link{scat_plot}()}.
#'
#' @export
#'
traj_plot <- function(traj) {
  p <- ggplot(traj, aes(x = .data$pos_x, y = .data$pos_y)) +
    geom_point(aes(colour = 1:length(traj$pos_x)),
               size=min(4,100/nrow(traj))) +
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
  suppressWarnings(plot(p))
}




#' Circular Scatterplot of Turn Angles and Step Lengths
#'
#' This function produces a circular scatterplot with the step lengths plotted
#' as distance from the center of a circle and the turn angles as angles
#' (polar coordinates).
#'
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{make_traj}()}. It must contain
#'  the columns \code{traj$angle} and \code{traj$steplength}.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object.
#'
#' @examples set.seed(123)
#'
#' traj <- make_traj(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull", list(shape=3)
#' )
#' plot1 <- circ_plot(traj)
#'
#' @references \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cop_scat_plot}()}, \code{\link{traj_plot}()},
#' \code{\link{cop_plot}()}, \code{\link{scat_plot}()}.
#'
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
  y_breaks <- pretty(c(0, max(traj$steplength, na.rm = TRUE)), n = 6)
  y_breaks <- c(y_breaks, tail(y_breaks, 1) + y_breaks[2])
  x_breaks <- seq(-0.75 * pi, pi, 0.25 * pi)

  #blow up the marginal density, so its zero is the the outermost gridline of the circular plot
  marginal_angle_dens <- cbind(x = marginal_angle_dens$x,
                               y = (0.6 * marginal_angle_dens$y + 1) * tail(y_breaks, 1)) %>%
    as.data.frame()

  #set positions of radius labels
  radius_label_pos <-  cbind(x = rep(c(-1,1),length(y_breaks))[2:(length(y_breaks) - 1)]*0.875 * pi,
                             y = y_breaks[2:(length(y_breaks) - 1)]) %>%
    as.data.frame()

  #set positions of angle labels
  max_rad <- max(y_breaks)
  angle_label_pos <-  cbind(x = x_breaks,
                            y = c(1.2, 1.2, 1.2, 1.15, 1.2, 1.2, 1.2, 1.15)*marginal_angle_dens$y[map(x_breaks,~which.min(abs(marginal_angle_dens$x-.x)))%>%unlist()]) %>%
    as.data.frame()

  #convert to cartesian coordinates to determine plot-margins
  angle_label_pos_cart <- angle_label_pos
  angle_label_pos_cart$x <- angle_label_pos$y*sin(angle_label_pos$x)
  angle_label_pos_cart$y <- angle_label_pos$y*cos(angle_label_pos$x)

  panel_size <- max(max(abs(angle_label_pos_cart$x)),max(abs(angle_label_pos_cart$y)))
  top_correction <- panel_size-abs(max(angle_label_pos_cart$y))
  right_correction <- panel_size-abs(max(angle_label_pos_cart$x))
  bottom_correction <- panel_size-abs(min(angle_label_pos_cart$y))
  left_correction <- panel_size-abs(min(angle_label_pos_cart$x))


  #set theme

  circ_plot_layers <- list(
    geom_hline( yintercept = y_breaks[seq(2,length(y_breaks),2)], colour = "darkgreen", size = 0.2),
    geom_hline( yintercept = y_breaks[seq(1,length(y_breaks),2)], colour = "darkred", size = 0.2),
    #geom_vline(xintercept = x_breaks, colour = "grey90", size = 0.2),
    coord_polar(start = pi, clip = "off"),

    # scale_x_continuous(limits = c(-pi, pi),
    # breaks = x_breaks,
    # labels = c(expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
    # 0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
    # )
    # ),

    # scale_y_continuous(breaks = y_breaks[1:(length(y_breaks) - 1)]),
    # geom_hline(
    #   yintercept = y_breaks[(length(y_breaks) - 1)],
    #   colour = "darkred",
    #   size = 0.3
    # ),

    geom_hline(
      yintercept = tail(y_breaks, 1),
      colour = "grey30",
      size = 0.2
    ),

    theme_bw(),

    theme(
      panel.grid = element_blank(),
      panel.border=element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
     # plot.margin = unit(c(-10-5*top_correction, -10-5*right_correction, -10-5*bottom_correction, -10-5*left_correction), "pt")
    ),

    geom_segment(
      aes(
        x = x_breaks,
        y = 0,
        xend = x_breaks,
        yend = y_breaks[(length(y_breaks) - 1)]
      ),
      colour = "grey30",
      size = 0.2
    ),

    geom_segment(
      aes(
        x = x_breaks,
        y = tail(y_breaks, 1),
        xend = angle_label_pos$x,
        yend = angle_label_pos$y
      ),
      colour = "grey30",
      size = 0.2
    ),

    #hjust= and vjust= don't seem to work with coord_polar()
    #neither does element_text(margin = margin())
    #The font size has different units for axes and labels, that's why I divide by 2.834646 below
    geom_label(
      data = radius_label_pos[seq(1,nrow(radius_label_pos),2),],
      aes(.data$x, .data$y, label = .data$y),
      size = 10 / 2.834646,
      color ="darkgreen",
      label.size = 0
    ),
    geom_label(
      data = radius_label_pos[seq(2,nrow(radius_label_pos),2),],
      aes(.data$x, .data$y, label = .data$y),
      size = 10 / 2.834646,
      color ="darkred",
      label.size = 0
    ),

    geom_label(
      data = angle_label_pos,
      aes(.data$x, .data$y),
      label = c(expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
                           0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
      ),
      size = 10 / 2.834646,
      label.size = 0
    )
  )


  # plot

  p <- ggplot() +
    circ_plot_layers +
    geom_point(data = traj, aes(y = .data$steplength, x = as.numeric(.data$angle)),alpha=0.5,size=0.3) +
    geom_line(data = marginal_angle_dens,
              aes(y = .data$y, x = .data$x),
              colour = "black",
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
  suppressWarnings(plot(p))
}


#' Scatterplot of Copula Values
#'
#' This function produces a scatterplot ('\code{\link[ggplot2]{ggplot}}' object) of
#' a sample from a copula. Either a sample is provided as input, or a sample of
#' 10000 points is drawn from a copula to quickly visualize it.
#'
#' @param input Either
#' a \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{make_traj}()}, which must contain columns
#' \code{traj$cop_u} and \code{traj$cop_v},
#' or a '\code{\linkS4class{cyl_copula}}' object or a '\code{\linkS4class{Copula}}' object
#' of the package '\pkg{copula}'.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object, the scatterplot.
#'
#' @examples set.seed(123)
#'
#' traj <- make_traj(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape=3)
#' )
#' cop_scat_plot(traj)
#' cop_scat_plot(cyl_quadsec(0.1))
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{traj_plot}()},
#' \code{\link{circ_plot}()}, \code{\link{cop_plot}()}, \code{\link{scat_plot}()}.
#'
#' @export
#'
cop_scat_plot <- function(input) {

  if ((any(is(input) == "cyl_copula"))||(any(is(input) == "Copula"))){
    sample <- rcylcop(10000,input)
    traj <- data.frame(cop_u=sample[,1], cop_v=sample[,2])
  } else if (is.data.frame(input)){
    if(!all(c("cop_u","cop_v") %in% colnames(input))){
      stop(error_sound(), "trajectory must contain the columns 'cop_u' and 'cop_v'")
    }
    traj <- input
  }
  else{
    stop(error_sound(), "Provide either a (cyl_)copula object or a trajectory data.frame")
  }

  plot_theme <- list(
    geom_point(size=0.01,alpha=0.5),
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
  return(suppressWarnings(cowplot::ggdraw(p)))
}



#' Surface Plot or Heat Map of the Distribution or the Density of a Copula
#'
#' This function plots the distribution or the density of a copula. It can produce
#' a surface plot using either functions from the '\pkg{rgl}' or from the
#' '\pkg{plotly}' package, or it can produce a heat map using functions from
#' '\pkg{ggplot2}'.
#'
#' @param copula '\code{\linkS4class{cyl_copula}}' or a '\code{\linkS4class{Copula}}' object
#' from the package '\pkg{copula}'.
#' @param type \link[base]{character} string describing what is plotted,
#' either \code{"pdf"} or \code{"cdf"}.
#' @param plot_type \link[base]{character} string describing what type of plot
#' is produced. Available plot types are:
#'   \code{"rgl"}: surface plot,
#'   \code{"plotly"}: interactive surface plot, or
#'   \code{"ggplot"}: heatmap
#' @param resolution \link[base]{numeric} value. The density or distribution
#' will be calculated at \code{resolution^2} points.
#' @param n_gridlines \link[base]{numeric} value giving the number of grid
#' lines drawn in u and v direction.
#'
#' @return Depending on \code{plot_type}, a '\code{\link[ggplot2]{ggplot}}' object
#' is returned, or a '\pkg{plotly}' visualization or '\pkg{rgl}' plot is produced.
#'
#' @examples
#' cop_plot(copula::frankCopula(2), type="pdf", plot_type="ggplot")
#' cop_plot(copula::frankCopula(2), type="cdf", plot_type="ggplot")
#' cop_plot(copula::frankCopula(2), type="pdf", plot_type="ggplot", resolution = 5)
#'
#' #opens a new window
#' cop_plot(cyl_quadsec(0.1), type="pdf", plot_type="rgl")
#' cop_plot(cyl_quadsec(0.1), type="pdf", plot_type="rgl", n_gridlines = 60)
#'
#' cop_plot(cyl_quadsec(0.1), type="pdf", plot_type="plotly", n_gridlines = 20)
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cop_scat_plot}()}, \code{\link{traj_plot}()},
#' \code{\link{circ_plot}()}, \code{\link{scat_plot}()}.
#'
#' @export
#'
cop_plot <- function(copula,
                     type = c("pdf", "cdf"),
                     plot_type = "rgl",
                     resolution = 50,
                     n_gridlines = 11) {

#check input

   if (type == "pdf")
    fun <- dcylcop
  else if (type == "cdf")
    fun <- pcylcop
  else
    stop(error_sound(), "Type must be either cdf or pdf")

  if (!any(is(copula) == "cyl_copula") &&
      !any(is(copula) == "Copula")) {
    stop(
      error_sound(),
      "copula must be a cyl_copula-object or a copula-object from the 'copula'-package"
    )
  }

  if (plot_type != "rgl" &&
      plot_type != "plotly" && plot_type != "ggplot") {
    stop(error_sound(), "plot_type must be rgl, plotly or ggplot")
  }
  if(cylcop.env$silent==F)  printCop(copula)


# Generate matrix of values

  u <- seq(0, 1, length.out = resolution)
  v <- seq(0, 1, length.out = resolution)

  #only necessary when we use the copula::dCopula generic of the copula package and not our own dcylcop generic
  # u[1] <- 0.00000001
  # u[length(u)] <- 0.99999999
  # v[1] <- 0.00000001
  # v[length(v)] <- 0.99999999

  u_grid <- seq(0, 1, length.out = n_gridlines)
  v_grid <- seq(0, 1, length.out = n_gridlines)
  #only necessary when we use the copula::dCopula generic of the copula package and not our own dcylcop generic
  # u_grid[1] <- 0.00000001
  # u_grid[length(u_grid)] <- 0.99999999
  # v_grid[1] <- 0.00000001
  # v_grid[length(v_grid)] <- 0.99999999

  if (plot_type == "rgl" || plot_type == "plotly") {

    mat <- as.matrix(expand.grid(u = u, v = v))

    #set time for progress bar
    time <- 0
    ptm <- proc.time()

    #calculate first and last slice of values to get an idea of the time it takes on average

    outp <- rep(NA, resolution^2)
    outp[1:resolution] <- fun(mat[1:resolution,],copula)
    outp[(nrow(mat)-resolution+1):nrow(mat)] <- fun(mat[(nrow(mat)-resolution+1):nrow(mat),],copula)


    if(cylcop.env$silent==F){
      #generate progressbar
      time <- (proc.time() - ptm)[3] * resolution / 2 %>% as.integer()
      if (time > 10) {
        message(
          "estimated time for calculation of surface: ",
          floor(time / 60),
          " minutes, ",
          round(time - 60 * floor(time / 60)),
          " seconds\n"
        )
        pb <- utils::txtProgressBar(min = 2, max = resolution - 1)
      }
    }

    #calculate rest of grid
    for (i in seq(2, resolution - 1)) {
      outp[((i-1)*resolution+1):(i*resolution)] <- fun(mat[((i-1)*resolution+1):(i*resolution),],copula)
      if (time > 10 && cylcop.env$silent==F){
        utils::setTxtProgressBar(pb, i)
        if(i==resolution/2) waiting_sound()
      }
    }
    if (time > 10 && cylcop.env$silent==F){
      close(pb)
      done_sound()
    }

    quant_995 <- quantile(outp,probs = 0.995)
    max_outp <- max(outp)
    if(max_outp>(2*quant_995)){
      outp[which(outp>quant_995)] <- quant_995
      warning(
        warning_sound(),
        paste0("Maximum of ", type, " is ", round(max_outp,2),
               "; for better visulatization,\n values larger than ",
               round(quant_995,2), " (99.5 percentile) are cut off.")
      )
    }

    mat <- matrix(outp, ncol = resolution, nrow=resolution, byrow = F)



    #if gridlines are at already calculated gridvalues, no need to recalculate them
    if ((resolution - 1) %% (n_gridlines - 1) == 0) {
      if(n_gridlines==0) mat_grid <- u_grid else{
        mat_grid <- mat[seq(1, resolution, (resolution - 1) / (n_gridlines - 1)),
            seq(1, resolution, (resolution - 1) / (n_gridlines - 1))]
      }
    }

    #calculate gridlines
    else{
      time <- time / (resolution ^ 2)

      #Calculate timing
      if (time * (n_gridlines ^ 2) > 10 && cylcop.env$silent==F) {
        message(
          "gridlines are not a multiple of the surface points and have to be newly calculated.\nThis will take ",
          floor(time * (n_gridlines ^ 2) / 60),
          " minutes, ",
          round((time * (n_gridlines ^ 2)) - 60 * floor(time * (n_gridlines ^ 2) / 60)),
          " seconds"
        )
      }

      #calculate values for gridlines
      mat_grid <- as.matrix(expand.grid(u = u_grid, v = v_grid))
      quant_995 <- quantile(outp,probs = 0.995)
      outp_grid <- fun(mat_grid,copula)
      if(max_outp>(2*quant_995)){
        outp_grid[which(outp_grid>quant_995)] <- quant_995
      }
      mat_grid <- matrix(outp_grid, ncol = n_gridlines, nrow=n_gridlines, byrow = F)
    }

    #------------Make surface plot with rgl------

    if (plot_type == "rgl") {
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
                    round(copula@parameters[i],3))
          }
        }
        title <-
          paste(type, "of", copula@name, "with", parameter_summary)
      }


      # Make plot

      p <-  persp3d(
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
      ) +      light3d(theta = 0, phi = 30) +
        title3d(main = title)+
        (if(n_gridlines>0){surface3d(
          u_grid,
          v_grid,
          mat_grid,
          color = "black",
          lit = FALSE,
          front = "lines",
          back = "lines"
        )} )
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
      zg <- c(mat_grid)
      gridlines_x <-
        mutate(gridlines_x, zg) %>% dplyr::rename(ug = u) %>% dplyr::rename(vg = v)

      #for aesthetics, add a small number so the gridlines are slightly above the surface
      gridlines_x$zg <- gridlines_x$zg + 0.01

      gridlines_y <- expand.grid(v = u_grid, u = v_grid)
      gridlines_y <- gridlines_y[c("u", "v")]
      zg <- c(t(mat_grid))
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
                    round(copula@parameters[i],3))
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
          showscale = FALSE,
          colorscale = col,
          opacity = 1,
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
      if(n_gridlines>0){
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
              opacity = 1,
              line = list(
                width = 0.5,
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
              opacity = 1,
              line = list(
                width = 0.9,
                color = "white",
                reverscale = FALSE
              )
            )
        }

        gridlines_x$zg <- gridlines_x$zg - 0.02
        gridlines_y$zg <- gridlines_y$zg - 0.02
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
              opacity = 1,
              line = list(
                width = 0.5,
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
              opacity = 1,
              line = list(
                width = 0.9,
                color = "white",
                reverscale = FALSE
              )
            )
        }
      }

      #turn warnings back to original state

      p <- p %>% plotly::layout(showlegend = FALSE)
      suppressWarnings(print(p))
    }
  }

#------------Make heatmap with ggplot------
  else{

# Calculate values on grid. Different grid layout, than plotly and rgl

    outp <- expand.grid(u = u, v = v)
    #outp$z <- pmap_dbl(outp,  ~ fun(c(.x, .y), copula))
    outp$z <-  fun(as.matrix(outp[,1:2]), copula)
    quant_995 <- quantile(outp$z,probs = 0.995)
    max_outp <- max(outp$z)
    if(max_outp>(2*quant_995)){
      outp$z[which(outp$z>quant_995)] <- quant_995
      warning(
        warning_sound(),
        paste0("Maximum of ", type, " is ", round(max_outp,2),
               "; for better visulatization,\n values larger than ",
               round(quant_995,2), " (99.5 percentile) are cut off.")
      )
    }

# Make plot

    p <- ggplot(outp, aes(v, u)) +
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
        panel.grid = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")

      ) +
      geom_raster(aes(fill = .data$z), alpha = 0.2) +
      geom_vline(
        xintercept = seq(0, 1, by = 0.125),
        colour = "grey60",
        size = 0.8,
        alpha = 0.25
      ) +
      geom_hline(
        yintercept = seq(0, 1, by = 0.125),
        colour = "grey60",
        size = 0.8,
        alpha = 0.25
      ) +
      scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
      scale_y_continuous(limits=c(0,1),expand=c(0,0))
    suppressWarnings(plot(p))
  }
}
