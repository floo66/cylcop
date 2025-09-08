#' Scatterplot of Turn Angles and Step Lengths
#'
#' This function produces a scatterplot ('\code{\link[ggplot2]{ggplot}}' object) of
#' the turn angles and step lengths.
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{traj_sim}()}. It must contain
#'  the columns \code{traj$angle} and \code{traj$steplength}.
#' @param theta (alternatively) \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable).
#' @param x (alternatively) \link[base]{numeric} \link[base]{vector} of
#' measurements of a linear variable.
#' @param periodic \link[base]{logical} value denoting whether the plot should
#' be periodically extended past -pi and pi.
#' @param margins \link[base]{logical} determining whether the marginal kernel
#' density estimates are computed and plotted. Alternatively, \code{margins} can
#' be a list of length 2 containing a kernel density estimate for \code{theta} and
#' a kernel density estimate for \code{x}. One entry must be of type
#' \code{'density.circular'} (as returned e.g. by \code{\link{fit_circ_np}(theta))},
#' and the other entry must be of type \code{"density"}
#' (as returned e.g. by \code{\link{fit_lin_np}(x))}.
#'
#' @details You can either specify \code{traj} or the angles and step lengths
#' (\code{theta} and \code{x}).
#' If \code{margins=T}, the code will attempt to find appropriate bandwidths for
#' the kernel density estimate autonomously, also taking into account computational time.
#' For more control over the actual method and parameters used to obtain the kernel
#' density estimates, you can calculate them "by hand" using e.g.
#' \code{\link{fit_circ_np}(theta))}
#' and \code{\link{fit_lin_np}(x))}.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object, the scatterplot.
#'
#' @examples set.seed(123)
#' traj <- traj_sim(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_lin = list(name = "weibull", coef = list(shape = 3))
#' )
#'
#' plot1 <- plot_joint_scat(traj)
#' plot2 <- plot_joint_scat(traj, periodic = TRUE)
#' plot3 <- plot_joint_scat(theta=traj$angle, x=traj$steplength, periodic = TRUE, margins=TRUE)
#'
#' bw <- opt_circ_bw(theta = traj$angle, method = "nrd",kappa.est = "trigmoments")
#' ang_dens <- fit_circ_np(theta=traj$angle, bandwidth=bw)
#' step_dens <- fit_lin_np(x=traj$steplength)
#' plot4 <- plot_joint_scat(traj, periodic = TRUE, margins=list(ang_dens, step_dens))
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_cop_scat}()}, \code{\link{plot_track}()},
#' \code{\link{plot_joint_circ}()}, \code{\link{plot_cop_surf}()}.
#'
#' @export
#'
plot_joint_scat <-
  function(traj = NULL,
           theta = NULL,
           x = NULL,
           periodic = FALSE,
           margins = FALSE) {
    #validate input
    tryCatch({
      check_arg_all(list(
        check_argument_type(traj,
                            type = c("data.frame", "list")),
        check_argument_type(traj,
                            type = "NULL")
      )
      , 2)
      check_arg_all(list(
        check_argument_type(theta,
                            type = "numeric"),
        check_argument_type(theta,
                            type = "NULL")
      )
      , 2)
      check_arg_all(list(
        check_argument_type(x,
                            type = "numeric"),
        check_argument_type(x,
                            type = "NULL")
      )
      , 2)
      check_arg_all(check_argument_type(periodic,
                                        type = "logical")
                    , 1)
      check_arg_all(list(
        check_argument_type(margins,
                            type = "logical"),
        check_argument_type(margins,
                            type = "list")
      )
      , 2)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    })
    if ((is.null(theta) &&
         !is.null(x)) || (!is.null(theta) && is.null(x))) {
      stop(error_sound(), "Specify angles AND step lengths!")
    }
    if (is.null(traj) && (is.null(theta) || is.null(x))) {
      stop(error_sound(),
           "Specify either the trajectory or angles and steplengths!")
    }
    if (!is.null(traj) && (!is.null(theta) || !is.null(x))) {
      stop(error_sound(),
           "Specify either the trajectory or angles and steplengths!")
    }
    if (length(theta) != length(x)) {
      stop(error_sound(), "x and theta must have the same length!")
    }


    if (!is.null(theta)) {
      traj <- data.frame(angle = theta, steplength = x)
    }

    if (!is.list(margins)) {
      calc_margins <- margins
      plot_margins <- margins
    } else{
      marginal_angle_dens <- NULL
      marginal_step_dens <- NULL

      for (i in 1:2) {
        if (any(is(margins[[i]]) == "density.circular")) {
          marginal_angle_dens <- margins[[i]]
        } else if (any(is(margins[[i]]) == "density")) {
          marginal_step_dens <- margins[[i]]
        } else {
          stop(
            error_sound(),
            "If a list of densities is provided with input 'margins', the elements
      of that list must be of type 'density.circular' or 'density'."
          )
        }
      }

      # Check if both variables are assigned
      if (is.null(marginal_angle_dens) || is.null(marginal_step_dens)) {
        stop(
          error_sound(),
          "If a list of densities is provided with input 'margins', one entry
      of that list must be of type 'density.circular' and the other of type 'density'."
        )
      }
      plot_margins <- TRUE
      calc_margins <- FALSE
    }

    plot_theme <- list(
      geom_hline(
        yintercept = 0,
        colour = "grey60",
        size = 0.2
      ),
      geom_point(shape = 16, alpha = 0.5),
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

    #circular marginal density estimated with von Mises kernel density
    if (calc_margins) {
      #rough bandwidth estimate with trigonometric moments. This can fail when
      #the function has no root (i.e. uniroot gives an error) in this case use
      #only a subset of the sample (for speed) and maximize the cross validation
      #log–likelihood with respect to the bandwidth
      bw <-  tryCatch({
        suppressWarnings(opt_circ_bw(traj$angle, method = "nrd", kappa.est = "trigmoments"))
      },
      error = function(e) {
        if (length(traj$angle) > 1000) {
          sample_angle <-
            traj$angle[sample(seq_along(traj$angle), 1000, replace = FALSE)]
        } else{
          sample_angle <- traj$angle
        }
        suppressWarnings(opt_circ_bw(sample_angle, method = "cv"))
      })

      marginal_angle_dens <- traj$angle %>%
        half2full_circ() %>%
        circular::circular(zero = 0, rotation = "counter") %>%
        circular::density.circular(
          #circular marginal density estimated with von Mises kernel density and rough bandwidth estimate
          bw = bw,
          kernel = "vonmises",
          na.rm = TRUE,
          control.circular = list(rotation = "counter", zero =
                                    0)
        )

      if(min(traj$steplength)>0){
        lowlim <- 0
      }else{
        lowlim <- -Inf
      }
      marginal_step_dens <-
        fit_lin_np(traj$steplength, limits=c(lowlim,Inf))

    }
    if (plot_margins) {
      marginal_angle_dens$x <- full2half_circ(marginal_angle_dens$x)
      marginal_angle_dens$y <- as.double(marginal_angle_dens$y)
      marginal_angle_dens <-
        cbind(x = marginal_angle_dens$x, y = marginal_angle_dens$y) %>%
        as.data.frame()
      marginal_step_dens <-
        cbind(x = marginal_step_dens$x, y = marginal_step_dens$y) %>%
        as.data.frame()
    }


    if (periodic == TRUE) {
      #calculate periodic images
      periodic_image <-
        traj %>%
        dplyr::filter((.data$angle <= -0.5 * pi) |
                        (.data$angle >= 0.5 * pi)) %>%
        mutate(angle = .data$angle - (sign(.data$angle) * 2 * pi)) %>%
        mutate(period = TRUE)
      traj <- mutate(traj, period = FALSE) %>%
        rbind(periodic_image)

      # main plot
      p <-
        ggplot(traj,
               aes(
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
                           ),expand=c(0.01,0.01)) +
        plot_theme +
        scale_color_manual(values = c("black", "grey60"))

      # add densities on the margins of the main plot
      if (any(plot_margins != FALSE)) {
        xmarg <- cowplot::axis_canvas(p, axis = "x") +
          geom_ribbon(
            data = marginal_step_dens,
            aes(
              x = .data$x,
              ymin = 0,
              ymax = .data$y
            ),
            alpha = 0.1,
            fill = "red",
            color = "black"
          )
        ymarg <-
          cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) +
          geom_ribbon(
            data = marginal_angle_dens,
            aes(
              x = .data$x,
              ymin = 0,
              ymax = .data$y
            ),
            alpha = 0.1,
            fill = "blue",
            color = "black"
          ) +
          coord_flip()

        p1 <-
          suppressWarnings(cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top"))
        p <-
          suppressWarnings(cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right"))
      }
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
      if (any(plot_margins != FALSE)) {
        xmarg <- cowplot::axis_canvas(p, axis = "x") +
          geom_ribbon(
            data = marginal_step_dens,
            aes(
              x = .data$x,
              ymin = 0,
              ymax = .data$y
            ),
            alpha = 0.1,
            fill = "red",
            color = "black"
          )
        ymarg <-
          cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE) +
          geom_ribbon(
            data = marginal_angle_dens,
            aes(
              x = .data$x,
              ymin = 0,
              ymax = .data$y
            ),
            alpha = 0.1,
            fill = "blue",
            color = "black"
          ) +
          geom_hline(yintercept = 0) +
          coord_flip()

        p1 <-
          suppressWarnings(cowplot::insert_xaxis_grob(p, xmarg, grid::unit(.2, "null"), position = "top"))
        p <-
          suppressWarnings(cowplot::insert_yaxis_grob(p1, ymarg, grid::unit(.2, "null"), position = "right"))
      }
    }

    return(suppressWarnings(cowplot::ggdraw(p)))
  }


#' Plot a Trajectory in Euclidean Space
#'
#' This function plots the locations of a trajectory or multiple trajectories.
#'
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{traj_sim}()}. It must contain
#'  the columns \code{traj$pos_x} and \code{traj$pos_y}. It is also possible to specify a
#'  \link[base]{list} of such data.frames containing multiple trajectories.
#' @param x_coord (alternatively) \link[base]{numeric} \link[base]{vector} of x-coordinates or
#' a \link[base]{list} of x-coordinate vectors of multiple trajectories.
#' @param y_coord (alternatively) \link[base]{numeric} \link[base]{vector} of y-coordinates or
#' a \link[base]{list} of y-coordinate vectors of multiple trajectories.
#'
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object.
#'
#' @examples set.seed(123)
#' traj <- traj_sim(50,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_lin = list(name = "weibull", coef = list(shape = 3))
#' )
#' plot1 <- plot_track(traj=traj)
#'
#' x_coord <- list(runif(1000),runif(20),runif(3))
#' y_coord <- list(runif(1000),runif(20),runif(3))
#'
#' plot2 <- plot_track(x_coord=x_coord, y_coord=y_coord)
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_cop_scat}()},
#' \code{\link{plot_joint_circ}()}, \code{\link{plot_cop_surf}()}, \code{\link{plot_joint_scat}()}.
#'
#' @export
#'
plot_track <- function(traj = NULL,
                      x_coord = NULL,
                      y_coord = NULL
                      ) {
  tryCatch({
    check_arg_all(list(
      check_argument_type(traj,
                          type = c("data.frame", "list")),
      check_argument_type(traj,
                          type = "NULL"),
      check_argument_type(traj,
                          type = "list")
    )
    ,
    3)
    check_arg_all(list(
      check_argument_type(x_coord,
                          type = "numeric"),
      check_argument_type(x_coord,
                          type = "NULL"),
      check_argument_type(x_coord,
                          type = "list")
    )
    ,
    3)
    check_arg_all(list(
      check_argument_type(y_coord,
                          type = "numeric"),
      check_argument_type(y_coord,
                          type = "NULL"),
      check_argument_type(y_coord,
                          type = "list")
    )
    ,
    3)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })
  if ((is.null(x_coord) &&
       !is.null(y_coord)) ||
      (!is.null(x_coord) && is.null(y_coord))) {
    stop(error_sound(), "Specify x-coordinates AND y-coordinates")
  }
  if (is.null(traj) && (is.null(x_coord) || is.null(y_coord))) {
    stop(error_sound(),
         "Specify either the trajectory or x-coordinates and y-coordinates")
  }
  if (!is.null(traj) && (!is.null(x_coord) || !is.null(y_coord))) {
    stop(error_sound(),
         "Specify either the trajectory OR x-coordinates and y-coordinates")
  }


  #validate input

  if (!is.null(traj)) {
    if (is.data.frame(traj)) {
      traj_lst <- list(traj)
    } else{
      traj_lst <- traj
    }
    for (i in seq_along(traj_lst)) {
      if (!all(c("pos_x", "pos_y") %in% names(traj_lst[[i]]))) {
        stop(
          error_sound(),
          "traj-data.frames must contain the columns traj$pos_x and traj$pos_y."
        )
      }
      traj_lst[[i]] <-
        traj_lst[[i]][, (names(traj_lst[[i]])) %in% c("pos_x", "pos_y")]
    }
  }

  if (!is.null(y_coord)) {
    if (!is.list(y_coord)) {
      y_coord_lst <- list(y_coord)
    } else{
      y_coord_lst <- y_coord
    }
    if (!is.list(x_coord)) {
      x_coord_lst <- list(x_coord)
    } else{
      x_coord_lst <- x_coord
    }
    if (length(y_coord_lst) != length(x_coord_lst)) {
      stop(error_sound(),
           "x_coord and y_coord must have the same length")
    }
    traj_lst <- vector(mode = "list", length = length(y_coord_lst))
    for (i in seq_along(y_coord_lst)) {
      if (length(y_coord_lst[[i]]) != length(x_coord_lst[[i]])) {
        stop(error_sound(),
             "x_coord and y_coord must have the same length")
      }
      traj_lst[[i]] <-
        data.frame(pos_x = x_coord_lst[[i]], pos_y = y_coord_lst[[i]])
    }
  }
  if (map(traj_lst, ~ any(is.na(.x))) %>% unlist() %>% any()) {
    stop(error_sound(), "There were NA's in the trajectory.")
  }




  if (length(traj_lst) == 1) {
    traj <- traj_lst[[1]]
    p <- ggplot(traj, aes(x = .data$pos_x, y = .data$pos_y)) +
      geom_point(aes(colour = seq_len(nrow(traj))),
                 size = min(4, 100 / nrow(traj))) +
      geom_path(aes(colour = seq_len(nrow(traj)))) +
      geom_point(
        data = dplyr::slice(traj, 1, nrow(traj)),
        #mark first and last point of trajectory
        aes(x = .data$pos_x, y = .data$pos_y),
        size = max(1.5,min(2*min(4, 100 / nrow(traj)),4)),
        color = "red"
      ) +
      scale_colour_gradientn(colours = inferno(1000, end = 0.9)) +
      theme_bw() +
      xlab("X-position") +
      ylab("Y-position") +
      labs(color = 'Step number') +
      theme(
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black")
      )
  } else{
    traj_lengths <- map(traj_lst,  ~ nrow(.x)) %>% unlist()
    track_id <- rep(seq_along(traj_lst), traj_lengths)
    traj <- do.call(rbind, traj_lst)
    traj$id <- as.factor(track_id)

    p <-
      ggplot(traj, aes(
        x = .data$pos_x,
        y = .data$pos_y,
        colour = .data$id
      )) +
      geom_point(size = min(4, 100 / nrow(traj)), alpha = 0.7) +
      geom_path(alpha = 0.7) +
      scale_colour_discrete(type = inferno(length(traj_lst), end = 0.7)) +
      theme_bw() +
      xlab("X-position") +
      ylab("Y-position") +
      labs(color = 'Track ID') +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size =
                                                         1))) +
      theme(
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black")
      )
  }
  return(suppressWarnings(cowplot::ggdraw(p)))
}




#' Circular Scatterplot of Turn Angles and Step Lengths
#'
#' This function produces a circular scatterplot with the step lengths plotted
#' as distance from the center of a circle and the turn angles as angles
#' (polar coordinates).
#'
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{traj_sim}()}. It must contain
#'  the columns \code{traj$angle} and \code{traj$steplength}.
#' @param theta (alternatively) \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable) or "circular" component of pseudo-observations.
#' @param x (alternatively) \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable) or "linear" component of pseudo-observations.
#' @param margin \link[base]{logical} determining whether the marginal kernel
#' density estimates are computed and plotted. Alternatively, \code{margins} can
#' be a \code{'density.circular'} object (as returned e.g. by \code{\link{fit_circ_np}(theta))}.
#'
#' @details You can either specify \code{traj} or the angels and step lengths
#' \code{theta} and \code{x}.
#'
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object.
#'
#' @examples set.seed(123)
#'
#' traj <- traj_sim(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = list(name="vonmises",coef=list(0, 1)),
#'   marginal_lin = list(name="weibull", coef=list(shape=3))
#' )
#' plot1 <- plot_joint_circ(traj, margin = T)
#'
#' @references \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_cop_scat}()}, \code{\link{plot_track}()},
#' \code{\link{plot_cop_surf}()}, \code{\link{plot_joint_scat}()}.
#'
#' @export
#'
plot_joint_circ <- function(traj = NULL,
                      theta = NULL,
                      x = NULL,
                      margin = FALSE) {
  #validate input
  tryCatch({
    check_arg_all(list(
      check_argument_type(traj,
                          type = c("data.frame", "list")),
      check_argument_type(traj,
                          type = "NULL")
    )
    , 2)
    check_arg_all(list(
      check_argument_type(theta,
                          type = "numeric"),
      check_argument_type(theta,
                          type = "NULL")
    )
    , 2)
    check_arg_all(list(
      check_argument_type(x,
                          type = "numeric"),
      check_argument_type(x,
                          type = "NULL")
    )
    , 2)
    check_arg_all(list(
      check_argument_type(margin,
                          type = "logical"),
      check_argument_type(margin,
                          type = "density.circular")
    )
    , 2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })

  #validate input

  if ((is.null(theta) &&
       !is.null(x)) || (!is.null(theta) && is.null(x))) {
    stop(error_sound(), "Specify angles AND step lengths")
  }
  if (is.null(traj) && (is.null(theta) || is.null(x))) {
    stop(error_sound(),
         "Specify either the trajectory or angles and steplengths")
  }
  if (!is.null(traj) && (!is.null(theta) || !is.null(x))) {
    stop(error_sound(),
         "Specify either the trajectory or angles and steplengths")
  }


  if (!is.null(theta)) {
    traj <- data.frame(angle = theta, steplength = x)
  }

  if (any(is.na(traj$angle)) || any(is.na(traj$steplength))) {
    traj <-
      traj[-c(union(which(is.na(traj$angle)), which(is.na(
        traj$steplength
      )))),]
  }


  if (!any(is(margin) == "density.circular")) {
    calc_margin <- margin
    plot_margin <- margin
  } else{
    calc_margin <- F
    plot_margin <- T
    marginal_angle_dens <- margin
  }
  if(calc_margin){
    #rough bandwidth estimate with trigonometric moments. This can fail when
    #the function has no root (i.e. uniroot gives an error) in this case use
    #only a subset of the sample (for speed) and maximize the cross validation
    #log–likelihood with respect to the bandwidth
    bw <-  tryCatch({
      suppressWarnings(opt_circ_bw(traj$angle, method = "nrd", kappa.est = "trigmoments"))
    },
    error = function(e) {
      if (length(traj$angle) > 1000) {
        sample_angle <-
          traj$angle[sample(seq_along(traj$angle), 1000, replace = FALSE)]
      } else{
        sample_angle <- traj$angle
      }
      suppressWarnings(opt_circ_bw(sample_angle, method = "cv"))
    })

    marginal_angle_dens <- traj$angle %>%
      half2full_circ() %>%
      circular::circular(zero = 0, rotation = "counter") %>%
      circular::density.circular(
        #circular marginal density estimated with von Mises kernel density and rough bandwidth estimate
        bw = bw,
        kernel = "vonmises",
        na.rm = TRUE,
        control.circular = list(rotation = "counter", zero =
                                  0)
      )
  }

  marginal_angle_dens$x <- full2half_circ(marginal_angle_dens$x)
  marginal_angle_dens$y <- as.double(marginal_angle_dens$y)



  #set the breaks for the gridlines of the plot
  y_breaks <-
    pretty(c(0, max(traj$steplength, na.rm = TRUE)), n = 6)
  y_breaks <- c(y_breaks, tail(y_breaks, 1) + y_breaks[2])
  x_breaks <- seq(-0.75 * pi, pi, 0.25 * pi)


  if(!plot_margin){
    marginal_angle_dens$x <- seq(-pi,0.75*pi,length.out=7)
    marginal_angle_dens$y <- rep(0,7)
  }

  #blow up the marginal density, so its zero is the the outermost gridline of the circular plot
  marginal_angle_dens <- cbind(x = marginal_angle_dens$x,
                               y = (0.6 * marginal_angle_dens$y + 1) * tail(y_breaks, 1)) %>%
    as.data.frame()



  #set positions of radius labels
  radius_label_pos <-
    cbind(x = rep(c(-1, 1), length(y_breaks))[2:(length(y_breaks) - 1)] *
            0.875 * pi,
          y = y_breaks[2:(length(y_breaks) - 1)]) %>%
    as.data.frame()

  #set positions of angle labels
  max_rad <- max(y_breaks)

  if(plot_margin){
    ang_offset <- c(1.3, 1.3, 1.3, 1.15, 1.3, 1.3, 1.3, 1.15)
  }else{
    ang_offset <- c(1, 1, 1, 0.95, 1, 1, 1, 0.95)
  }

  angle_label_pos <-  cbind(
    x = x_breaks,
    y = ang_offset *
      marginal_angle_dens$y[map(x_breaks,  ~ which.min(abs(marginal_angle_dens$x -
                                                             .x))) %>% unlist()]
  ) %>%
    as.data.frame()

  #convert to cartesian coordinates to determine plot-margins
  angle_label_pos_cart <- angle_label_pos
  angle_label_pos_cart$x <-
    angle_label_pos$y * sin(angle_label_pos$x)
  angle_label_pos_cart$y <-
    angle_label_pos$y * cos(angle_label_pos$x)

  panel_size <-
    max(max(abs(angle_label_pos_cart$x)), max(abs(angle_label_pos_cart$y)))
  top_correction <- panel_size - abs(max(angle_label_pos_cart$y))
  right_correction <- panel_size - abs(max(angle_label_pos_cart$x))
  bottom_correction <- panel_size - abs(min(angle_label_pos_cart$y))
  left_correction <- panel_size - abs(min(angle_label_pos_cart$x))


  #set theme
  if(plot_margin){
    y_intercept_lim <- length(y_breaks)
  }else{
    y_intercept_lim <- length(y_breaks)-1
  }
  circ_plot_layers <- list(
    geom_hline(
      yintercept = y_breaks[seq(2, y_intercept_lim, 2)],
      colour = "darkgreen",
      size = 0.2
    ),
    geom_hline(
      yintercept = y_breaks[seq(1, y_intercept_lim, 2)],
      colour = "darkred",
      size = 0.2
    ),
    coord_polar(start = pi, clip = "off"),

    theme_bw(),

    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
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
        y = if(plot_margin){tail(y_breaks, 1)}else{tail(y_breaks, 2)[1]},
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
      data = radius_label_pos[seq(1, nrow(radius_label_pos), 2),],
      aes(.data$x, .data$y, label = .data$y),
      size = 10 / 2.834646,
      color = "darkgreen",
      label.size = 0
    ),
    geom_label(
      data = radius_label_pos[seq(2, nrow(radius_label_pos), 2),],
      aes(.data$x, .data$y, label = .data$y),
      size = 10 / 2.834646,
      color = "darkred",
      label.size = 0
    ),

    geom_label(
      data = angle_label_pos,
      aes(.data$x, .data$y),
      label = c(
        expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
                   0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
      ),
      size = 10 / 2.834646,
      label.size = 0
    )
  )


  # plot
  if(plot_margin){
    margin_layer <-
      list(
        geom_line(
          data = marginal_angle_dens,
          aes(y = .data$y, x = .data$x),
          colour = "black",
          size = 0.2
        ),
        geom_ribbon(
          data = marginal_angle_dens,
          aes(
            x = .data$x,
            ymin = tail(y_breaks, 1),
            ymax = .data$y
          ),
          alpha = 0.3,
          fill = "blue"
        ),
        geom_hline(
          yintercept = tail(y_breaks, 1),
          colour = "grey30",
          size = 0.2
        )
      )
  }else{
    margin_layer <- NULL
  }


  p <- ggplot() +
    circ_plot_layers +
    geom_point(
      data = traj,
      aes(y = .data$steplength, x = as.numeric(.data$angle)),
      alpha = 0.5,
      size = 0.5
    ) +
    margin_layer
  return(suppressWarnings(cowplot::ggdraw(p)))

}


#' Scatterplot of Copula Values
#'
#' This function produces a scatterplot ('\code{\link[ggplot2]{ggplot}}' object) of
#' a sample from a copula.
#'
#' @param traj a \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{traj_sim}()}, which must contain the columns
#' \code{traj$cop_u} and \code{traj$cop_v}.
#' @param u (alternatively) \link[base]{numeric} \link[base]{vector} of first
#' components of pseudo-observations or draws from a copula.
#' @param v (alternatively) \link[base]{numeric} \link[base]{vector} of second
#' components of pseudo-observations or draws from a copula.
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object, the scatterplot.
#'
#'
#' @examples set.seed(123)
#' traj <- traj_sim(100,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_lin = list(name = "weibull", coef = list(shape = 3))
#' )
#' plot_cop_scat(traj = traj)
#'
#' sample <- rcylcop(100,cyl_quadsec(0.1))
#' plot_cop_scat(u = sample[,1], v = sample[,2])
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_track}()},
#' \code{\link{plot_joint_circ}()}, \code{\link{plot_cop_surf}()}, \code{\link{plot_joint_scat}()}.
#'
#' @export
#'
plot_cop_scat <- function(traj = NULL,
                          u = NULL,
                          v = NULL) {
  tryCatch({
    check_arg_all(list(
      check_argument_type(traj,
                          type = c("data.frame", "list")),
      check_argument_type(traj,
                          type = "NULL")
    )
    , 2)
    check_arg_all(list(
      check_argument_type(u,
                          type = "numeric",
                          lower = 0,
                          upper =1),
      check_argument_type(u,
                          type = "NULL")
    )
    , 2)
    check_arg_all(list(
      check_argument_type(v,
                          type = "numeric",
                          lower = 0,
                          upper =1),
      check_argument_type(v,
                          type = "NULL")
    )
    , 2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })

  if ((is.null(u) &&
       !is.null(v)) || (!is.null(u) && is.null(v))) {
    stop(error_sound(), "Specify pseudo-observations u AND v!")
  }
  if (is.null(traj) && (is.null(u) || is.null(v))) {
    stop(error_sound(),
         "Specify either the trajectory or pseudo-observations!")
  }
  if (!is.null(traj) && (!is.null(u) || !is.null(v))) {
    stop(error_sound(),
         "Specify either the trajectory or pseudo-observations!")
  }
  if (length(u) != length(v)) {
    stop(error_sound(), "u and v must have the same length!")
  }


  if (!is.null(u)) {
    traj <- data.frame(cop_u = u, cop_v = v)
  }

  if (!all(c("cop_u", "cop_v") %in% colnames(traj))) {
    stop(error_sound(),
         "trajectory must contain the columns 'cop_u' and 'cop_v'")
  }

  if (any(is.na(traj$cop_u)) || any(is.na(traj$cop_v))) {
    traj <-
      traj[-c(union(which(is.na(traj$cop_u)), which(is.na(
        traj$cop_v
      )))),]
  }

  plot_theme <- list(
    geom_point(size = 0.5, alpha = 0.5),
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
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    plot_theme +
    coord_fixed()
  return(suppressWarnings(cowplot::ggdraw(p)))
}



#' Surface Plot or Heat Map of the Distribution or the Density of a Copula
#'
#' This function plots the distribution or the density of a copula. It can produce
#' a 3 dimensional surface plot using either functions from the '\pkg{rgl}' or from the
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
#' if(interactive()){
#'  plot_cop_surf(copula::frankCopula(2),
#'    type="pdf",
#'    plot_type="ggplot",
#'    resolution = 5
#'  )
#'  plot_cop_surf(copula::frankCopula(2),
#'    type="cdf",
#'    plot_type="ggplot",
#'    resolution = 5
#'  )
#'
#' #opens a new window
#'   plot_cop_surf(cyl_quadsec(0.1),
#'     type="pdf",
#'     plot_type="rgl"
#'   )
#'   plot_cop_surf(cyl_quadsec(0.1),
#'     type="pdf",
#'     plot_type="rgl",
#'     n_gridlines = 60
#'   )
#'
#'   plot_cop_surf(cyl_quadsec(0.1),
#'     type="pdf",
#'     plot_type="plotly",
#'     n_gridlines = 10,
#'     resolution = 10
#'   )
#' }
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_cop_scat}()}, \code{\link{plot_track}()},
#' \code{\link{plot_joint_circ}()}, \code{\link{plot_joint_scat}()}.
#'
#' @export
#'
plot_cop_surf <- function(copula,
                     type = "pdf",
                     plot_type = "rgl",
                     resolution = 51,
                     n_gridlines = 11) {
  #validate input
  tryCatch({
    check_arg_all(list(
      check_argument_type(copula,
                          type = "Copula"),
      check_argument_type(copula,
                          type = "cyl_copula")
    )
    , 2)
    check_arg_all(check_argument_type(
      argument = type,
      type = "character",
      values = c("pdf", "cdf"),
      length = 1
    )
    ,
    1)
    check_arg_all(check_argument_type(
      plot_type,
      type = "character",
      values = c("rgl", "plotly", "ggplot"),
      length = 1
    )
    ,
    1)
    check_arg_all(
      check_argument_type(
        resolution,
        type = "numeric",
        integer = T,
        length = 1,
        lower = 1
      )
      ,
      1
    )
    check_arg_all(
      check_argument_type(
        n_gridlines,
        type = "numeric",
        integer = T,
        length = 1,
        lower = 0
      )
      ,
      1
    )
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })

  #check input

  if (type == "pdf")
    fun <- dcylcop
  else if (type == "cdf")
    fun <- pcylcop
  else
    stop(error_sound(), "Type must be either cdf or pdf")

  if (cylcop.env$silent == F)
    printCop(copula)


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

    outp <- rep(NA, resolution ^ 2)
    outp[seq_len(resolution)] <-
      fun(mat[seq_len(resolution),], copula)
    outp[(nrow(mat) - resolution + 1):nrow(mat)] <-
      fun(mat[(nrow(mat) - resolution + 1):nrow(mat),], copula)


    if (cylcop.env$silent == F) {
      #generate progressbar
      time <-
        (proc.time() - ptm)[3] * resolution / 2 %>% as.integer()
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
      outp[((i - 1) * resolution + 1):(i * resolution)] <-
        fun(mat[((i - 1) * resolution + 1):(i * resolution),], copula)
      if (time > 10 && cylcop.env$silent == F) {
        utils::setTxtProgressBar(pb, i)
        if (i == resolution / 2)
          waiting_sound()
      }
    }
    if (time > 10 && cylcop.env$silent == F) {
      close(pb)
      done_sound()
    }

    quant_995 <- quantile(outp, probs = 0.995)
    max_outp <- max(outp)
    if (max_outp > (2 * quant_995)) {
      outp[which(outp > quant_995)] <- quant_995
      warning(
        warning_sound(),
        paste0(
          "Maximum of ",
          type,
          " is ",
          round(max_outp, 2),
          "; for better visulatization,\n values larger than ",
          round(quant_995, 2),
          " (99.5 percentile) are cut off."
        )
      )
    }

    mat <-
      matrix(outp,
             ncol = resolution,
             nrow = resolution,
             byrow = F)



    #if gridlines are at already calculated gridvalues, no need to recalculate them
    if ((resolution - 1) %% (n_gridlines - 1) == 0) {
      if (n_gridlines == 0)
        mat_grid <- u_grid
      else{
        mat_grid <-
          mat[seq(1, resolution, (resolution - 1) / (n_gridlines - 1)),
              seq(1, resolution, (resolution - 1) / (n_gridlines - 1))]
      }
    }

    #calculate gridlines
    else{
      time <- time / (resolution ^ 2)

      #Calculate timing
      if (time * (n_gridlines ^ 2) > 10 && cylcop.env$silent == F) {
        message(
          "gridlines are not a multiple of the surface points and have to be newly calculated.\nThis will take ",
          floor(time * (n_gridlines ^ 2) / 60),
          " minutes, ",
          round((time * (
            n_gridlines ^ 2
          )) - 60 * floor(time * (
            n_gridlines ^ 2
          ) / 60)),
          " seconds"
        )
      }

      #calculate values for gridlines
      mat_grid <- as.matrix(expand.grid(u = u_grid, v = v_grid))
      quant_995 <- quantile(outp, probs = 0.995)
      outp_grid <- fun(mat_grid, copula)
      if (max_outp > (2 * quant_995)) {
        outp_grid[which(outp_grid > quant_995)] <- quant_995
      }
      mat_grid <-
        matrix(outp_grid,
               ncol = n_gridlines,
               nrow = n_gridlines,
               byrow = F)
    }

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
            paste(
              parameter_summary,
              ", ",
              copula@param.names[i],
              " = ",
              round(copula@parameters[i], 3)
            )
        }
      }
      #rgl doesn't like \n character, so need to split title in 2
      if("cyl_rect_combine" %in%is(copula)){
        base_cop <- get_cop_name(copula@cop.lo)
        title1 <- paste(type, "of rectangular patchwork of")
        title2 <- base_cop
      }else{
        title1 <-
          paste(type, "of", copula@name, "with")
        title2 <- parameter_summary
      }
    }


    #------------Make surface plot with rgl------

    if (plot_type == "rgl") {
      # Generate title

      # Make plot

      p <-  rgl::persp3d(
        v,
        u,
        t(mat),
        col = inferno(50)[cut(mat, breaks = 50)],
        fit = "smooth",
        lit = TRUE,
        alpha = 0.9,
        xlab = "v",
        ylab = "u",
        zlab = if (type == "pdf")
          "c(u,v)"
        else
          "C(u,v)",
        expand = 0
      ) +      rgl::light3d(theta = 0, phi = 30) +
        rgl::title3d(main = paste(title1,title2)) +
        (if (n_gridlines > 0) {
          rgl::surface3d(
            v_grid,
            u_grid,
            t(mat_grid),
            color = "black",
            lit = FALSE,
            front = "lines",
            back = "lines"
          )
        })
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

      gridlines_y <- expand.grid(v = u_grid, u = v_grid)
      gridlines_y <- gridlines_y[c("u", "v")]
      zg <- c(t(mat_grid))
      gridlines_y <-
        mutate(gridlines_y, zg) %>% dplyr::rename(ug = u) %>% dplyr::rename(vg = v)


      # Make plot

      p <-
        plot_ly(
          x = v,
          y = u,
          z = mat,
          showscale = FALSE,
          colorscale = col,
          opacity = 1,
          lighting = list(ambient = 0.9, specular = 0)
        ) %>%
        plotly::layout(
          title = paste0(title1,"\n",title2),
          scene = list(
            xaxis = list(title = "v"),
            yaxis = list(title = "u"),
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
        ) %>% add_surface(
          hovertemplate = paste0('u: %{x:.3f}<br>',

                                'v: %{y:.3f}<br>',

                                if (type == "pdf"){
                                  'c(u,v):  %{z:.3f}'
                                  }else{
                                  'C(u,v):  %{z:.3f}'
                                })
        )


      # Add gridlines
      if (n_gridlines > 0) {

        gridline_fun <- function(plot, data, offset=0.1){
          #for aesthetics, add a small number so the gridlines are slightly above the surface
          data$zg <- data$zg+offset
          add_trace(
            plot,
            data = data,
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
          p <- gridline_fun(p, gridlines_x[i:(i + n_gridlines - 1), ], offset = 0.02)
          p <- gridline_fun(p, gridlines_y[i:(i + n_gridlines - 1), ], offset = 0.02)
          p <- gridline_fun(p, gridlines_x[i:(i + n_gridlines - 1), ], offset = -0.02)
          p <- gridline_fun(p, gridlines_y[i:(i + n_gridlines - 1), ], offset = -0.02)
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
    outp$z <-  fun(as.matrix(outp[, 1:2]), copula)
    quant_995 <- quantile(outp$z, probs = 0.995)
    max_outp <- max(outp$z)
    if (max_outp > (2 * quant_995)) {
      outp$z[which(outp$z > quant_995)] <- quant_995
      warning(
        warning_sound(),
        paste0(
          "Maximum of ",
          type,
          " is ",
          round(max_outp, 2),
          "; for better visulatization,\n values larger than ",
          round(quant_995, 2),
          " (99.5 percentile) are cut off."
        )
      )
    }

    # Make plot

    gridline_size <- if(n_gridlines <= 11){0.8}else if(n_gridlines>=11 && n_gridlines <21){0.6}else{0.3}

    p <- ggplot(outp, aes(v, u)) +
      geom_raster(
        aes(fill = .data$z),
        interpolate = T,
        hjust = 0,
        vjust = 0
      ) +
      coord_fixed() +
      theme_bw() +
      labs(fill = if (type == "pdf")
        "c(u,v)"
        else
          "C(u,v)") +
      scale_fill_gradientn(colours = inferno(1000), limits = c(0, max(outp$z))) +
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
        xintercept = seq(0, 1, length.out=n_gridlines),
        colour = "grey60",
        size = gridline_size,
        alpha = 0.25
      ) +
      geom_hline(
        yintercept = seq(0, 1, length.out=n_gridlines),
        colour = "grey60",
        size = gridline_size,
        alpha = 0.25
      ) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
    suppressWarnings(plot(p))
  }
}


#' Circular Boxplot of Turn Angles and Step Lengths
#'
#' This function produces circular boxplots (a '\code{\link[ggplot2]{ggplot}}' object)
#' of the turn angles corresponding to specific quantiles of the step lengths.
#'
#' The step lengths are split into quantiles. For each quantile a boxplot of the
#' corresponding turn angles is produced and wrapped around the circle.
#' The turn angle values are plotted
#' as scatter plot overlaying the boxplot.  Outliers are plotted in red.
#'
#' @param traj \link[base]{data.frame} containing the trajectory produced by e.g.
#' \code{\link{traj_sim}()}. It must contain
#'  the columns \code{traj$angle} and \code{traj$steplength}.
#' @param theta (alternatively) \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable) or "circular" component of pseudo-observations.
#' @param x (alternatively) \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable) or "linear" component of pseudo-observations.
#' @param levels \link[base]{integer} value between 1 and 15, the number of
#' quantiles into which the step lengths are split.
#' @param marginal_lin named \link[base]{list} (for parametric estimates) or
#' a '\code{\link[stats]{density}}' object (for kernel density estimates).
#' The outputs of functions \code{\link{fit_lin_np}()} or
#' \code{\link{fit_lin_param}()} can be used here directly.
#'  If \code{marginal_lin} is specified, the limits of the quantiles of the step lengths
#'  are determined from that distribution instead of from the data specified with
#' \code{traj$steplength} or \code{x}.
#' @param spacing \link[base]{numeric} value between 0 and 10 determining the
#' spacing between the boxplots.
#' @param legend_pos \link[base]{character} string denoting the position of the legend (limits
#' of the step length quantiles). Either \code{"left"}, \code{"right"}, \code{"top"}, or
#' \code{"bottom"}
#'
#' @details
#' The distance of the plotted points from the center
#'  of the circle is of no meaning and is randomized.
#' Only their angular position is relevant.
#' The median of the turn angles is defined as the center of the shortest arc
#' that connects all points. The length of the whiskers is 1.5 times the interquartile range.
#'
#' You can either specify \code{traj} or the angels (\code{theta})
#' and linear masurements (code{x}).
#' If entered "by hand", the named list describing the marginal linear distribution
#' (for \code{marginal_lin}) must contain 2 entries:
#' \enumerate{
#'    \item{\code{name}:
#' a \link[base]{character} string denoting the name of the linear distribution,
#' i.e. the name of its
#'    distribution function without the "p",
#'   e.g. "norm" for normal distribution.}
#'    \item{\code{coef}: a named list containing the parameters of the distribution
#'    given in \code{"name"}.}
#' }
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object, the circular boxplot.
#'
#' @examples set.seed(1234)
#'
#' traj <- traj_sim(100,
#'   copula = cyl_rect_combine(copula::frankCopula(6)),
#'   marginal_circ = list(name= "vonmises", coef=list(0, 2)),
#'   marginal_lin = list(name = "weibull", coef=list(shape=3))
#' )
#'
#' plot1 <- plot_joint_box(traj)
#' plot2 <- plot_joint_box(traj,
#'   marginal_lin=list(name = "weibull", coef=list(shape=3))
#')
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_cop_scat}()}, \code{\link{plot_track}()},
#' \code{\link{plot_joint_circ}()}, \code{\link{plot_cop_surf}()}.
#'
#' @export
#'
plot_joint_box <-
  function(traj = NULL,
           theta = NULL,
           x = NULL,
           levels = 5,
           marginal_lin = NULL,
           spacing = 0.3,
           legend_pos = "right",
           min_a=3,
           max_a=4353) {
    #validate input
    tryCatch({
      check_arg_all(list(
        check_argument_type(traj,
                            type = c("data.frame", "list")),
        check_argument_type(traj,
                            type = "NULL")
      )
      , 2)
      check_arg_all(list(
        check_argument_type(theta,
                            type = "numeric"),
        check_argument_type(theta,
                            type = "NULL")
      )
      , 2)
      check_arg_all(list(
        check_argument_type(x,
                            type = "numeric"),
        check_argument_type(x,
                            type = "NULL")
      )
      , 2)
      check_arg_all(
        check_argument_type(
          levels,
          type = "numeric",
          length = 1,
          integer = T,
          lower = 1,
          upper = 15
        )
        ,
        1
      )
      check_arg_all(list(
        check_argument_type(marginal_lin,
                            type = "list"),
        check_argument_type(marginal_lin,
                            type = "density"),
        check_argument_type(marginal_lin,
                            type = "NULL")
      )
      ,
      3)
      check_arg_all(check_argument_type(
        spacing,
        type = "numeric",
        length = 1,
        lower = 0,
        upper = 10
      )
      ,
      1)
      check_arg_all(check_argument_type(
        legend_pos,
        type = "character",
        length = 1,
        values = c("left", "right", "top", "bottom")
      )
      ,
      1)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    })

    expand <- 1
    # IQR_fac: positive value determining the lengths of the whiskers of
    # the boxplots as factor of their interquartile ranges. All points falling outside the whiskers
    # are considered outliers.
    IQR_fac <- 1.5
    levels <- as.integer(levels)

    if (IQR_fac < 0) {
      stop("The length of the boxplot whiskers must be larger than 0")
    }

    if ((is.null(theta) &&
         !is.null(x)) || (!is.null(theta) && is.null(x))) {
      stop(error_sound(), "Specify angles AND step lengths")
    }
    if (is.null(traj) && (is.null(theta) || is.null(x))) {
      stop(error_sound(),
           "Specify either the trajectory or angles and steplengths")
    }
    if (!is.null(traj) && (!is.null(theta) || !is.null(x))) {
      stop(error_sound(),
           "Specify either the trajectory or angles and steplengths")
    }
    if (!is.null(marginal_lin)) {
      if ("density" %in% is(marginal_lin)) {
        parameter_lin <- marginal_lin
        marginal_lin <- "dens"
      } else{
        if (!all(c("name", "coef") %in% names(marginal_lin))) {
          stop(
            error_sound(),
            "marginal_lin must be a density object or list containing
             the entries 'name' and 'coef'."
          )
        }
        if (!is.character(marginal_lin$name)) {
          stop(error_sound(),
               "In marginal_lin: name must be of type character.")
        }
        if (!is.list(marginal_lin$coef)) {
          stop(error_sound(),
               "In marginal_lin: coef must be of type list.")
        }
        parameter_lin <- marginal_lin$coef
        marginal_lin <- marginal_lin$name
      }
    } else{
      parameter_lin <- NULL
      marginal_lin <- NULL
    }

    if (!is.null(theta)) {
      traj <- data.frame(angle = theta, steplength = x)
    }

    if (any(is.na(traj$angle)) || any(is.na(traj$steplength))) {
      traj <-
        traj[-c(union(which(is.na(traj$angle)), which(is.na(
          traj$steplength
        )))),]
    }

    plot_data <-
      calc_plot_data(
        traj = traj,
        marginal_lin = marginal_lin,
        parameter_lin = parameter_lin,
        IQR_fac = IQR_fac,
        spacing = spacing,
        levels = levels,
        expand = expand,
        min_a = min_a,
        max_a = max_a
      )

    make_plot(
      box = plot_data$box,
      box_border = plot_data$box_border,
      medians = plot_data$medians,
      whiskers = plot_data$whiskers,
      endbars = plot_data$endbars,
      outliers = plot_data$outliers,
      dat_jitter = plot_data$dat_jitter,
      quant = plot_data$quant,
      spacing = spacing,
      levels = levels,
      legend_pos = legend_pos,
      expand = expand
    )
  }

calc_plot_data <-
  function(traj,
           marginal_lin,
           parameter_lin,
           IQR_fac,
           spacing,
           levels,
           expand,
           min_a,
           max_a) {
    box_heights <-
      cbind(seq(spacing, (levels - 1 + (spacing * levels)), length.out = levels),
            seq((1 + spacing), (levels + (spacing * levels)), length.out = levels))

    quant <-
      quantile(traj$steplength, seq(0, 1, length.out = (levels + 1))[-1])

    if (!is.null(marginal_lin)) {
      r <- get_marg(marginal_lin)$r
      marginal_dist_sample <-
        do.call(r, c(n = 100000, parameter_lin))
      quant <-
        quantile(marginal_dist_sample, seq(0, 1, length.out = (levels + 1))[-1])
    }

    traj$quant <- 1
    for (i in 2:levels) {
      traj$quant[traj$steplength > quant[(i - 1)]] <- i
    }
    traj$quant <- as.factor(traj$quant)
    quant_ang <- matrix(ncol = 3, nrow = levels)
    for (i in seq_len(levels)) {
      angle <- half2full_circ(traj$angle[traj$quant == i])
      quant_ang[i, 2] <-
        circular::median.circular(circular::circular(angle, modulo = "2pi", zero = 0))
      opposite_median <- (quant_ang[i, 2] + pi) %% (2 * pi)
      if (opposite_median < quant_ang[i, 2]) {
        angle_above_median1 <- which(angle > quant_ang[i, 2])
        angle_above_median2 <- which(angle < opposite_median)
        angle[angle_above_median2] <-
          angle[angle_above_median2] + (2 * pi)
        angle_above_median <-
          c(angle_above_median1, angle_above_median2)
        quant_ang[i, 3] <- median(angle[angle_above_median])
        angle[angle_above_median2] <-
          angle[angle_above_median2] - (2 * pi)
        if (quant_ang[i, 3] > (2 * pi))
          quant_ang[i, 3] <- quant_ang[i, 3] - (2 * pi)
        quant_ang[i, 1] <- median(angle[-angle_above_median])

      }
      if (opposite_median >= quant_ang[i, 2]) {
        angle_below_median1 <- which(angle < quant_ang[i, 2])
        angle_below_median2 <- which(angle > opposite_median)
        angle[angle_below_median2] <-
          angle[angle_below_median2] - (2 * pi)
        angle_below_median <-
          c(angle_below_median1, angle_below_median2)
        quant_ang[i, 1] <- median(angle[angle_below_median])
        angle[angle_below_median2] <-
          angle[angle_below_median2] + (2 * pi)
        if (quant_ang[i, 1] < 0)
          quant_ang[i, 1] <- quant_ang[i, 1] + (2 * pi)
        quant_ang[i, 3] <- median(angle[-angle_below_median])
      }
      quant_ang[i,] <- full2half_circ(quant_ang[i,])
      if (quant_ang[i, 1] > quant_ang[i, 3]) {
        temp <- quant_ang[i, 1]
        quant_ang[i, 1] <- quant_ang[i, 3]
        quant_ang[i, 3] <- temp
      }
    }

    box <- vector(mode = "list", length = levels)
    box_border_x <- vector(mode = "list", length = levels)
    box_border_y <- vector(mode = "list", length = levels)
    medians <- vector(mode = "list", length = levels)
    whiskers <- vector(mode = "list", length = levels)
    endbars <- vector(mode = "list", length = levels)
    outliers <- vector(mode = "list", length = levels)
    dat_jitter <- vector(mode = "list", length = levels)
    clockwise <- rep(T, levels)
    for (i in seq_len(levels)) {
      if (quant_ang[i, 1] < 0 &&
          quant_ang[i, 2] > 0 &&
          quant_ang[i, 2] > quant_ang[i, 3]) {
        box[[i]] <-
          rbind(
            cbind(-pi, quant_ang[i, 1], box_heights[i, 1], box_heights[i, 2]),
            cbind(pi, quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2]),
            cbind(quant_ang[i, 3], quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2])
          )
        clockwise[i] <- F
      }
      else if (quant_ang[i, 2] < 0 &&
               quant_ang[i, 3] > 0 &&
               quant_ang[i, 1] >= quant_ang[i, 2]) {
        box[[i]] <-
          rbind(
            cbind(quant_ang[i, 1], quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2]),
            cbind(quant_ang[i, 2],-pi, box_heights[i, 1], box_heights[i, 2]),
            cbind(pi, quant_ang[i, 3], box_heights[i, 1], box_heights[i, 2])
          )
        clockwise[i] <- F
      }
      else if (quant_ang[i, 3] < 0 &&
               quant_ang[i, 2] > 0 &&
               quant_ang[i, 1] >= quant_ang[i, 2]) {
        box[[i]] <-
          rbind(
            cbind(quant_ang[i, 1], quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2]),
            cbind(quant_ang[i, 2], pi, box_heights[i, 1], box_heights[i, 2]),
            cbind(-pi, quant_ang[i, 3], box_heights[i, 1], box_heights[i, 2])
          )
        clockwise[i] <- T
      }
      else if (quant_ang[i, 1] > 0 &&
               quant_ang[i, 2] < 0 &&
               quant_ang[i, 2] < quant_ang[i, 3]) {
        box[[i]] <-
          rbind(
            cbind(pi, quant_ang[i, 1], box_heights[i, 1], box_heights[i, 2]),
            cbind(-pi, quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2]),
            cbind(quant_ang[i, 3], quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2])
          )
        clockwise[i] <- T
      }
      else{
        box[[i]] <-
          rbind(
            cbind(quant_ang[i, 1], quant_ang[i, 2], box_heights[i, 1], box_heights[i, 2]),
            cbind(quant_ang[i, 2], quant_ang[i, 3], box_heights[i, 1], box_heights[i, 2])
          )
        if (quant_ang[i, 1] < quant_ang[i, 3]) {
          clockwise[i] <- T
        }
        else{
          clockwise[i] <- F
        }
      }
      box_border_x[[i]] <-
        cbind(rep(box[[i]][, 1], 2),
              rep(box[[i]][, 2], 2),
              c(box[[i]][, 3], box[[i]][, 4]),
              c(box[[i]][, 3], box[[i]][, 4]))
      box_border_y[[i]] <-
        rbind(c(rep(quant_ang[i, 1], 2), box_heights[i, 1], box_heights[i, 2]),
              c(rep(quant_ang[i, 3], 2), box_heights[i, 1], box_heights[i, 2]))
      inter_quart_range <-
        IQR_fac * abs(add_angles(quant_ang[i, 3], (-1 * quant_ang[i, 1])))

      whiskers_start <- c(quant_ang[i, 1], quant_ang[i, 3])
      if (clockwise[i] == T) {
        whiskers_max <-
          add_angles(whiskers_start,
                     c(-inter_quart_range, inter_quart_range))
      }
      else{
        whiskers_max <-
          add_angles(whiskers_start,
                     c(inter_quart_range,-inter_quart_range))
      }
      opposite_median <- add_angles(quant_ang[i, 2], pi)
      if (inter_quart_range >= (2 * pi)) {
        whiskers_max <- rep(opposite_median, 2)
      }
      if (!is.null(
        within_angles(
          angle1 = whiskers_start[2],
          angle2 = whiskers_max[2],
          data = opposite_median,
          clockwise = clockwise[i]
        )
      )) {
        whiskers_max[2] <- opposite_median
      }
      if (!is.null(
        within_angles(
          angle1 = whiskers_start[1],
          angle2 = whiskers_max[1],
          data = opposite_median,
          clockwise = !clockwise[i]
        )
      )) {
        whiskers_max[1] <- opposite_median
      }

      angles <- traj$angle[traj$quant == i]
      ind1 <- within_angles(
        angle1 = whiskers_start[1],
        angle2 = whiskers_max[1],
        data = angles,
        clockwise = !clockwise[i]
      )
      whiskers_max[1] <-
        angles[ind1][which.min(abs(add_angles(whiskers_max[1], (-1 * angles[ind1]))))]
      ind2 <- within_angles(
        angle1 = whiskers_start[2],
        angle2 = whiskers_max[2],
        data = angles,
        clockwise = clockwise[i]
      )
      whiskers_max[2] <-
        angles[ind2][which.min(abs(add_angles(whiskers_max[2], (-1 * angles[ind2]))))]
      ind3 <- within_angles(
        angle1 = quant_ang[i, 1],
        angle2 = quant_ang[i, 3],
        data = angles,
        clockwise = clockwise[i]
      )
      avg_height <- (box_heights[i, 1] + box_heights[i, 2]) / 2

      angles_a <- traj$angle[min_a:max_a][traj$quant[min_a:max_a] == i]
      ind1_a <- within_angles(
        angle1 = whiskers_start[1],
        angle2 = whiskers_max[1],
        data = angles_a,
        clockwise = !clockwise[i]
      )
      ind2_a <- within_angles(
        angle1 = whiskers_start[2],
        angle2 = whiskers_max[2],
        data = angles_a,
        clockwise = clockwise[i]
      )
      ind3_a <- within_angles(
        angle1 = quant_ang[i, 1],
        angle2 = quant_ang[i, 3],
        data = angles_a,
        clockwise = clockwise[i]
      )

      outliers[[i]] <-
        cbind(c(angles[-c(ind1, ind2, ind3)]), rep(avg_height, length(angles[-c(ind1, ind2, ind3)])))
      dat_jitter[[i]] <-
        cbind(c(angles_a[c(ind1_a, ind2_a, ind3_a)]), runif(length(angles_a[c(ind1_a, ind2_a, ind3_a)]), (avg_height -
                                                                                                            0.3), (avg_height + 0.3)), i)
      endbars[[i]] <-
        cbind(c(whiskers_max), c(whiskers_max), rep((avg_height - 0.4), 2), rep((avg_height +
                                                                                   0.4), 2))
      if (whiskers_start[1] < whiskers_max[1] &&
          clockwise[i] == T) {
        whiskers[[i]] <- rbind(
          c(whiskers_start[1],-pi, avg_height, avg_height),
          c(pi, whiskers_max[1], avg_height, avg_height)
        )
      }
      else if (whiskers_start[1] > whiskers_max[1] &&
               clockwise[i] == F) {
        whiskers[[i]] <- rbind(
          c(whiskers_start[1], pi, avg_height, avg_height),
          c(-pi, whiskers_max[1], avg_height, avg_height)
        )
      }
      else{
        whiskers[[i]] <-
          rbind(c(whiskers_start[1], whiskers_max[1], avg_height, avg_height))
      }
      if (whiskers_start[2] > whiskers_max[2] &&
          clockwise[i] == T) {
        whiskers[[i]] <- rbind(whiskers[[i]],
                               rbind(
                                 c(whiskers_start[2], pi, avg_height, avg_height),
                                 c(-pi, whiskers_max[2], avg_height, avg_height)
                               ))
      }
      else if (whiskers_start[2] < whiskers_max[2] &&
               clockwise[i] == F) {
        whiskers[[i]] <- rbind(whiskers[[i]],
                               rbind(
                                 c(whiskers_start[2],-pi, avg_height, avg_height),
                                 c(pi, whiskers_max[2], avg_height, avg_height)
                               ))
      }
      else{
        whiskers[[i]] <- rbind(whiskers[[i]],
                               rbind(
                                 c(whiskers_start[2], whiskers_max[2], avg_height, avg_height)
                               ))
      }

      medians[[i]] <-
        rbind(c(rep(quant_ang[i, 2], 2), box_heights[i, 1], box_heights[i, 2]))
    }
    box <- do.call(rbind, args = box)
    label <- as.factor(c(box[, 4]))
    levels(label) <- seq(1, levels)
    box <-
      data.frame(
        xmin = c(box[, 1]),
        xmax = c(box[, 2]),
        ymin = c(box[, 3]),
        ymax = c(box[, 4]),
        label = label
      )
    box_border_x <- do.call(rbind, args = box_border_x)
    box_border_y <- do.call(rbind, args = box_border_y)
    box_border <- rbind(box_border_x, box_border_y)
    box_border <-
      data.frame(
        x = c(box_border[, 1]),
        xend = c(box_border[, 2]),
        y = c(box_border[, 3]),
        yend = c(box_border[, 4])
      )
    medians <- do.call(rbind, args = medians)
    medians <-
      data.frame(
        x = c(medians[, 1]),
        xend = c(medians[, 2]),
        y = c(medians[, 3]),
        yend = c(medians[, 4])
      )
    whiskers <- do.call(rbind, args = whiskers)
    whiskers <-
      data.frame(
        x = c(whiskers[, 1]),
        xend = c(whiskers[, 2]),
        y = c(whiskers[, 3]),
        yend = c(whiskers[, 4])
      )
    outliers <- do.call(rbind, args = outliers)
    outliers <-
      data.frame(x = c(outliers[, 1]), y = c(outliers[, 2]))
    endbars <- do.call(rbind, args = endbars)
    endbars <-
      data.frame(
        x = c(endbars[, 1]),
        xend = c(endbars[, 2]),
        y = c(endbars[, 3]),
        yend = c(endbars[, 4])
      )
    dat_jitter <- do.call(rbind, args = dat_jitter)
    dat_jitter <-
      data.frame(
        x = c(dat_jitter[, 1]),
        y = c(dat_jitter[, 2]),
        size = (c(dat_jitter[, 3]) / levels * 0.2 * expand)
      )
    out <- list(
      box = box,
      box_border = box_border,
      medians = medians,
      whiskers = whiskers,
      endbars = endbars,
      outliers = outliers,
      dat_jitter = dat_jitter,
      quant = quant
    )
    return(out)
  }

make_plot <-
  function(box,
           box_border,
           medians,
           whiskers,
           endbars,
           outliers,
           dat_jitter,
           quant,
           spacing,
           levels,
           legend_pos,
           expand) {
    box_heights <-
      cbind(seq(spacing, (levels - 1 + (spacing * levels)), length.out = levels),
            seq((1 + spacing), (levels + (spacing * levels)), length.out = levels))
    #set the breaks for the gridlines of the plot
    y_breaks <-
      sort(c(c(box_heights), (max(box_heights) + spacing)))
    x_breaks <- seq(-0.75 * pi, pi, 0.25 * pi)
    #set theme
    circ_plot_layers <- list(
      coord_polar(start = pi, clip = "off"),

      scale_x_continuous(
        limits = c(-pi, pi),
        breaks = x_breaks,
        labels = c(
          expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
                     0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
        )
      ),

      scale_y_continuous(breaks = y_breaks[seq_along(y_breaks)]),
      geom_hline(
        yintercept = y_breaks[(length(y_breaks))],
        colour = "grey60",
        size = 0.3 * expand
      ),

      theme_bw(),
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = expand * 12, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = expand * 12, colour = "black"),
        legend.text = element_text(size = expand * 12, colour = "black"),
        legend.position = legend_pos,
        legend.text.align = 0
      ),
      if (legend_pos == "bottom") {
        guides(fill = guide_legend(nrow = levels, byrow = TRUE))
      },
      # scale_fill_discrete(name = paste0("Step length \n",levels,"-quantile"),
      #                     labels = paste0(seq(1,levels),": (",format(round(c(0.00,quant[1:(levels-1)]), digits=2), nsmall = 2),
      #                                     ",",format(round(c(quant[1:levels]), digits=2), nsmall = 2),"]")),
      scale_fill_viridis(
        name = paste0("Step Length \n", levels, "-Quantile"),
        labels = parse(text = paste0(
          seq(1, levels),
          c('^st', '^nd', '^rd', rep('^th', levels))[seq_len(levels)],
          ' ~ ": (',
          format(round(c(0.0, quant[seq_len(levels -
                                              1)]), digits = 1), nsmall = 1),
          ', ',
          format(round(c(quant[seq_len(levels)]), digits =
                         1), nsmall = 1),
          ']"'
        )),
        # labels = str2expression(paste0(seq(1,levels),c("^st","^nd","^rd",rep("^th",levels))[1:levels],"u")),
        alpha = 1,
        begin = 0.2,
        end = 0.95,
        direction = -1,
        discrete = TRUE,
        option = "inferno"
      ),
      geom_segment(
        aes(
          x = x_breaks,
          y = 0,
          xend = x_breaks,
          yend = y_breaks[(length(y_breaks))]
        ),
        colour = "grey90",
        size = 0.3 * expand
      ),
      #for position of axis labels
      geom_hline(
        yintercept = y_breaks[(length(y_breaks))] + 0.5 * expand,
        colour = "white",
        alpha = 0,
        size = 0.3 * expand
      )
    )


    # plot

    ggplot() +
      circ_plot_layers +
      geom_rect(
        data = box,
        aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data$label
        ),
        alpha = 0.6
      ) +
      geom_segment(data = box_border,
                   aes(
                     x = .data$x,
                     y = .data$y,
                     xend = .data$xend,
                     yend = .data$yend
                   )) +
      geom_segment(
        data = medians,
        aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        ),
        size = 0.5 * expand
      ) +
      geom_segment(
        data = whiskers,
        aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        ),
        size = 0.4 * expand
      ) +
      geom_segment(
        data = endbars,
        aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        ),
        size = 0.4 * expand
      ) +
      geom_point(
        data = outliers,
        aes(x = .data$x,
            y = .data$y),
        color = "red",
        size = 0.6 * expand
      ) +
      geom_point(
        data = dat_jitter,
        aes(
          x = .data$x,
          y = .data$y,
          size = .data$size * 0.001
        ),
        color = "black",
        alpha = 0.3,
        show.legend = FALSE
      ) +
      scale_size_continuous(range = c(0.05, 1))
  }
add_angles <- function(angle1, angle2) {
  angle_sum <- angle1 + angle2
  for (i in seq_along(angle_sum)) {
    if (angle_sum[i] > pi)
      angle_sum[i] <- -pi + (angle_sum[i] - pi)
    else if (angle_sum[i] < -pi)
      angle_sum[i] <- pi - (-pi - angle_sum[i])
  }
  return(angle_sum)
}

within_angles <- function(angle1, angle2, data, clockwise) {
  if (clockwise == TRUE) {
    if (angle2 > angle1) {
      within_angle <- which(data >= angle1 & data <= angle2)
    }
    else if (angle2 < angle1) {
      within_angle <-
        c(within_angles(angle1, pi, data, T),
          within_angles(-pi, angle2, data, T))
    }
  }
  if (clockwise == FALSE) {
    if (angle2 < angle1) {
      within_angle <- which(data >= angle2 & data <= angle1)
    }
    else if (angle2 > angle1) {
      within_angle <-
        c(within_angles(angle1,-pi, data, F),
          within_angles(pi, angle2, data, F))
    }
  }
  if (!any(within_angle))
    within_angle <- NULL
  return(within_angle)
}


#' Circular Histogram of Turn Angles
#'
#' This function produces a circular histogram of turn angles, i.e. angles
#' on the the half-circle between -pi and pi.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable) or "circular" component of pseudo-observations.
#' They must be on the half-circle, i.e. theta must be in \eqn{[-\pi, \pi)}.
#' @param nbars \link[base]{numeric} \link[base]{integer}, the number of bins (bars)
#' in the histogram.
#'
#'
#' @return A '\code{\link[ggplot2]{ggplot}}' object.
#'
#' @examples set.seed(123)
#'
#' theta <- cylcop::rvonmisesmix(n = 100,
#'   mu = c(0, pi),
#'   kappa = c(5, 2),
#'   prop = c(4, 2)
#' )
#' plot1 <- plot_circ_hist(theta)
#'
#' @references \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{plot_joint_scat}()}.
#'
#' @export
#'
plot_circ_hist <- function(theta, nbars = 20) {
  # validate input
  tryCatch({
    check_arg_all(check_argument_type(
      theta,
      type = "numeric",
      lower = -pi,
      upper = pi
    )
    ,
    1)
    check_arg_all(check_argument_type(
      nbars,
      type = "numeric",
      length = 1,
      integer = T,
      lower = 2
    )
    ,
    1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })
  theta <- theta / pi
  theta_hist <-
    graphics::hist(theta,
                   breaks = seq(-1, 1, length.out = (nbars + 1)),
                   plot = FALSE)
  df_rect <-
    data.frame(
      xmin = theta_hist$breaks[1:(length(theta_hist$breaks) - 1)],
      xmax = theta_hist$breaks[2:length(theta_hist$breaks)],
      ymax = theta_hist$counts,
      ymin = rep(0, nbars)
    )

  breaks <- pretty(0:max(theta_hist$counts))
  max_plot <- max(breaks)
  min_plot <- -0.2 * max_plot

  df_line_x_major <- data.frame(
    x = seq(-0.75, 1, 0.25),
    xend = seq(-0.75, 1, 0.25),
    y = rep((0.75 * min_plot), 8),
    yend = rep((max_plot * 1.1), 8)
  )
  df_line_x_minor <- df_line_x_major
  df_line_x_minor$x <- df_line_x_minor$x - 0.125
  df_line_x_minor$xend <- df_line_x_minor$xend - 0.125


  p <- ggplot() +
    geom_segment(
      data = df_line_x_major,
      aes(
        x = .data$x,
        y = .data$y,
        xend = .data$xend,
        yend = .data$yend
      ),
      colour = "grey40",
      size = 0.5
    ) +
    geom_segment(
      data = df_line_x_minor,
      aes(
        x = .data$x,
        y = .data$y,
        xend = .data$xend,
        yend = .data$yend
      ),
      colour = "grey60",
      size = 0.2
    ) +
    geom_rect(
      data = df_rect,
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      color = "black",
      size = 0.5,
      fill = "grey60"
    ) +
    coord_polar(start = pi) +
    xlab("theta") + theme_bw() +
    scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
    scale_y_continuous(breaks = breaks, limits = c(min_plot, (max_plot *
                                                                1.2))) +
    geom_label(
      data = data.frame(
        x = c(-0.75,-0.5,-0.25, 0, 0.25, 0.5, 0.75, 1),
        y = rep((max_plot * 1.2), 8)
      ),
      aes(.data$x, .data$y),
      label = c(
        expression(-0.75 * pi,-0.5 * pi,-0.25 * pi, 0,
                   0.25 * pi, 0.5 * pi, 0.75 * pi, pi)
      ),
      size = 10 / 2.834646,
      label.size = 0
    ) +

    geom_hline(yintercept = breaks[seq(2, length(breaks), 2)],
               colour = "darkgreen",
               size = 0.2) +
    geom_hline(yintercept = breaks[seq(1, length(breaks), 2)],
               colour = "darkred",
               size = 0.2) +
    geom_label(
      data = data.frame(x = rep(-0.95, length(breaks[seq(2, length(breaks), 2)])), y =
                          breaks[seq(2, length(breaks), 2)]),
      aes(.data$x, .data$y),
      label = breaks[seq(2, length(breaks), 2)],
      size = 10 / 2.834646,
      label.size = 0,
      fill = "white",
      alpha = 0.7,
      color = "darkgreen"
    ) +
    geom_label(
      data = data.frame(x = rep(0.95, length(breaks[seq(1, length(breaks), 2)])), y =
                          breaks[seq(1, length(breaks), 2)]),
      aes(.data$x, .data$y),
      label = breaks[seq(1, length(breaks), 2)],
      size = 10 / 2.834646,
      label.size = 0,
      fill = "white",
      alpha = 0.7,
      color = "darkred"
    ) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  return(suppressWarnings(cowplot::ggdraw(p)))
}
