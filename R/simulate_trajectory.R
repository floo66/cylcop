
#' Generate a Trajectory with Correlated Step Lengths and Turn Angles
#'
#' The function draws values from a circular-linear bivariate distribution of
#' turn angles and step lengths specified by the marginal distributions and a
#' circular-linear copula. From the start
#' point (0,0) and the second (potentially user specified)
#' point, a trajectory is then built with these turn angles and step lengths.
#'
#' @param n \link[base]{integer}, number of trajectory steps to generate.
#' @param copula '\code{\linkS4class{cyl_copula}}' object.
#' @param marginal_circ named \link[base]{list} (for parametric estimates) or
#' a '\code{\link[circular]{density.circular}}' object (for kernel density estimates).
#' The output of function \code{\link{fit_angle}()} can be used here directly for
#' both cases.
#' @param marginal_lin named \link[base]{list} (for parametric estimates) or
#' a '\code{\link[stats]{density}}' object (for kernel density estimates).
#' The output of function \code{\link{fit_steplength}()} can be used here directly for
#' both cases.
#' @param ignore_first \link[base]{logical} value. If \code{ignore_first = TRUE} (default),
#' a trajectory of length \code{n+2} is generated and the first two steps of that
#'  trajectory are removed.
#' @param pos_2 (optional) \link[base]{numeric} \link[base]{vector} of length 2
#'  containing the coordinates of the second point in the trajectory.
#'  The first point is always at (0,0). If
#' no value is specified, the second point is obtained by going in a random direction
#' from the first point for a distance drawn from the marginal step length distribution.
#'
#' @details
#' Samples are drawn from the circular-linear copula and then transformed
#' using the quantile functions of the marginal circular and the marginal linear
#' distribution. To generate draws from any bivariate joint distribution (not
#' necessarily a circular-linear one) without also producing a trajectory,
#' the function \code{\link{rjoint}()} can be used.
#'
#' If entered "by hand", the named lists describing the parametric distributions
#' (\code{marginal_circ} and \code{marginal_lin}) must contain 2 entries:
#' \enumerate{
#'    \item{\code{name}:
#' a \link[base]{character} string denoting the name of the distribution.
#' For the circular distribution, it can be \code{"vonmises"}, \code{"vonmisesmix"}, or
#'   \code{"wrappedcauchy"}. For the linear distribution, it must be a
#'   string denoting the name of a linear distribution in the environment, i.e. the name of its
#'    distribution function without the "p",
#'   e.g. "norm" for normal distribution}
#'    \item{\code{coef}: For the circular distribution \code{coef} is a (named) \link[base]{list} of
#' parameters of the circular
#' marginal distribution as taken by the functions
#' \code{\link[circular]{qvonmises}()}, \code{\link{qvonmisesmix}()},
#' or \code{\link{qwrappedcauchy}()}. For the linear distribution, \code{coef} is
#' a named list containing the parameters of the distribution given in \code{"name"}.}
#' }
#'
#'
#' @return A \link[base]{data.frame} containing the trajectory. It has 6 columns
#' containing the x and y coordinates, the turn angles, the step lengths, and
#' the values sampled from the copula.
#'
#' @examples require(circular)
#' set.seed(123)
#'
#' traj <- traj_sim(n = 5,
#' copula = cyl_quadsec(0.1),
#' marginal_circ = list(name="vonmises",coef=list(0, 1)),
#' marginal_lin = list(name="weibull",coef=list(shape=3))
#' )
#'
#' traj
#'
#' angles <- rvonmisesmix(100,
#'   mu = c(0, pi),
#'   kappa = c(2, 3),
#'   prop = c(0.4, 0.6)
#' )
#' angles <- full2half_circ(angles)
#' bw <- opt_circ_bw(theta = angles, method = "nrd", kappa.est = "trigmoments")
#' marg_ang <- fit_angle(theta = angles, parametric = FALSE, bandwidth = bw)
#'
#' steplengths <- rlnorm(100, 0, 0.3)
#' marg_stepl <- fit_steplength(x = steplengths, parametric = "lnorm")
#'
#' traj_sim(n = 5,
#' copula = cyl_quadsec(0.1),
#' marginal_circ = marg_ang,
#' marginal_lin = marg_stepl,
#' ignore_first = FALSE,
#' pos_2 = c(5,5)
#' )
#'
#' @seealso \code{\link{traj_get}()},
#' \code{\link{fit_steplength}()}, \code{\link{fit_angle}()},
#' \code{\link{plot_track}()}, \code{\link{plot_cop_scat}()},
#' \code{\link{plot_joint_scat}()}, \code{\link{plot_joint_circ}()}.
#' @export
#'
traj_sim <-
  function(n,
           copula,
           marginal_circ,
           marginal_lin,
           ignore_first = TRUE,
           pos_2 = NULL) {

    #validate input
    tryCatch({
      check_arg_all(check_argument_type(n, type="numeric", length = 1, lower = 1),1)
      n <- as.integer(n)
      check_arg_all(check_argument_type(copula, type="cyl_copula"),1)
      check_arg_all(list(check_argument_type(marginal_circ,
                                        type="list"),
                    check_argument_type(marginal_circ,
                                        type="density.circular"))
                    ,2)
      check_arg_all(list(check_argument_type(marginal_lin,
                                        type="list"),
                    check_argument_type(marginal_lin,
                                        type="density"))
      ,2)
      check_arg_all(check_argument_type(ignore_first,
                                        type="logical")
                    ,1)
      check_arg_all(list(check_argument_type(pos_2,
                                             type="numeric",
                                             length = 2),
                         check_argument_type(pos_2,
                                             type="NULL"))
                    ,2)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )

    if("density.circular" %in% is(marginal_circ)){
      parameter_circ <- marginal_circ
      marginal_circ <- "dens"
    }else{
      if(!all(c("name", "coef") %in% names(marginal_circ))){
        stop(error_sound(),
             "marginal_circ must be a density.circular object or list containing
             the entries 'name' and 'coef'."
             )
      }
      if(!is.character(marginal_circ$name)){
        stop(error_sound(),
             "In marginal_circ: name must be of type character."
        )
      }
      if(!any(c("vonmises","wrappedcauchy", "vonmisesmix") %in% marginal_circ$name)){
        stop(error_sound(),
             "In marginal_circ: name must be \"vonmises\", \"wrappedcauchy\", or \"vonmisesmix\"."
        )
      }
      if(!is.list(marginal_circ$coef)){
        stop(error_sound(),
             "In marginal_circ: coef must be of type list."
        )
      }
      parameter_circ <- marginal_circ$coef
      marginal_circ <- marginal_circ$name
    }


    if("density" %in% is(marginal_lin)){
      parameter_lin <- marginal_lin
      marginal_lin <- "dens"
    }else{
      if(!all(c("name", "coef") %in% names(marginal_lin))){
        stop(error_sound(),
             "marginal_lin must be a density object or list containing
             the entries 'name' and 'coef'."
        )
      }
      if(!is.character(marginal_lin$name)){
        stop(error_sound(),
             "In marginal_lin: name must be of type character."
        )
      }
      if(!is.list(marginal_lin$coef)){
        stop(error_sound(),
             "In marginal_lin: coef must be of type list."
        )
      }
      parameter_lin <- marginal_lin$coef
      marginal_lin <- marginal_lin$name
    }



    #-----checks preparations and get parameters-------------------------------------------

    #get the  marginal distributions, densities, etc. into one list
    marg_circ <- get_marg(marginal_circ)
    marg_lin <- get_marg(marginal_lin)
    if(cylcop.env$silent==F){
      #echo what functions and parameters are used
      printCop(copula)
      if(marginal_circ != "dens"){
        message(
          sprintf(
            "%-18s %-18s parameters: %s",
            "\nCircular marginal:",
            marginal_circ,
            print_param(marg_circ$d, parameter_circ)
          ))
      }else{
        message("\nangular non-parametric density")
      }
      if(marginal_lin != "dens"){
        message(
          sprintf(
            "%-18s %-18s parameters: %s",
            "Linear marginal:",
            marginal_lin,
            print_param(marg_lin$r, parameter_lin)
          ),
          "\n")
      }else{
        message("\nlinear non-parametric density")
      }
    }
    #Test if the precision of the approximation for the wrapped cauchy marginal distribution is sufficient
    if(marginal_circ=="wrappedcauchy")
      tryCatch(do.call(marg_circ$p, c(0, parameter_circ)),
               warning = function(w) {
                 rlang::abort(conditionMessage(w))
               })

    #Give arbitrary first 2 positions and generate otherwise empty trajectory
    #prevp2 is position at step-2 and prevp1 is position at step-1
    if(!any(is.null(pos_2)) && length(pos_2)!=2){
      stop(error_sound(),
           "pos_2 must be a vector of length 2.")
    }
    if(any(is.null(pos_2))){
      dist <- do.call(marg_lin$r,c(1, parameter_lin))
      ang <- runif(1,-pi,pi)
      pos_2 <- c(dist*cos(ang), dist*sin(ang))
  }
    prevp2 <- c(0, 0)
    prevp1 <- pos_2

    if(ignore_first){
      n <- n+2
    }
    traj <-
      data.frame(
        pos_x = c(0, pos_2[1], rep(NA, n - 2)),
        pos_y = c(0, pos_2[2], rep(NA, n - 2)),
        angle = c(NA, 1, rep(NA, n - 2)),
        steplength = rep(NA, n),
        cop_u = rep(NA, n),
        cop_v = rep(NA, n)
      )


    #-----actually produce trajectory-------------------------------------------
    if(n<=100){
      step_start <- 3
      step_end <- n
    }else if (n<1002 & n>100){
      step_start <- c(3,101)
      step_end <- c(100,n)
    }else{
      step_start <- c(3,101,seq(1001,(n-1),1000))
      step_end <- unique(c(100,seq(1000,n,1000),n))
    }
    #start timer
    ptm <- proc.time()
    time <- 0

    for(i in seq_along(step_start)){
      cop_sample <- rcylcop((step_end[i]-step_start[i]+1), copula)
      if(marginal_lin!="dens") {step_vec <- do.call(marg_lin$q, c(list(p=cop_sample[,2]), parameter_lin))}
      else{step_vec <- do.call(marg_lin$q,list(cop_sample[,2], parameter_lin))}
      if(marginal_circ!="dens"){suppressWarnings(angle_vec <- do.call(marg_circ$q, c(list(p=cop_sample[,1]), parameter_circ)))}
      else{angle_vec <- do.call(marg_circ$q,list(cop_sample[,1], parameter_circ))}
      traj[step_start[i]:step_end[i],3] <- angle_vec
      traj[step_start[i]:step_end[i],4] <- step_vec
      traj[step_start[i]:step_end[i],5] <- cop_sample[,1]
      traj[step_start[i]:step_end[i],6] <- cop_sample[,2]

      step_in_batch <- 1
      x_pos_vec <- rep(0,length(step_start[i]:step_end[i]))
      y_pos_vec <- x_pos_vec
      for(step in step_start[i]:step_end[i]){

        cop_uv <- cop_sample[step_in_batch,,drop=F]

        if(any(cop_uv==1)){
          cop_uv[which(cop_uv>0.99999999999)] <- 0.99999999999
        }
        if(any(cop_uv==0)){
          cop_uv[which(cop_uv<0.00000000001)] <- 0.00000000001
        }
        #convert to correlated marg_lin and marg_circ distributions using the inverse cdf's
        steplength <- step_vec[step_in_batch]
        angle <- angle_vec[step_in_batch]


        #take the step and add it to trajectory
        point <- angstep2xy(angle, steplength, prevp1, prevp2)
        x_pos_vec[step_in_batch] <- point[1]
        y_pos_vec[step_in_batch] <- point[2]
        prevp2 <- prevp1
        prevp1 <- point

        step_in_batch <- step_in_batch+1
      }
      traj[step_start[i]:step_end[i],1] <- x_pos_vec
      traj[step_start[i]:step_end[i],2] <- y_pos_vec

      if(cylcop.env$silent==F){
        #get time for 100 steps of trajectory, if trajectory is at least that long, to set up the progress bar
        if (step_end[i] == 1000) {
          time <- (proc.time() - ptm)[3] %>% as.double()
          if ((time / 1000 * n) > 20) {
            message(
              "Estimated time to generate trajectory: ",
              floor((time / 1000 * n) / 60),
              " minutes, ",
              round((time / 1000 * n) - 60 * floor((time / 1000 * n) / 60)),
              " seconds"
            )
            pb <- utils::txtProgressBar(min = 1, max = n)
          }
        }
        if ((time / 1000 * n) > 20 && step_end[i] >= 1000){
          utils::setTxtProgressBar(pb, step_end[i])
          if(i==which(step_end>=n/2)[1]) waiting_sound()
        }
      }
    }
    if(cylcop.env$silent==F){
      #close progress bar
      if ((time / 1000 * n) > 20 && n >= 1000){
        close(pb)
        done_sound()
      }
    }

    if(ignore_first){
      traj <- traj[3:n,]
    }
    return(traj)
  }


#' Get a Trajectory from Coordinates
#'
#' The function calculates step lengths and turn angles from x- and y-coordinates
#' and calculates pseudo-observations from those step lengths and turn angles.
#'
#' @param x_coords \link[base]{vector} of \link[base]{numeric} values
#' containing the x-coordinates of a trajectory.
#' @param y_coords \link[base]{vector} of \link[base]{numeric} values
#' containing the y-coordinates of a trajectory.
#'
#'
#' @return A \link[base]{data.frame} containing the trajectory. It has 6 columns
#' containing the x and y coordinates, the turn angles, the step lengths, and
#' the pseudo-observations.
#'
#' @examples
#' set.seed(123)
#'
#' traj <- traj_sim(n = 5,
#' copula = cyl_quadsec(0.1),
#' marginal_circ = list(name="vonmises",coef=list(0, 1)),
#' marginal_lin = list(name="weibull",coef=list(shape=3))
#' )
#'
#' traj_from_coords <- traj_get(traj[,1], traj[,2])
#'
#'
#' @seealso \code{\link{traj_sim}()}.
#' @export
#'
traj_get <-
  function(x_coords, y_coords) {

#    validate input
    tryCatch({
      check_arg_all(check_argument_type(x_coords,
                                        type="numeric"),1)
      check_arg_all(check_argument_type(y_coords,
                                        type="numeric"),1)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )

    if(is.matrix(x_coords) || is.matrix(y_coords)){
      stop(
        error_sound(),
        "x_coords and y_coords must be vectors, not matrices."
      )
    }

    if(any(is.na(x_coords)) || any(is.na(y_coords))){
      stop(
        error_sound(),
        "x_coords and y_coords contain NA's."
      )
    }

    if(length(x_coords) != length(y_coords)){
      stop(
        error_sound(),
        "x_coords and y_coords must have the same length."
      )
    }

    if(length(x_coords)<3){
      stop(
        error_sound(),
        "Provide at least 3 coordinates."
      )
    }

    steplength <- x_coords
    steplength[1] <- NA
    steplength[2] <- sqrt((x_coords[2]-x_coords[1])^2+
                            (y_coords[2]-y_coords[1])^2)
    angle <- x_coords
    angle[1:2] <- NA
    for (i in 1:(length(x_coords)-2)) {
      steplength[i+1] <- sqrt((x_coords[i+2]-x_coords[i+1])^2+
                                (y_coords[i+2]-y_coords[i+1])^2)
      bear1 <- bearing(c(x_coords[i],y_coords[i]),
                       c(x_coords[i+1],y_coords[i+1]),
                       fullcirc = TRUE)
      bear2 <- bearing(c(x_coords[i+1],y_coords[i+1]),
                       c(x_coords[i+2],y_coords[i+2]),
                       fullcirc = TRUE)
      angle[i+2] <- full2half_circ((bear2-bear1)%% (2 * pi))
    }

    pseudo_obs <- copula::pobs(cbind(angle[-c(1,2)], steplength[-c(1,2)]),
                               ties.method = "average")


    traj <-
      data.frame(
        pos_x = x_coords,
        pos_y = y_coords,
        steplength = steplength,
        angle = angle,
        cop_v = c(NA,NA,pseudo_obs[,2]),
        cop_u = c(NA,NA,pseudo_obs[,1])
      )


    return(traj)
  }
