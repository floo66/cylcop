
#' Generate a Trajectory with Correlated Step Lengths and Turn Angles
#'
#' The function draws values from a circular-linear bivariate distribution of
#' turn angles and step lengths specified by the marginal distributions and a
#' circular-linear copula. Samples are drawn from the copula and then transformed
#' using the quantile functions of the marginal distributions. From the start
#' point (0,0) and the second (user specified)
#' point, a trajectory is then built with these turn angles and step lengths.
#'
#' @param n \link[base]{integer}, number of trajectory steps to generate.
#' @param copula '\code{\linkS4class{cyl_copula}}' object.
#' @param marginal_circ \link[base]{character} string denoting the name of the circular
#'   distribution. It can be \code{"vonmises"}, \code{"mixedvonmises"},
#'   \code{"wrappedcauchy"}, or \code{"dens"} (for kernel density estimate).
#' @param parameter_circ (named) \link[base]{list} of parameters of the circular
#' marginal distribution as taken by the functions
#' \code{\link[circular]{qvonmises}()}, \code{\link{qmixedvonmises}()},
#' or \code{\link{qwrappedcauchy}()}. If \code{marginal_circ = "dens"},
#' \code{parameter_circ} must be a named \link[base]{list}, containing information
#' on the kernel density estimate, which can be obtained
#' using \code{\link{fit_angle}(...,parametric = FALSE)}.
#' @param marginal_lin \link[base]{character} string denoting the name of the
#' linear distribution, i.e. the name of its distribution function without the "p",
#'   e.g. "norm" for normal distribution, or \code{"dens"} (for kernel density estimate).
#' @param parameter_lin  (named) \link[base]{list} of parameters of the linear
#' marginal distribution. For \code{marginal_lin = "dens"}, \code{parameter_lin}
#'  must be a named \link[base]{list}, containing information
#' on the kernel density estimate, which can be obtained
#' using \code{\link{fit_steplength}(...,parametric = FALSE)}.
#' @param pos_2 \link[base]{numeric} \link[base]{vector} containing the coordinates
#' of the second point in the trajectory. The first point is always at (0,0).
#'
#' @return A \link[base]{data.frame} containing the trajectory. It has 6 columns
#' containing the x and y coordintates, the step lengths, the turn angles, and
#' the values sampled from the copula.
#'
#' @examples require(circular)
#' set.seed(123)
#'
#' traj <- make_traj(5,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape=3)
#' )
#' traj
#'
#' angles <- circular::rmixedvonmises(100,
#'   mu1 = circular::circular(0),
#'   mu2 = circular::circular(pi),
#'   kappa1 = 2,
#'   kappa2 = 3,
#'   prop = 0.4
#' )
#' angles <- full2half_circ(angles)
#' bw <- opt_circ_bw(theta = angles,loss = "adhoc", kappa.est = "trigmoments")
#' dens <- fit_angle(theta = angles, parametric = FALSE, bandwidth = bw)
#' make_traj(5,
#'   copula = cyl_quadsec(0.1),
#'   marginal_circ = "dens",
#'   parameter_circ = dens,
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape=3),
#'   pos_2 = c(5,5)
#' )
#'
#' @seealso \code{\link{fit_steplength}()}, \code{\link{fit_angle}()},
#' \code{\link{traj_plot}()}, \code{\link{cop_scat_plot}()},
#' \code{\link{scat_plot}()}, \code{\link{circ_plot}()}.
#' @export
#'
make_traj <-
  function(n,
           copula,
           marginal_circ = c("vonmises","wrappedcauchy", "mixedvonmises", "dens"),
           parameter_circ,
           marginal_lin,
           parameter_lin,
           pos_2 = c(1, 0)) {

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
    prevp2 <- c(0, 0)
    prevp1 <- pos_2
    traj <-
      data.frame(
        pos_x = c(0, pos_2[1], rep(NA, n - 2)),
        pos_y = c(0, pos_2[2], rep(NA, n - 2)),
        steplength = c(NA, 1, rep(NA, n - 2)),
        angle = rep(NA, n),
        cop_v = rep(NA, n),
        cop_u = rep(NA, n)
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

    for(i in 1:length(step_start)){
      cop_sample <- rcylcop((step_end[i]-step_start[i]+1), copula)
      if(marginal_lin!="dens") {step_vec <- do.call(marg_lin$q, c(list(p=cop_sample[,2]), parameter_lin))}
      else{step_vec <- do.call(marg_lin$q,list(cop_sample[,2], parameter_lin))}
      if(marginal_circ!="dens"){suppressWarnings(angle_vec <- do.call(marg_circ$q, c(list(p=cop_sample[,1]), parameter_circ)))}
      else{angle_vec <- do.call(marg_circ$q,list(cop_sample[,1], parameter_circ))}
      traj[step_start[i]:step_end[i],3] <- step_vec
      traj[step_start[i]:step_end[i],4] <- angle_vec
      traj[step_start[i]:step_end[i],5] <- cop_sample[,2]
      traj[step_start[i]:step_end[i],6] <- cop_sample[,1]

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


        # with hnorm, we only get "right-turns" so we flip them to left turns with a probability of 0.5
        if (marginal_circ == "hnorm" && runif(1) < 0.5) {
          angle <- (-1) * angle
        }


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
    return(traj)
  }
