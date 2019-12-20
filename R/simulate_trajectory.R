
#' Generate a trajectory
#'
#' @param n Number of trajectory steps to generate.
#' @param copula A \code{cyl_copula} object.
#' @param marginal_circ A character string denoting the name of the circular
#'   distribution. "vonmises","wrappedcauchy", "dens" (for kernel density
#'   estimate), or "hnorm".
#' @param parameter_circ A list of parameters of the circular marginal dsitribution as taken by the functions
#'   \code{qvonmises()},\code{qwrappedcauchy()}, or \code{qhorm()}. For \code{marginal_circ = "dens"}, \code{parameter_circ} must be the density object,
#'   which can be obtained using \code{fit_angle(...,parametric = FALSE)}.
#' @param marginal_lin A character string denoting the name of the linear distribution. I.e. the name of its distribution function without the "p"
#'   e.g. "norm" for normal distribution.
#' @param parameter_lin  A list of parameters of the circular marginal dsitribution. For \code{marginal_circ = "dens"}, \code{parameter_circ} must be the density object,
#'   which can be obtained using \code{fit_steplength(...,parametric = FALSE)}.
#'
#' @return A data.frame containing the trajectory.
#' @export
#'
make_traj <-
  function(n,
           copula,
           marginal_circ = c("vonmises","wrappedcauchy", "dens", "hnorm"),
           parameter_circ = list(0, 1),
           marginal_lin,
           parameter_lin = list()) {

#-----checks preparations and get parameters-------------------------------------------

    #get the  marginal distributions, densities, etc. into one list
    marg_circ <- get_marg(marginal_circ)
    marg_lin <- get_marg(marginal_lin)

    #echo what functions and parameters are used
    cat("\n")
    printCop(copula)
    cat("\n")
    cat(
      sprintf(
        "%-18s %-20s parameters: %-60s",
        "Circular marginal:",
        marginal_circ,
        print_param(marg_circ$d, parameter_circ)
        ),
      "\n"
    )
    cat(
      sprintf(
        "%-18s %-20s parameters: %-60s",
        "Linear marginal:",
        marginal_lin,
        print_param(marg_lin$r, parameter_lin)
      ),
      "\n\n"
    )

    #Test if the precision of the approximation for the wrapped cauchy marginal distribution is sufficient
    if(marginal_circ=="wrappedcauchy")
      tryCatch(do.call(marg_circ$p, c(0, parameter_circ)),
              warning = function(w) {
                rlang::abort(conditionMessage(w))
              })

    #give warnings regarding approximations, with hnorm as marginal
    if (any(is(copula) == "cyl_gauss") && marginal_circ != "hnorm") {
      warning(cylcop::warning_sound(),
        "gaussian copula only gives at least approximatley correct results when used with hnorm as
        marginal circular distribution"
      )
    }
    if (!any(is(copula) == "cyl_gauss") && marginal_circ == "hnorm") {
      warning(cylcop::warning_sound(), "hnorm as marginal circular distribution only really makes sense with gauss copula")
    }
    if (marginal_circ == "hnorm") {
      prob <- 1 - do.call(extraDistr::phnorm, c(pi, parameter_circ))
      if (prob > 0.0001) {
        warning(cylcop::warning_sound(),
          "The probability to draw a value from hnorm that is larger than pi is ",
          round(100 * prob, 5),
          "%"
        )
      }
    }

    #Give arbitrary first 2 positions and generate otherwise empty trajectory
    #prevp2 is position at step-2 and prevp1 is position at step-1
    prevp2 <- c(0, 0)
    prevp1 <- c(1, 0)
    traj <-
      data.frame(
        pos_x = c(0, 1, rep(NA, n - 2)),
        pos_y = c(0, 0, rep(NA, n - 2)),
        steplength = c(NA, 1, rep(NA, n - 2)),
        angle = rep(NA, n),
        cop_v = rep(NA, n),
        cop_u = rep(NA, n)
      )


#-----actually produce trajectory-------------------------------------------

    #start timer
    ptm <- proc.time()
    time <- 0

    for (i in 3:n) {

      #get time for 100 steps of trajectory, if trajectory is at least that long, to set up the progress bar
      if (i == 100) {
        time <- (proc.time() - ptm)[3] %>% as.double()
        if ((time / 100 * n) > 5) {
          cat(
            "estimated time to generate trajectory: ",
            floor((time / 100 * n) / 60),
            " minutes, ",
            (time / 100 * n) - 60 * floor((time / 100 * n) / 60),
            " seconds\n"
          )
          pb <- utils::txtProgressBar(min = 1, max = n)
        }
      }
      if ((time / 100 * n) > 5 && i >= 100){
        utils::setTxtProgressBar(pb, i)
       if(i==n/2) cylcop::waiting_sound()
      }

      #get a sample of 2 correlated uniformly distributed random variables.
      #Subtract a small number to avoid to get "exactly" 1 (within machine precision) and therefore a steplength of INF
      cop_uv <- rCopula(1, copula)
      if(any(cop_uv==1)){
        cop_uv[which(cop_uv>0.99999999999)] <- 0.99999999999
      }
      if(any(cop_uv==0)){
        cop_uv[which(cop_uv<0.00000000001)] <- 0.00000000001
      }
      #convert to correlated marg_lin and marg_circ distributions using the inverse cdf's
      steplength <- do.call(marg_lin$q, c(cop_uv[2], parameter_lin))
      if (marginal_circ == "vonmises") {
        #We supress warnings, so we are not annoyed by the warning of the conversion to circular
        #As usual the circular package is a huge pile of equine manure, but I figured the following:
        #To get the quantile function of the full circle [0,2pi) distribution, use
        # angle <- suppressWarnings(do.call(qvonmises, c(cop_uv[1],parameter_circ, list(from=0))) )
        #you can then convert it to angles on the half circle. Or you use the the quantile function
        #of the half circle [-pi,pi) distribution directly, as below. Both approaches obviouysly give different results,
        #but you can convert them into each other by changing the copula for symmetry reasons. See text
        angle <-
          suppressWarnings(do.call(marg_circ$q,
                                   c(cop_uv[1], parameter_circ, list(from = -pi)
                                     )) - 2 * pi)
      }
      else{
        angle <- do.call(marg_circ$q, c(cop_uv[1], parameter_circ))
      }

      #with hnorm, we only get "right-turns" so we flip them to left turns with a probability of 0.5
      if (marginal_circ == "hnorm" && runif(1) < 0.5) {
        angle <- (-1) * angle
      }

      #take the step and add it to trajectory
      point <- angstep2xy(angle, steplength, prevp1, prevp2)
      traj[i, ] <- c(point, steplength, angle, cop_uv[2], cop_uv[1])
      prevp2 <- prevp1
      prevp1 <- point
    }

    #close progress bar
    if ((time / 100 * n) > 5 && n >= 100){
      close(pb)
      cylcop::done_sound()
    }
    return(traj)
    }
