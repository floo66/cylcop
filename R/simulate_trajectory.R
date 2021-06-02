
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
#' @param parameter_lin  A list of parameters of the linear marginal dsitribution. For \code{marginal_circ = "dens"}, \code{parameter_circ} must be the density object,
#'   which can be obtained using \code{fit_steplength(...,parametric = FALSE)}.
#' @return A data.frame containing the trajectory.
#' @export
#'
make_traj <-
  function(n,
           copula,
           marginal_circ = c("vonmises","wrappedcauchy", "mixedvonmises", "dens", "hnorm"),
           parameter_circ = list(0, 1),
           marginal_lin,
           parameter_lin = list()) {
    #-----checks preparations and get parameters-------------------------------------------

    #get the  marginal distributions, densities, etc. into one list
    marg_circ <- get_marg(marginal_circ)
    marg_lin <- get_marg(marginal_lin)
    if(cylcop.env$silent==F){
      #echo what functions and parameters are used
      printCop(copula)
      message(
        sprintf(
          "%-18s %-18s parameters: %s",
          "\nCircular marginal:",
          marginal_circ,
          print_param(marg_circ$d, parameter_circ)
        ))
      message(
        sprintf(
          "%-18s %-18s parameters: %s",
          "Linear marginal:",
          marginal_lin,
          print_param(marg_lin$r, parameter_lin)
        ),
        "\n")
    }
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
      cop_sample <- rCopula((step_end[i]-step_start[i]+1), copula)
      step_vec <- do.call(marg_lin$q, c(list(p=cop_sample[,2]), parameter_lin))
      angle_vec <- do.call(marg_circ$q, c(list(p=cop_sample[,1]), parameter_circ))
      traj[step_start[i]:step_end[i],3] <- step_vec
      traj[step_start[i]:step_end[i],4] <- angle_vec
      traj[step_start[i]:step_end[i],5] <- cop_sample[2,]
      traj[step_start[i]:step_end[i],6] <- cop_sample[1,]

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
        if (step_end[i] == 100) {
          time <- (proc.time() - ptm)[3] %>% as.double()
          if ((time / 100 * n) > 20) {
            message(
              "Estimated time to generate trajectory: ",
              floor((time / 100 * n) / 60),
              " minutes, ",
              round((time / 100 * n) - 60 * floor((time / 100 * n) / 60)),
              " seconds"
            )
            pb <- utils::txtProgressBar(min = 1, max = n)
          }
        }
        if ((time / 100 * n) > 20 && step_end[i] >= 100){
          utils::setTxtProgressBar(pb, step_end[i])
          if(i==which(step_end>=n/2)[1]) cylcop::waiting_sound()
        }
      }
    }
    if(cylcop.env$silent==F){
      #close progress bar
      if ((time / 100 * n) > 20 && n >= 100){
        close(pb)
        cylcop::done_sound()
      }
    }
    return(traj)
  }
