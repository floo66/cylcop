#' Density, Distribution, Random Number Generation and Quantiles of Kernel Density Estimates
#'
#' Calculate the density (\code{ddens()}), the distribution (\code{pdens()}),
#'  the quantiles (\code{qdens()}) and generate random
#' samples (\code{rdens()}) of a kernel density estimate as returned by
#' \code{\link{fit_angle}()} or \code{\link{fit_steplength}()}.
#'
#' @param density a '\code{\link[stats]{density}}' object (for linear kernel density estimates)
#' or a '\code{\link[circular]{density.circular}}' object
#' (for circular kernel density estimates) containing information about the kernel density
#'   estimate. These objects can be obtained using
#'   \code{\link{fit_angle}(..., parametric = FALSE)} or
#'   \code{\link{fit_steplength}(..., parametric = FALSE)}.
#' @param x \link[base]{numeric} \link[base]{vector} giving the points where
#' the density or distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rdens()}.
#'
#' @returns
#' \code{ddens()} and \code{pdens()} give a \link[base]{vector} of length \code{length(x)} containing
#' the density or distribution function at the corresponding values of \code{x}.
#' \code{qdens()} gives a \link[base]{vector} of length \code{length(p)} containing
#' the quantiles at the corresponding values of \code{p}. The function \code{rdens()}
#' generates a \link[base]{vector} of length \code{n} containing the random samples.
#'
#' @examples set.seed(123)
#'
#' steps <- rweibull(10, shape=3)
#' dens <- fit_steplength(x = steps, parametric = FALSE)
#' ddens(c(0.1,0.3), dens)
#' pdens(c(0.1,0.3), dens)
#' qdens(c(0.1,0.3), dens)
#' rdens(4, dens)
#'
#' angles <- full2half_circ(
#'   circular::rvonmises(10, mu = circular::circular(0), kappa = 2)
#' )
#' dens <- fit_angle(theta = angles, parametric = FALSE)
#' ddens(c(0.1,0.3), dens)
#' pdens(c(0.1,0.3), dens)
#' qdens(c(0.1,0.3), dens)
#' rdens(4, dens)
#'
#' @seealso \code{\link{fit_angle}()}, \code{\link{fit_steplength}()},
#' \code{\link{fit_steplength}()}.
#'
#' @name dens
#' @aliases ddens pdens qdens rdens
#'
#
NULL





# Density of Kernel Density Estimate
#'
#' @rdname dens
#' @export
#'
rdens <- function(n, density) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      length=1,
                                      integer=T,
                                      lower=1)
                  ,1)
    check_arg_all(list(check_argument_type(density,
                                           type="density.circular"),
                       check_argument_type(density,
                                           type="density"))
    ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  # Draw a points from the domain of the KDE
  sample <- sample(as.double(density$x),
         n,
         replace = TRUE,
         prob = density$y)

  #Now draw from the kernel distributions at these points,
  #with kappa (for vonmises) or sd (for Gaussian) equal to the bandwidth
  kernel<-get_marg(density$kernel)
  if("density.circular" %in% is(density)){
   out <-  map(sample, ~full2half_circ(kernel$r(1, circular::circular(.x), density$bw))) %>% unlist()
  }else{
    out <- map(sample, ~kernel$r(1, .x, density$bw)) %>% unlist()
  }
  return(out)
}



# Probability density from a kernel density estimate
#
#' @rdname dens
#' @export
#'
ddens <- function(x, density) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(list(check_argument_type(density,
                                           type="density.circular"),
                       check_argument_type(density,
                                           type="density"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  dens_interpol<-stats::approxfun(density$x,density$y)
  map_dbl(x,~dens_interpol(.x))
}



# Distribution function from a kernel density estimate
#'
#' @rdname dens
#' @export
#'
pdens <- function(x, density) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(list(check_argument_type(density,
                                           type="density.circular"),
                       check_argument_type(density,
                                           type="density"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  # could do something like this (below) to precisely calculate the cdf, but this is quite slow,
  # better use linear interpolation
  # rowMeans(pnorm(outer(x,na.omit(traj$steplength),"-"),0,density$bw))

  density$x <- as.double(density$x)
  density$y <- as.double(density$y)
  delx <-
    (density$x[length(density$x)] - density$x[1]) / length(density$x)

  #possible speedup by saving preious values of sum and just adding a term, already fast enough, though
  #right Riemann sum, then set the first value, cdf(min(x))=0
  cdf<-map_dbl(seq(1,length(density$y)-1,1),~sum(density$y[1:.x]) * delx)
  cdf <- c(0,cdf,1)

  #add one more x-value, its cdf is 1
  density$x<-c(density$x,2*density$x[length(density$x)]-density$x[length(density$x)-1])
  cdf_interpol <- stats::approxfun(c(density$x,Inf),c(cdf,1))
  map_dbl(x,~cdf_interpol(.x))
}



# Quantiles of a kernel density estimate
#'
#' @rdname dens
#' @export
#'
qdens <- function(p, density) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(list(check_argument_type(density,
                                           type="density.circular"),
                       check_argument_type(density,
                                           type="density"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  density$x <- as.double(density$x)
  density$y <- as.double(density$y)
  delx <-
    (density$x[length(density$x)] - density$x[1]) / length(density$x)

  bin_search <- function(p) {
    l <- 0L
    r <- length(density$x)
    while (abs(r - l) > 1) {
      m <- floor((l + r) / 2)
      prob <- sum(density$y[1:m]) * delx
      if (prob < p)
        l <- m
      else if (prob > p)
        r <- m
    }

    cdf<-function(index){
      if(index==0L) return(0)
      else if(index==length(density$y)) return(1)
      else return(sum(density$y[1:index]) * delx)
      }

    #add one more x-value, its cdf is 1
    density$x<-c(density$x,2*density$x[length(density$x)]-density$x[length(density$x)-1])

    #linear interpolation, faster than precise calculation, see pdens()
    q_interpol <- stats::approxfun(c(cdf(l), cdf(r)),c(density$x[l+1],density$x[r+1]))
    return(q_interpol(p))
  }
  map_dbl(p,  ~ bin_search(.x))
}
