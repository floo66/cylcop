#' Density, distribution, random number generation and quantiles of the wrapped
#' Cauchy distribution
#'
#' Distribution function and quantiles of the wrapped Cauchy distribution can
#' not be obtained analytically. They are therefore missing in the
#' \code{circular} package and are obtained here numerically.
#' \code{rwrappedcauchy()} and \code{dwrappedcauchy()} don't need a numerical
#' approximation and are just here for consistency in parametrization with the
#' other wrappedcauchy functions. You could also convert scale to rho
#' (rho=exp(-scale)) and use \code{circular::rwrappedcauchy()}. In fact, for
#' \code{dwrappedcauchy()}, you usually SHOULD convert scale to rho
#' (rho=exp(-scale)) and use \code{circular::dwrappedcauchy()} which does not
#' make a numerical approximation and is faster than
#' \code{cylcop::dwrappedcauchy()}.
#' @param location The numeric value of the mean of the distribution.
#' @param n Number of random samples to be generated.
#' @param scale The numeric value of the parameter tuning the spread of the
#'   density.
#' @name wrappedcauchy
#' @export
#'
rwrappedcauchy <- function(n,
                           location = 0,
                           scale = 1
                           ) {
  cauchy_sample <- rcauchy(n,location = location, scale = scale)
  temp <- cbind(cos(cauchy_sample), sin(cauchy_sample))

  # returns x-y coordinates of points on a circle of radius 1, convert to angles:
  angle <- sign(temp[, 2]) * acos(temp[, 1])
}



#---------------Wrapped Cauchy density----------------
#
#' @rdname wrappedcauchy
#' @param K The limit used to approximate the distribution.
#' @param theta A numeric vector hoding the angles at which to calculate the density or distribution.
#' @param check_prec Logical, whether to check if precision of numerical approximation
#' with the current parameters is higher than 99%
#' @export
#'
dwrappedcauchy <- function(theta,
                           location = 0,
                           scale = 1,
                           K = 100,
                           check_prec = F) {
  warning(cylcop::warning_sound(),
          "dwrappedcauchy() gives a numerical approximation. For a faster, analytical result use circular::dwrappedcauchy()"
          )

 if(check_prec) check_precision(scale, K)

  pdf<-rep(0,length(theta))
  for (k in -K:K) pdf<-pdf+dcauchy(theta+2*pi*k,location = location, scale = scale)
  return(pdf)
}



#---------------Wrapped Cauchy distribution function----------------
#
#' @rdname wrappedcauchy
#'
#' @export
#'
pwrappedcauchy <- function(theta,
                           location = 0,
                           scale = 1,
                           K = 100,
                           check_prec = F) {

  if(check_prec) check_precision(scale, K)

  cdf<-rep(0,length(theta))
  for (k in -K:K) cdf<-cdf+
    pcauchy(theta+2*pi*k,location = location, scale = scale)-
    pcauchy(-pi+2*pi*k,location = location, scale = scale)
  return(cdf)
}



#---------------Quantiles of the wrapped Cauchy distribution----------------
#
#' @rdname wrappedcauchy
#' @param p vector of probabilities
#'
#' @export
#'
qwrappedcauchy <- function(p,
                           location = 0,
                           scale = 1,
                           K = 100,
                           check_prec = F) {

  if(any(p<0) || any(p>1)){
    stop(cylcop::error_sound(), "The entries of p must be between 0 and 1.")
  }
  if(check_prec)  check_precision(scale, K)

  qq=p
  for (i in 1:length(p))
  {
    ff=function (theta) {pwrappedcauchy2(theta, location = location, scale= scale, K=K)-p[i]}
    qq[i]=tryCatch(uniroot(ff,lower=-pi,upper=pi)$root,
                   # Becasue of the numerical approximation, the domain of qwrappedcauchy() does not go up to exactly 1
                   # if it throws an error, it is becasue p is larger than 1-epsilon, where epsilon gets smaller with increasing K
                   error = function(e) {
                     return(pi)
                   }
    )
  }
  return(qq)
}


#' Check precision of \code{wrappedcauchy} functions
#'
#' Check that the error in the numerical approximation is smaller than 1%
#'
#' @param scale The numeric value of the parameter tuning the spread of the
#'   density.
#' @param K The limit used to approximate the distribution.
#'
#' @return NULL, warning message if the error is larger than 1%.
#' @export
#'
check_precision <- function(scale, K){

  #the probability of drawing an angle smaller than pi from a wrapped cauchy with mean 0 should be 1
     cdf<-rep(0,length(theta))
   for (k in -K:K) cdf<-cdf+
       pcauchy(theta+2*pi*k,location = location, scale = scale)-
       pcauchy(-pi+2*pi*k,location = location, scale = scale)
   precision <- cdf
  if (precision < 0.99) {
    warning(
      cylcop::warning_sound(),
      "With this value of K, the error of the approximation is larger than ",
      round((1 - precision) * 100, 1),
      "%"
    )
  }
   return(NULL)
}
