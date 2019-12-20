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
#' @param ... Additional arguments (NOT USED ???)
#' @name wrappedcauchy
#' @export
#'
rwrappedcauchy <- function(n,
                           location = 0,
                           scale = 1,
                           ...) {
  temp <- rwrappedg(n = n,
                    "cauchy",
                    location = location,
                    scale = scale)

  # returns x-y coordinates of points on a circle of radius 1, conver to angles:
  angle <- sign(temp[, 2]) * acos(temp[, 1])
}



#---------------Wrapped Cauchy density----------------
#
#' @rdname wrappedcauchy
#' @param K The limit used to approximate the distribution.
#' @param theta Anumeric vector hoding the angles at which to calculate the density or distribution.
#'
#' @export
#'
dwrappedcauchy <- function(theta,
                           location = 0,
                           scale = 1,
                           K = 100) {
  warning(cylcop::warning_sound(),
          "dwrappedcauchy() gives a numerical approximation. For a faster, analytical result use circular::dwrappedcauchy()"
          )

  check_precision(scale, K)

  dwrappedg(theta,
            "cauchy",
            K = K,
            location = location,
            scale = scale)
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
                           K = 100) {

  check_precision(scale, K)

  pwrappedg(theta,
            "cauchy",
            K = K,
            location = location,
            scale = scale)
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
                           K = 100) {

  if(any(p<0) || any(p>1)){
    stop(cylcop::error_sound(), "The entries of p must be between 0 and 1.")
  }
  check_precision(scale, K)

  q <- map_dbl(p, function(p) {
    tryCatch(
      qwrappedg(
        p,
        "cauchy",
        K = K,
        location = location,
        scale = scale
      ),

      # Becasue of the numerical approximation, the domain of qwrappedg() does not go up to exactly 1
      # if it throws an error, it is becasue p is larger than 1-epsilon, where epsilon gets smaller with increasing K
      error = function(e) {
        return(pi)
      }
    )
  })
  return(q)
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
   precision <-
    pwrappedg(pi,
              "cauchy",
              K = K,
              location = 0,
              scale = scale)
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
