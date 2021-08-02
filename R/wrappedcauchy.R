#' Density, Distribution, Quantiles and Random Number Generation for the Wrapped
#' Cauchy Distribution
#'
#' The distribution function \code{pwrappedcauchy()} and quantiles
#'  \code{qwrappedcauchy()}of the wrapped Cauchy distribution can
#' not be obtained analytically. They are therefore missing in the
#' '\pkg{circular}' package and are obtained here numerically.
#' Random number generation \code{rwrappedcauchy()} and density
#' \code{dwrappedcauchy()} don't need a numerical
#' approximation and are provided here just for consistency in parametrization
#'  with the other wrapped Cauchy functions. One could also convert \code{scale} to \code{rho}
#' (\code{rho = exp(-scale)}) and use \code{circular::\link[circular]{rwrappedcauchy}(rho)}.
#'  In fact, for the density, one usually SHOULD convert \code{scale} to \code{rho}
#' (\code{rho = exp(-scale)}) and use \code{circular::\link[circular]{dwrappedcauchy}(rho)}
#'  which does not make a numerical approximation and is therefore faster than
#' \code{cylcop::dwrappedcauchy()}.
#'
#' The density is calculated by wrapping the Cauchy distribution \eqn{K} times
#'  around the circle in each direction and summing the density at each point of
#'   the circle. E.g. the density of the wrapped Cauchy distribution
#'   at angle \eqn{\theta} is calculated as  the sum of the Cauchy density at
#'   \eqn{\theta+2 \pi k}, where the integer \eqn{k} goes from \eqn{-K} to \eqn{K}.
#'   The distribution function is obtained similarly and the quantiles are calculated
#'   by numerical inversion.
#'
#' @param location \link[base]{numeric} value, the mean of the distribution.
#' @param scale \link[base]{numeric} value, the parameter tuning the spread of the
#' density.
#' @param theta \link[base]{numeric} \link[base]{vector} giving the angles where
#' the density or distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rwrappedcauchy()}.
#' @param K \link[base]{integer} value, the number of "wraps" used in each direction
#' to approximate the distribution.
#' @param check_prec \link[base]{logical}, whether to check if the precision of
#' the numerical approximation with the current parameters is higher than 99\%.
#'
#' @returns \code{dwrappedcauchy()}) and \code{pwrappedcauchy()}) give a
#' \link[base]{vector} of length \code{length(theta)} containing
#' the density or distribution function at the corresponding values of \code{theta}.
#' \code{qwrappedcauchy()} gives a \link[base]{vector} of length \code{length(p)}
#'containing the quantiles at the corresponding values of \code{p}.
#'\code{rwrappedcauchy()} generates a \link[base]{vector} of length \code{n}
#' containing the random samples, i.e. angles in \eqn{[-\pi, \pi)}.
#'
#' @examples set.seed(123)
#'
#' rwrappedcauchy(10, location = 0, scale =3)
#'
#' dwrappedcauchy(c(0.1, pi), location = pi, scale =2)
#' circular::dwrappedcauchy(circular::circular(c(0.1,pi)), mu = circular::circular(pi), rho =exp(-2))
#'
#' prob <- pwrappedcauchy(c(0.1, pi), location = pi, scale =2)
#' prob
#' qwrappedcauchy(prob, location = pi, scale =2)
#'
#' @name wrappedcauchy
#'
#' @seealso \code{circular::\link[circular]{dwrappedcauchy}()},
#' \code{circular::\link[circular]{rwrappedcauchy}()}.
#'
#' @aliases rwrappedcauchy pwrappedcauchy rdwrappedcauchy qwrappedcauchy
#'
NULL

# Random numbers
#'
#' @rdname wrappedcauchy
#' @export
#'
rwrappedcauchy <- function(n,
                           location = 0,
                           scale = 1) {
  cauchy_sample <- rcauchy(n, location = location, scale = scale)
  temp <- cbind(cos(cauchy_sample), sin(cauchy_sample))

  # returns x-y coordinates of points on a circle of radius 1, convert to angles:
  angle <- sign(temp[, 2]) * acos(temp[, 1])
  angle
}



#---------------Wrapped Cauchy density----------------
#
#' @rdname wrappedcauchy
#' @export
#'
dwrappedcauchy <- function(theta,
                           location = 0,
                           scale = 1,
                           K = 100,
                           check_prec = FALSE) {
  warning(
    warning_sound(),
    "dwrappedcauchy() gives a numerical approximation. For a faster, analytical result use circular::dwrappedcauchy()"
  )

  if (check_prec)
    check_precision(scale, K)

  pdf <- rep(0, length(theta))
  for (k in-K:K)
    pdf <- pdf + dcauchy(theta + 2 * pi * k, location = location, scale = scale)
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
                           check_prec = FALSE) {
  if (check_prec)
    check_precision(scale, K)

  cdf <- rep(0, length(theta))
  for (k in-K:K)
    cdf <- cdf +
    pcauchy(theta + 2 * pi * k, location = location, scale = scale) -
    pcauchy(-pi + 2 * pi * k, location = location, scale = scale)
  return(cdf)
}



#---------------Quantiles of the wrapped Cauchy distribution----------------
#
#' @rdname wrappedcauchy
#' @export
#'
qwrappedcauchy <- function(p,
                           location = 0,
                           scale = 1,
                           K = 100,
                           check_prec = FALSE) {
  if (any(p < 0) || any(p > 1)) {
    stop(error_sound(),
         "The entries of p must be between 0 and 1.")
  }
  if (check_prec)
    check_precision(scale, K)

  qq = p
  for (i in 1:length(p))
  {
    ff = function (theta) {
      pwrappedcauchy(theta,
                     location = location,
                     scale = scale,
                     K = K) - p[i]
    }
    qq[i] = tryCatch(
      uniroot(ff, lower = -pi, upper = pi)$root,
      # Becasue of the numerical approximation, the domain of qwrappedcauchy() does not go up to exactly 1
      # if it throws an error, it is becasue p is larger than 1-epsilon, where epsilon gets smaller with increasing K
      error = function(e) {
        return(pi)
      }
    )
  }
  return(qq)
}


# Check precision of \code{wrappedcauchy} functions
#
# Check that the error in the numerical approximation is smaller than 1%
#
# @param scale The numeric value of the parameter tuning the spread of the
#   density.
# @param K The limit used to approximate the distribution.
#
# @return NULL, warning message if the error is larger than 1%.
#
#
check_precision <- function(scale, K) {
  #the probability of drawing an angle smaller than pi from a wrapped cauchy with mean 0 should be 1
  cdf <- 0
  for (k in-K:K)
    cdf <- cdf +
      pcauchy(pi + 2 * pi * k, location = 0, scale = scale) -
      pcauchy(-pi + 2 * pi * k, location = 0, scale = scale)
  precision <- cdf
  if (precision < 0.99) {
    warning(
      warning_sound(),
      "With this value of K, the error of the approximation is larger than ",
      round((1 - precision) * 100, 1),
      "%"
    )
  }
  return(NULL)
}
