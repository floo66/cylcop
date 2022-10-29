#' Density, Distribution, Quantiles and Random Number Generation for the Wrapped
#' Cauchy Distribution
#'
#' The distribution function (\code{pwrappedcauchy()}) and quantiles
#'  (\code{qwrappedcauchy()}) of the wrapped Cauchy distribution cannot
#'  be obtained analytically. They are therefore missing in the
#' '\pkg{circular}' package and are obtained here numerically.
#' Random number generation (\code{rwrappedcauchy()}) and density
#' (\code{dwrappedcauchy()}) don't need a numerical
#' approximation and are provided here for consistency in parametrization
#'  with the other wrapped Cauchy functions.
#'
#' One could alternatively  convert \code{scale} to \code{rho} via
#' \code{rho = exp(-scale)} and use
#' \code{circular::\link[circular]{rwrappedcauchy}(theta, mu=location rho=rho)} or
#' \code{circular::\link[circular]{dwrappedcauchy}(theta, mu=location rho=rho)}.
#'
#' The wrapped Cauchy cdf, for which there is no analytical expression,
#' is calculated by wrapping the Cauchy distribution \eqn{K} times
#'  around the circle in each direction and summing the Cauchy cdfs at each point of
#'   the circle. Let \eqn{\Omega} follow a Cauchy distribution and
#'   \eqn{\Theta} a wrapped Cauchy distribution, where \eqn{\Theta} can take values
#'   \eqn{\theta \in [-\pi,\pi)}.
#'    \eqn{Pr(\Theta \le \theta)} is approximated as
#'    \deqn{\sum^K_{k=-K}Pr(\Omega \le \theta+2\pi k)-Pr(\Omega \le -\pi+2\pi k).}
#'    The quantiles are calculated by numerical inversion.
#'
#' @param location \link[base]{numeric} value, the mean of the distribution.
#' @param scale \link[base]{numeric} value, the parameter tuning the spread of the
#' density. It must be non-negative.
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
#' @return
#' \itemize{
#' \item{\code{dwrappedcauchy()}}{ gives a \link[base]{vector} of length \code{length(theta)}
#'  containing the density at \code{theta}.}
#' \item{\code{pwrappedcauchy()}}{ gives a
#' \link[base]{vector} of length \code{length(theta)} containing
#' the distribution function at the corresponding values of \code{theta}.}
#' \item{\code{qwrappedcauchy()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rwrappedcauchy()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples, i.e. angles in \eqn{[-\pi, \pi)}.}
#'}
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
#' @aliases rwrappedcauchy pwrappedcauchy dwrappedcauchy qwrappedcauchy
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
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(location,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
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
                           scale = 1) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(theta,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(location,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  rho <- exp(-scale)
  pdf <-
    (1 - rho ^ 2) / ((2 * pi) * (1 + rho ^ 2 - 2 * rho * cos(theta - location)))

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
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(theta,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(location,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(K,
                                      type="numeric",
                                      lower=1,
                                      integer=T)
                  ,1)
    check_arg_all(check_argument_type(check_prec,
                                      type="logical")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
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
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(location,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(K,
                                      type="numeric",
                                      lower=1,
                                      integer=T)
                  ,1)
    check_arg_all(check_argument_type(check_prec,
                                      type="logical")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  if (any(p < 0) || any(p > 1)) {
    stop(error_sound(),
         "The entries of p must be between 0 and 1.")
  }
  if (check_prec)
    check_precision(scale, K)

  qq = p
  for (i in seq_along(p))
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
