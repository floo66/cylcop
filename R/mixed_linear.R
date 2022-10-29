#' Density, Distribution, Quantiles and Random Number Generation for the mixed normal
#' distribution
#'
#' The number of components in the mixed normal distribution is specified by the length
#' of the parameter vectors. The quantiles are numerically obtained from the distribution function using
#' monotone cubic splines.
#'
#' @param x \link[base]{numeric} \link[base]{vector} giving the points where
#' the density function is evaluated.
#' @param q \link[base]{numeric} \link[base]{vector} giving the quantiles where
#' the distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rnormmix()}.
#' @param mu \link[base]{numeric} \link[base]{vector} holding the means of the components.
#' @param sigma \link[base]{numeric} \link[base]{vector} holding the standard
#' deviations of the components.
#' @param prop \link[base]{numeric} \link[base]{vector}, holding the mixing proportions
#' of the components.
#'
#' @return
#' \itemize{
#' \item{\code{dnormmix()}}{ gives a \link[base]{vector} of length \code{length(x)}
#'  containing the density at \code{x}.}
#' \item{\code{pnormmix()}}{ gives a
#' \link[base]{vector} of length \code{length(q)} containing
#' the distribution function at the corresponding values of \code{q}.}
#' \item{\code{qnormmix()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rnormmix()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples.}
#'}
#'
#'
#' @examples
#'
#' rnormmix(10, mu = c(0, 3, 7), sigma = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' dnormmix(c(0, 2, 1), mu = c(0, 3), sigma = c(2, 2), prop = c(0.6, 0.4))
#'
#' prob <- pnormmix(c(0.1, 7), mu = c(0, 3, 7), sigma = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#' prob
#' qnormmix(prob, mu = c(0, 3, 7), sigma = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' @name normmix
#'
#'
#' @aliases rnormmix pnormmix dnormmix qnormmix
#'
NULL



#' Density, Distribution, Quantiles and Random Number Generation for the mixed Weibull
#' distribution
#'
#' The number of components in the mixed Weibull distribution is specified by the length
#' of the parameter vectors. The quantiles are numerically obtained from the distribution function using
#' monotone cubic splines.
#'
#' @param x \link[base]{numeric} \link[base]{vector} giving the points where
#' the density function is evaluated.
#' @param q \link[base]{numeric} \link[base]{vector} giving the quantiles where
#' the distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rweibullmix()}.
#' @param shape \link[base]{numeric} \link[base]{vector} holding the shape parameter
#' of the components.
#' @param scale \link[base]{numeric} \link[base]{vector} holding the scale parameter
#' of the components.
#' @param prop \link[base]{numeric} \link[base]{vector}, holding the mixing proportions
#' of the components.
#'
#' @return
#' \itemize{
#' \item{\code{dweibullmix()}}{ gives a \link[base]{vector} of length \code{length(x)}
#'  containing the density at \code{x}.}
#' \item{\code{pweibullmix()}}{ gives a
#' \link[base]{vector} of length \code{length(q)} containing
#' the distribution function at the corresponding values of \code{q}.}
#' \item{\code{qweibullmix()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rweibullmix()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples.}
#'}
#'
#' @examples
#'
#' rweibullmix(10, shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' dweibullmix(c(0, 2, 1), shape = c(1, 3), scale = c(2, 2), prop = c(0.6, 0.4))
#'
#' prob <- pweibullmix(c(0.1, 7), shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#' prob
#' qweibullmix(prob, shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' @name weibullmix
#'
#'
#' @aliases rweibullmix pweibullmix dweibullmix qweibullmix
#'
NULL



#' Density, Distribution, Quantiles and Random Number Generation for the mixed log-normal
#' distribution
#'
#' The number of components in the mixed log-normal distribution is specified by the length
#' of the parameter vectors. The quantiles are numerically obtained from the distribution function using
#' monotone cubic splines.
#'
#' @param x \link[base]{numeric} \link[base]{vector} giving the points where
#' the density function is evaluated.
#' @param q \link[base]{numeric} \link[base]{vector} giving the quantiles where
#' the distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rlnormmix()}.
#' @param meanlog \link[base]{numeric} \link[base]{vector} holding the means
#' of the components on the log scale.
#' @param sdlog \link[base]{numeric} \link[base]{vector} holding the standard
#' deviations of the components on the log scale.
#' @param prop \link[base]{numeric} \link[base]{vector}, holding the mixing proportions
#' of the components.
#'
#' @return
#' \itemize{
#' \item{\code{dlnormmix()}}{ gives a \link[base]{vector} of length \code{length(x)}
#'  containing the density at \code{x}.}
#' \item{\code{plnormmix()}}{ gives a
#' \link[base]{vector} of length \code{length(q)} containing
#' the distribution function at the corresponding values of \code{q}.}
#' \item{\code{qlnormmix()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rlnormmix()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples.}
#'}
#'
#' @examples
#'
#' rlnormmix(10, meanlog = c(1, 3, 7), sdlog = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' dlnormmix(c(0, 2, 1), meanlog = c(1, 3), sdlog = c(2, 2), prop = c(0.6, 0.4))
#'
#' prob <- plnormmix(c(0.1, 7), meanlog = c(1, 3, 7), sdlog = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#' prob
#' qlnormmix(prob, meanlog = c(1, 3, 7), sdlog = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' @name lnormmix
#'
#'
#' @aliases rlnormmix plnormmix dlnormmix qlnormmix
#'
NULL



#' Density, Distribution, Quantiles and Random Number Generation for the mixed gamma
#' distribution
#'
#' The number of components in the mixed gamma distribution is specified by the length
#' of the parameter vectors. The quantiles are numerically obtained from the distribution function using
#' monotone cubic splines.
#'
#' @param x \link[base]{numeric} \link[base]{vector} giving the points where
#' the density function is evaluated.
#' @param q \link[base]{numeric} \link[base]{vector} giving the quantiles where
#' the distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rgammamix()}.
#' @param shape \link[base]{numeric} \link[base]{vector} holding the shape parameter
#' of the components.
#' @param scale \link[base]{numeric} \link[base]{vector} holding the scale parameter
#' of the components.
#' @param rate \link[base]{numeric} \link[base]{vector} an alternative way to specify the scale
#' (\code{scale = 1 / rate}).
#' @param prop \link[base]{numeric} \link[base]{vector}, holding the mixing proportions
#' of the components.
#'
#' @return
#' \itemize{
#' \item{\code{dgammamix()}}{ gives a \link[base]{vector} of length \code{length(x)}
#'  containing the density at \code{x}.}
#' \item{\code{pgammamix()}}{ gives a
#' \link[base]{vector} of length \code{length(q)} containing
#' the distribution function at the corresponding values of \code{q}.}
#' \item{\code{qgammamix()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rgammamix()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples.}
#'}
#'
#' @examples
#'
#' rgammamix(10, shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' dgammamix(c(0, 2, 1), shape = c(1, 3), rate = c(2, 2), prop = c(0.6, 0.4))
#'
#' prob <- pgammamix(c(0.1, 7), shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#' prob
#' qgammamix(prob, shape = c(1, 3, 7), scale = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' @name gammamix
#'
#'
#' @aliases rgammamix pgammamix dgammamix qgammamix
#'
NULL



# Random numbers
#'
#' @rdname normmix
#' @export
#'
rnormmix <- function(n, mu, sigma, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(sigma,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rlinearmix(n,param1=mu, param2=sigma, prop=prop, func=rnorm)
}

# Density
#'
#' @rdname normmix
#' @export
#'
dnormmix <- function(x, mu, sigma, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(sigma,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=x,param1=mu, param2=sigma, prop=prop, func=dnorm)
}

# Distribution
#'
#' @rdname normmix
#' @export
#'
pnormmix <- function(q, mu, sigma, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(q,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(sigma,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=q,param1=mu, param2=sigma, prop=prop, func=pnorm)
}

# Quantiles
#'
#' @rdname normmix
#' @export
#'
qnormmix <- function(p, mu, sigma, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(sigma,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  qlinearmix(p, param1=mu, param2=sigma, prop=prop, lowlim=-Inf, npts=1000, envir=cylcop_normmix.env, func=pnorm)
}



# Random numbers
#'
#' @rdname weibullmix
#' @export
#'
rweibullmix <- function(n, shape, scale, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rlinearmix(n,param1=shape, param2=scale, prop=prop, func=rweibull)
}

# Density
#'
#' @rdname weibullmix
#' @export
#'
dweibullmix <- function(x, shape, scale, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=x,param1=shape, param2=scale, prop=prop, func=dweibull)
}

# Distribution
#'
#' @rdname weibullmix
#' @export
#'
pweibullmix <- function(q, shape, scale, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(q,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=q,param1=shape, param2=scale, prop=prop, func=pweibull)
}

# Quantiles
#'
#' @rdname weibullmix
#' @export
#'
qweibullmix <- function(p, shape, scale, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  qlinearmix(p, param1=shape, param2=scale, prop=prop, lowlim=0, npts=1000, envir=cylcop_weibullmix.env, func=pweibull)
}



# Random numbers
#'
#' @rdname lnormmix
#' @export
#'
rlnormmix <- function(n, meanlog, sdlog, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(meanlog,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(sdlog,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rlinearmix(n,param1=meanlog , param2=sdlog, prop=prop, func=rlnorm)
}

# Density
#'
#' @rdname lnormmix
#' @export
#'
dlnormmix <- function(x, meanlog , sdlog, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(meanlog,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(sdlog,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=x,param1=meanlog , param2=sdlog, prop=prop, func=dlnorm)
}

# Distribution
#'
#' @rdname lnormmix
#' @export
#'
plnormmix <- function(q, meanlog , sdlog, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(q,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(meanlog,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(sdlog,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  pdlinearmix(point=q,param1=meanlog , param2=sdlog, prop=prop, func=plnorm)
}

# Quantiles
#'
#' @rdname lnormmix
#' @export
#'
qlnormmix <- function(p, meanlog , sdlog, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(meanlog,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(sdlog,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  qlinearmix(p, param1=meanlog , param2=sdlog, prop=prop, lowlim=0, npts=1000, envir=cylcop_lnormmix.env, func=plnorm)
}



# Random numbers
#'
#' @rdname gammamix
#' @export
#'
rgammamix <- function(n, shape, rate=1, scale=1/rate, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(rate,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rate <- 1/scale
  rlinearmix(n,param1=shape, param2=rate, prop=prop, func=rgamma)
}


# Density
#'
#' @rdname gammamix
#' @export
#'
dgammamix <- function(x, shape, rate=1, scale=1/rate, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(rate,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rate <- 1/scale
  pdlinearmix(point=x,param1=shape, param2=rate, prop=prop, func=dgamma)
}


# Distribution
#'
#' @rdname gammamix
#' @export
#'
pgammamix <- function(q, shape, rate=1, scale=1/rate, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(q,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(rate,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rate <- 1/scale
  pdlinearmix(point=q,param1=shape, param2=rate, prop=prop, func=pgamma)
}


# Quantiles
#'
#' @rdname gammamix
#' @export
#'
qgammamix <- function(p, shape, rate=1, scale=1/rate, prop){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(shape,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(rate,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(scale,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  rate <- 1/scale
  qlinearmix(p, param1=shape, param2=rate, prop=prop, lowlim=0, npts=1000, envir=cylcop_gammamix.env, func=pgamma)
}


cylcop_gammamix.env <- new.env(parent = emptyenv())
cylcop_lnormmix.env <- new.env(parent = emptyenv())
cylcop_weibullmix.env <- new.env(parent = emptyenv())
cylcop_normmix.env <- new.env(parent = emptyenv())


rlinearmix <- function(n, param1, param2, prop, func){
  if(!is.numeric(prop)||!is.numeric(param1)||!is.numeric(param2)){
    stop(
      "parameter(s), and proportion must be numeric vectors"
    )
  }

  if(length(prop)!=length(param1)||length(prop)!=length(param2)||length(param1)!=length(param2)){
    stop(
      "parameter(s), and proportion must have the same length"
    )
  }

  param1 <- as.vector(param1, mode = "numeric")
  param2 <- as.vector(param2, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)

  dist_component <- sample(1:length(prop), size = n, replace = TRUE, prob = prop)
  draws_per_component <- tabulate(dist_component, nbins = length(prop))
  x <- rep(0, n)
  for (i in seq_along(prop)) {
    x[dist_component == i] <- func(draws_per_component[i], param1[i], param2[i])
  }
  return(x)
}

pdlinearmix <- function(point, param1, param2, prop, func){
  if(!is.numeric(prop)||!is.numeric(param1)||!is.numeric(param2)){
    stop(
      "parameter(s), and proportion must be numeric vectors"
    )
  }

  if(length(prop)!=length(param1)||length(prop)!=length(param2)||length(param1)!=length(param2)){
    stop(
      "parameter(s), and proportion must have the same length"
    )
  }
  param1 <- as.vector(param1, mode = "numeric")
  param2 <- as.vector(param2, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)

  out <- 0
  for (i in seq_along(prop)) {
    out <- out + prop[i]*func(point,param1[i],param2[i])
  }
  out
}

qlinearmix <- function(p, param1, param2, prop, lowlim, npts, envir, func){
  if(!is.numeric(prop)||!is.numeric(param1)||!is.numeric(param2)){
    stop(
      "parameter(s), and proportion must be numeric vectors"
    )
  }

  if(length(prop)!=length(param1)||length(prop)!=length(param2)||length(param1)!=length(param2)){
    stop(
      "parameter(s), and proportion must have the same length"
    )
  }
  param1 <- as.vector(param1, mode = "numeric")
  param2 <- as.vector(param2, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)


  if(all(p>0.999999)) return(rep(Inf,length(p)))
  if(all(p<0.000001)) return(rep(lowlim,length(p)))

  #Interpolation with splines. When the qnormmix() is called
  #for the first time after loading the package, the cdf values are calculated and
  #stored in a  new environment. This takes some time. All subsequent calls to qnormmix()
  #use the stored cdf values and are very fast.

  if(lowlim>=0){
    min_val1 <- 0
    min_val2 <- 0
  }else{
    curr <- Inf
    min_val1 <- -0.1
    min_val2 <- -0.1
    while(curr>0.000001){
      min_val1 <- min_val1*2
      curr <- pdlinearmix(min_val1, param1, param2, prop, func)
      if(curr>0.01) min_val2 <- min_val1
    }
    if(curr==0){
      incr <- -0.01*min_val1
      while(curr<0.000001){
        min_val1 <- min_val1+incr
        curr <- pdlinearmix(min_val1, param1, param2, prop, func)
      }
    }
  }

  max_val1 <- 0.1
  max_val2 <- 0.1
  curr <- -Inf
  while(curr<0.999999){
    max_val1 <- max_val1*2
    curr <- pdlinearmix(max_val1, param1, param2, prop, func)
    if(curr<0.99) max_val2 <- max_val1
  }

  if(curr==1){
    decr <- 0.01*max_val1
    while(1-curr<0.000001){
      max_val1 <- max_val1-decr
      curr <- pdlinearmix(max_val1, param1, param2, prop, func)
    }
  }

  if(min_val1==min_val2){
    pts <- unique(c(seq(min_val1, max_val2, length.out = round(0.95*npts)),
                    seq(max_val2, max_val1, length.out = round(0.05*npts))))
  }else if(max_val1==max_val2){
    pts <- unique(c(seq(min_val1, min_val2, length.out = round(0.05*npts)),
                    seq(min_val2, max_val2, length.out = round(0.95*npts))))
  }else{
    pts <- unique(c(seq(min_val1, min_val2, length.out = round(0.025*npts)),
                    seq(min_val2, max_val2, length.out = round(0.95*npts)),
                    seq(max_val2, max_val1, length.out = round(0.025*npts))))
  }



  interpolate_from_prev_calc <- function(param1, param2, prop){
    if(any(c(param1, param2, prop)!=envir$parameter_vals)){
      stop("not the parameters in the environment")
    }
    result <- stats::spline(envir$cdf, pts, method = "hyman", xout = p)$y
    return(result)
  }

  result <- tryCatch(interpolate_from_prev_calc(param1, param2, prop),
                     error=function(e){
                       #cylcop.env$cdf doesn't exist yet, i.e. first time calling the function
                       #Or it corresponds to other parameter values than the ones the function is currently called with
                       assign("parameter_vals",c(param1, param2, prop),envir=envir)
                       assign("cdf",
                              pdlinearmix(pts, param1, param2, prop, func),
                              envir=envir)
                       result <- stats::spline(envir$cdf, pts, method = "hyman", xout = p)$y
                       return(result)
                     })
  result[which(p>0.999999)] <- Inf
  result[which(p<0.000001)] <- lowlim
  return(result)
}

