#' @include cyl_cop_class.R
NULL


#' An S4 class of bivariate Gaussian copulae
#'
#' Could use \code{normalCopula} from \code{copula}-package instead, just here for legacy
#' reasons, consistency in output format etc. It is not in any way periodic, just the
#' usual Gaussian copula.
#'
#' @section Objects from the Class:
#' Objects are created by \code{cyl_gauss()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
#'
#' @section Extends:
#' Class \code{cyl_gauss} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass("cyl_gauss", contains = "cyl_copula")

#' Construction of \code{cyl_gauss} objects
#'
#' @param rho A numeric value giving the correlation of the copula
#'
#' @export
#'
#' @examples
#' cyl_gauss(rho=0.5)
cyl_gauss <- function(rho = 0.9) {
  lowbnd = -1
  upbnd = 1
  if(rho < -1 || rho > 1) stop(cylcop::error_sound(), "parameter out of bounds")

  new(
    "cyl_gauss",
    name = "Gaussian copula",
    parameters = rho,
    param.names = "rho",
    param.lowbnd = lowbnd,
    param.upbnd = upbnd
  )
}


#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_gauss"), function(n, copula) {
  rho <- copula@parameters[1]
  sigma <- matrix(c(1, rho, rho, 1), 2)
  normdraw <- mvtnorm::mvrnorm(n = n,
                      mu = c(0, 0),
                      Sigma = sigma)
  unif <- pnorm(normdraw)
  return(unif)
})

#' Calcualte density
#' @rdname  Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_gauss"), function(u, copula) {
  rho <- copula@parameters[1]
  v <- u[2]
  u <- u[1]
  sigma <- matrix(c(1, rho, rho, 1), 2)
  marg <- c(qnorm(u, mean = 0, sd = 1), qnorm(v, mean = 0, sd = 1))
  # dmvnorm(x = marg, mean = c(0,0), sigma = sigma) /
  #  (dnorm(marg[1]) * dnorm(marg[2]))
  # Division and multiplication of small doubles, convert to log, avoid under-/ overflow
  (mvtnorm::dmvnorm(
    x = marg,
    mean = c(0, 0),
    sigma = sigma,
    log = TRUE
  ) -
      (dnorm(marg[1], log = TRUE) + dnorm(marg[2], log = TRUE))) %>%
    exp()
})

#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_gauss"), function(u, copula) {
  rho <- copula@parameters[1]
  v <- u[2]
  u <- u[1]
  sigma <- matrix(c(1, rho, rho, 1), 2)
  marg <- c(qnorm(u, mean = 0, sd = 1), qnorm(v, mean = 0, sd = 1))
  as.vector(mvtnorm::pmvnorm(
    upper = marg,
    mean = c(0, 0),
    sigma = sigma
  ))
})
