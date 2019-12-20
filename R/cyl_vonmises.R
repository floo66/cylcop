#' @include cyl_cop_class.R
NULL


#' An S4 class of bivariate vonMises copulas
#'
#' This class contains circular-linear copulas that are based on the approach by Johnson and
#' Wehrly 1978 with a vonMises periodic function. They are periodic in the circular
#' dimension, u, but not symmetric with respect to u=0.5. I.e. there is no
#' symmetry between positive and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by \code{cyl_vonmises()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
#' @slot flip A logical indicating whether the copula should be rotated 90 degrees to have negative correlation.
#'
#' @section Extends:
#' Class \code{cyl_vonmises} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass("cyl_vonmises", contains = "cyl_copula", slots = "flip")




#' Construction of \code{cyl_vonmises} objects
#'
#' @param mu A numeric value giving the mean of the vonMises function used to construct the copula.
#' @param kappa A numeric value giving the concentration of the vonMises function used to construct the copula.
#' @param flip A logical indicating whether the copula should be rotated 90 degrees to have negative correlation.
#'
#' @export
#'
#' @examples
#' cyl_vonmises(mu=pi, kappa=10, flip = TRUE)
#'
cyl_vonmises <- function(mu = 0,
                         kappa = 1,
                         flip = FALSE) {

  lowbnd = c(-Inf, 0)
  upbnd = c(Inf, Inf)

  new(
    "cyl_vonmises",
    name = "vonMises copula",
    parameters = c(mu, kappa),
    param.names = c("mu", "kappa"),
    param.lowbnd = lowbnd,
    param.upbnd = upbnd,
    flip = flip
  )
}



#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_vonmises"), function(n, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  u <- runif(n)
  w <- runif(n)
  #Calcualte the inverse of the conditional distribution of V given u, C_u(v) and
  #evaluate it at w
  #this is the inverse of pvonmises(v, mu=2*pi*u), which is the integrand in the copula cdf function below,
  #(1/(2*pi))*qvonmises(pvonmises(v, mu=2*pi*u),mu=2*pi*u)=v
  v <- (1 / (2 * pi)) *
    map2_dbl(u, w, ~ circular::qvonmises(.y,
                               mu = circular::circular(2 * pi * .x - mu),
                               kappa = kappa,
                               from = circular::circular(0)
             ))
  if (copula@flip)
    u <- 1 - u
  cop_uv <- cbind(u, v)
  return(cop_uv)
})



#' Calcualte density
#' @rdname Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_vonmises"), function(u, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  #drop=F, so it aslo works with single numbers
  v <- u[, 2, drop = F]
  if (copula@flip)
    u <- 1 - u[, 1, drop = F]
  else
    u <- u[, 1, drop = F]
  #Johnson and Wehrly's 1978 equation
  circular::dvonmises(x = circular::circular(2 * pi * (u - v)),
            mu = circular::circular(mu),
            kappa = kappa) %>%
    as.double() * 2 * pi
})



#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_vonmises"), function(u, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  #Slow version with numerical integration
  # integrand <- Vectorize(function(x, v, mu, kappa) pvonmises(1.99999999*pi*v-x, mu = mu, kappa = kappa, from = -x)) #1.99999 instead of 2, because function is defined on
  # out <- integrate(f = integrand, lower = 0, upper = 2*pi*u, v = v, mu = mu, kappa = kappa) #works only if f accepts vecor input
  # (1/(2*pi))*out$value

  #Fast version with sum of Bessel functions
  ppvm <- function(x, kappa) {
    tol <- 1e-20
    flag <- TRUE
    j <- 1
    sum <- 0
    while (flag) {
      term <- (besselI(x = kappa, nu = j,expon.scaled = FALSE)
               * cos(j * x)) / (j * j)
      sum <- sum + term
      j <- j + 1
      if (all(abs(term) < tol))
        flag <- FALSE
    }
    return(sum / (pi * besselI(x = kappa, nu = 0, expon.scaled = FALSE)))
  }

  if (!copula@flip) {
    cdf <- 1 / (2 * pi) *
      (
        v * 2 * pi * u + ppvm(2 * pi * (u - v) - mu, kappa) -
          ppvm(2 * pi * u - mu, kappa) -
          ppvm(-2 * pi * v - mu, kappa) +
          ppvm(-mu, kappa)
      )
  }
  else{
    #if the copula is flipped, instead of integrating from (0,0) to (u,v) we need to integrate from  (1-u,0) to (1,v).
    #Use the C-volume of the unflipped copula for that, it recursively calls this pCopula function.
    unflipped <- copula
    unflipped@flip <- FALSE
    cdf <- prob(unflipped, l = c((1 - u), 0), u = c(1, v))
  }
  return(c(cdf))
})



#-----Change attributes of existing cyl_vonmises object.-------------------------------------------
#
#' @rdname setCopParam
#' @export
setMethod("setCopParam", "cyl_vonmises", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
