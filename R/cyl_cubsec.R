#' @include cyl_cop_class.R
NULL



#' An S4 class of copulas with cubic sections
#'
#' This class contains bivariate circular-linear copulas with cubic sections in the linear dimension.
#' They are periodic in the circular dimension, u, and symmetric with respect to u=0.5. I.e. there is
#' symmetry between positive and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by \code{cyl_cubsec()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
#'
#' @section Extends:
#' Class \code{cyl_cubsec} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass("cyl_cubsec", contains = "cyl_copula")



#' Construction of \code{cyl_cubsec} objects
#'
#' @param a The numeric value of the first parameter of the copula. It must be in [- 1 / (2 * pi)), 1 / (2 * pi))]
#' @param b The numeric value of the second parameter of the copula. It must be in [- 1 / (2 * pi)), 1 / (2 * pi))]
#'
#' @export
#'
#' @examples
#' cyl_cubsec(a = 0.1, b = -0.1)
#'
cyl_cubsec <- function(a = 1 / (2 * pi),
                       b = 1 / (2 * pi)) {


  lowbnd = c(-1 / (2 * pi),-1 / (2 * pi))
  upbnd = c(1 / (2 * pi), 1 / (2 * pi))

  new(
    "cyl_cubsec",
    name = "Cub. sect. copula",
    parameters = c(a, b),
    param.names = c("a", "b"),
    param.lowbnd = lowbnd,
    param.upbnd = upbnd
  )
}



#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_cubsec"), function(n, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  cop_uv <- matrix(nrow = n,
                   ncol = 2,
                   dimnames = list(c(), c("u", "v")))
  for (i in 1:n) {
    u <- runif(1)
    w <- runif(1)
    alpha_prime <- a * 2 * pi * cos(2 * pi * u)
    beta_prime <- -b * 2 * pi * cos(2 * pi * u)

    #It is easier, and not really slower, to find the inverse of the Conditional distribution of V given u
    #numerically instead of analytically

    # alpha_prime <- complex(alpha_prime)
    # beta_prime <- complex(beta_prime)
    # cdf_cond_inv <- function(v, alpha_prime, beta_prime){
    #   v <- complex(v)
    #   f1 <- -18*alpha_prime^2-2*alpha_prime^3+
    #     27*alpha_prime*beta_prime+3*alpha_prime^2*beta_prime-
    #     9*beta_prime^2+3*alpha_prime*beta_prime^2-
    #     2*beta_prime^3+27*alpha_prime^2*v-54*alpha_prime*beta_prime*v+
    #     27*beta_prime^2*v
    #   f2 <- (f1+sqrt(4*(3*alpha_prime-alpha_prime^2-3*beta_prime+
    #                       alpha_prime*beta_prime-beta_prime^2)^3+f1^2))^(1/3)
    #
    #   -((-2*alpha_prime+beta_prime)/
    #       (3*(alpha_prime-beta_prime)))-
    #     (2^(1/3)*(3*alpha_prime-alpha_prime^2-3*beta_prime+alpha_prime*beta_prime-beta_prime^2))/
    #     (3*(alpha_prime-beta_prime)*f2)+(1/(3*2^(1/3)*(alpha_prime-beta_prime)))*f2
    # }

    #Calculate numerically the inverse of the conditional distribution of V given u, C_u(v) and
    #evaluate it at w

    #C_u(v)
    cdf_cond <- function(v) {
      v + v * alpha_prime + v ^ 2 * (beta_prime - 2 * alpha_prime) + v ^ 3 * (alpha_prime -
                                                                                beta_prime)
    }
    #Invert it
    cdf_cond_inv <- GoFKernel::inverse(f = cdf_cond,
                            lower = 0,
                            upper = 1)
    #evaluate it at w
    v <- cdf_cond_inv(w)
    cop_uv[i, 1] <- u
    cop_uv[i, 2] <- cdf_cond_inv(w)
  }
  return(cop_uv)
})



#' Calcualte density
#' @rdname Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_cubsec"), function(u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  alpha_prime <- function(a, u) {
    a * 2 * pi * cos(2 * pi * u)
  }
  beta_prime <- function(b, u) {
    -b * 2 * pi * cos(2 * pi * u)
  }
  pdf <-
    1 + alpha_prime(a, u) + 2 * v * (beta_prime(b, u) - 2 * alpha_prime(a, u)) + 3 *
    v * v * (-beta_prime(b, u) + alpha_prime(a, u))
  return(c(pdf))
})



#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_cubsec"), function(u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  alpha <- function(a, u) {
    a * sin(2 * pi * u)
  }
  beta <- function(b, u) {
    -b * sin(2 * pi * u)
  }
  cdf <-
    u * v + v * (1 - v) * (alpha(a, u) * (1 - v) + beta(b, u) * v)
  return(c(cdf))
})



#-----Change attributes of existing cyl_cubsec object.-------------------------------------------
#
#' @rdname setCopParam
#' @export
setMethod("setCopParam", "cyl_cubsec", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
