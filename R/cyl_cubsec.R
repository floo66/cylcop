#' @include cyl_cop_class.R
NULL



#' An S4 Class of Bivariate Copulas with Cubic Sections
#'
#' This class contains bivariate circular-linear copulas with cubic sections in the linear dimension.
#' They are periodic in the circular dimension, u, and symmetric with respect to u=0.5. I.e.
#' the can capture correlation in data where there is symmetry between positive and negative angles.
#' These copulas are described by two parameters, \code{a} and \code{b}.
#'
#' @section Objects from the Class:
#' Objects are created by \code{\link{cyl_cubsec}()}.
#'
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the parameter values.
#' @slot param.names \link[base]{character} \link[base]{vector} holding the
#' parameter names.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower
#'  bounds of the parameters.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper
#'  bounds of the parameters.
#'
#' @section Extends:
#' Class '\code{cyl_cubsec}' extends class '\code{\linkS4class{cyl_copula}}'.
#'
#' @references \insertRef{Nelsen1997}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
#'
setClass("cyl_cubsec", contains = "cyl_copula")



#' Construction of '\code{cyl_cubsec}' Objects
#'
#' Constructs a circular-linear copula with cubic sections of class
#'  '\code{\linkS4class{cyl_cubsec}}'.
#'
#' @param a \link[base]{numeric} value of the first parameter of the copula.
#' It must be in \eqn{[- 1 / (2 \pi)), 1 / (2 \pi))]}.
#' @param b \link[base]{numeric} value of the second parameter of the copula.
#' It must be in \eqn{[- 1 / (2 \pi)), 1 / (2 \pi))]}.
#'
#' @return An \R object of class '\code{\linkS4class{cyl_cubsec}}'.
#'
#' @export
#'
#' @examples
#' cop <- cyl_cubsec(a = 0.1, b = -0.1)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot")
#'
#' @references \insertRef{Nelsen1997}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#'
cyl_cubsec <- function(a = 1 / (2 * pi),
                       b = 1 / (2 * pi)) {
  lowbnd = c(-1 / (2 * pi), -1 / (2 * pi))
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
#' @rdname Cylcop
# @describeIn cyl_cubsec-class Generate random samples.
#' @export
setMethod("rcylcop", signature("numeric", "cyl_cubsec"), function(n, copula) {
  u <- runif(n)
  w <- runif(n)
  mat <- matrix(ncol = 2, c(u, w))
  v <- cylcop::ccylcop(mat, copula, cond_on = 1, inverse = TRUE)
  cop_uv <- cbind(u, v)
  return(cop_uv)
})



#' Calcualte density
#' @rdname Cylcop
# @describeIn cyl_cubsec-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_cubsec"), function(u, copula) {
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
#' @rdname Cylcop
# @describeIn cyl_cubsec-class Calculate the distribution.
#' @export
setMethod("pcylcop", signature("matrix", "cyl_cubsec"), function(u, copula) {
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



#' Condtional copula
#' @rdname ccylcop
# @describeIn cyl_cubsec-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_cubsec"), function(u,
                                                       copula,
                                                       cond_on = 2,
                                                       inverse = FALSE) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  u_orig <- matrix(ncol = 2, u)
  length <- nrow(u)
  v <- u_orig[, 2, drop = F]
  u <- u_orig[, 1, drop = F]

  if (cond_on == 2) {
    alpha <- a * sin(2 * pi * u)
    beta <- -b * sin(2 * pi * u)
    if (inverse == F) {
      result <-
        u + (-alpha + beta) * (1 - v) * v + (1 - v) * (alpha * (1 - v) + beta *
                                                         v) - v * (alpha * (1 - v) + beta * v)
    }
    if (inverse == TRUE) {
      result <-  numerical_inv_conditional_cop(u_orig, copula, cond_on = 2)
    }
  }
  else if (cond_on == 1) {
    if (inverse == F) {
      alpha_prime <- a * 2 * pi * cos(2 * pi * u)
      beta_prime <- -b * 2 * pi * cos(2 * pi * u)
      result <-
        v + v * alpha_prime + v ^ 2 * (beta_prime - 2 * alpha_prime) +
        v ^ 3 * (alpha_prime - beta_prime)
    }
    if (inverse == TRUE) {
      if (abs(a + b) < 0.001) {
        #numerical instability in the analytical solution, use numeric approximation instead
        result <-
          numerical_inv_conditional_cop(u_orig, copula, cond_on = 1)
      } else{
        result <- map2_dbl(u, v, function(u, w) {
          alpha_prime <- a * 2 * pi * cos(2 * pi * u)
          beta_prime <- -b * 2 * pi * cos(2 * pi * u)
          if (abs(w - 1) < 0.0000001)
            return(1.0)
          if (abs(w) < 0.0000001)
            return(0.0)
          term1 <-
            -18 * alpha_prime ^ 2 - 2 * alpha_prime ^ 3 + 27 * alpha_prime * beta_prime + 3 *
            alpha_prime ^ 2 * beta_prime - 9 * beta_prime ^ 2 + 3 * alpha_prime * beta_prime ^
            2 - 2 * beta_prime ^ 3 + 27 * alpha_prime ^ 2 * w - 54 * alpha_prime * beta_prime *
            w + 27 * beta_prime ^ 2 * w
          term2 <-
            sqrt(as.complex(
              4 * (
                3 * alpha_prime - alpha_prime ^ 2 - 3 * beta_prime + alpha_prime * beta_prime - beta_prime ^
                  2
              ) ^ 3 +
                (
                  -18 * alpha_prime ^ 2 - 2 * alpha_prime ^ 3 + 27 * alpha_prime * beta_prime + 3 *
                    alpha_prime ^ 2 * beta_prime -
                    9 * beta_prime ^ 2 + 3 * alpha_prime *
                    beta_prime ^ 2 - 2 * beta_prime ^ 3 +
                    27 * alpha_prime ^ 2 * w - 54 * alpha_prime *
                    beta_prime * w + 27 * beta_prime ^ 2 * w
                ) ^ 2
            ))

          if (abs(Re(term1 + term2)) < 10 ^ (-11)) {
            #avoid numerical instability
            mat <- matrix(ncol = 2, c(u, w))
            solution <-
              numerical_inv_conditional_cop(mat, copula, cond_on = 1)
            return(solution)
          }

          term3 <- (term1 + term2) ^ (1 / 3)

          solution1 <-
            -((beta_prime - 2 * alpha_prime) / (3 * (alpha_prime - beta_prime))) -
            (complex(real = 1, imaginary = sqrt(3)) * term3 / (6 * 2 ^ (1 /
                                                                          3) * (alpha_prime - beta_prime))) +
            (
              complex(real = 1, imaginary = -sqrt(3)) * (
                3 * alpha_prime - alpha_prime ^ 2 - 3 * beta_prime + alpha_prime * beta_prime - beta_prime ^
                  2
              ) / (3 * 2 ^ (2 / 3) * (alpha_prime - beta_prime) * term3)
            )

          solution2 <-
            -((beta_prime - 2 * alpha_prime) / (3 * (alpha_prime - beta_prime))) -
            (complex(real = 1, imaginary = -sqrt(3)) * term3 / (6 * 2 ^ (1 /
                                                                           3) * (alpha_prime - beta_prime))) +
            (
              complex(real = 1, imaginary = sqrt(3)) * (
                3 * alpha_prime - alpha_prime ^ 2 - 3 * beta_prime + alpha_prime * beta_prime - beta_prime ^
                  2
              ) / (3 * 2 ^ (2 / 3) * (alpha_prime - beta_prime) * term3)
            )

          solution3 <-
            -((beta_prime - 2 * alpha_prime) / (3 * (alpha_prime - beta_prime))) +
            (term3 / (3 * 2 ^ (1 / 3) * (alpha_prime - beta_prime))) -
            ((
              2 ^ (1 / 3) * (
                3 * alpha_prime - alpha_prime ^ 2 - 3 * beta_prime + alpha_prime * beta_prime - beta_prime ^
                  2
              )
            ) / (3 * (alpha_prime - beta_prime) * term3))
          solution <- c(solution1, solution2, solution3)
          real <- Re(solution)
          imaginary <- abs(Im(solution))
          ind <- which(real >= 0 & real <= 1)
          ind2 <- which(imaginary[ind] == min(imaginary[ind]))
          return(real[ind][ind2])
        })
      }

    }
  }
  else
    stop("cond_on must be either 1 or 2")
  return(result %>% as.numeric())
})



#-----Change attributes of existing cyl_cubsec object.-------------------------------------------
#
#'
#' @rdname setCopParam
# @describeIn cyl_cubsec-class Change attributes of existing object.
#' @export
setMethod("setCopParam", "cyl_cubsec", function(copula, param_val, param_name) {
  if (is.null(param_name))
    param_name <- copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
