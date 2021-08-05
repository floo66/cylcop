#' @include cyl_cop_class.R
NULL



#' An S4 Class of Bivariate Copulas with Quadratic Sections
#'
#' This class contains bivariate circular-linear copulas with quadratic sections
#'  in the linear dimension. They are periodic in the circular dimension, u,
#'  and symmetric with respect to u=0.5. I.e the can capture correlation in data
#'  where there is symmetry between positive and negative angles. These copulas
#'  are described by one parameter, \code{a}.
#'
#' @section Objects from the Class:
#' Objects are created by \code{\link{cyl_quadsec}()}.
#'
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the
#' parameter value.
#' @slot param.names \link[base]{character} \link[base]{vector} holding the
#' parameter name.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower
#'  bound of the parameter.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper
#'  bound of the parameter.
#'
#' @section Extends:
#' Class '\code{cyl_quadsec}' extends class '\code{\linkS4class{cyl_copula}}'.
#'
#' @references \insertRef{Quesada-Molina1995}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
#'
setClass("cyl_quadsec", contains = "cyl_copula")



#' Construction of '\code{cyl_quadsec}' Objects
#'
#' Constructs a circular-linear copula with cubic sections of class
#'  '\code{\linkS4class{cyl_quadsec}}'.
#'
#' @param a \link[base]{numeric} value of the parameter of the copula. It must be in
#' \eqn{[- 1 / (2 \pi)), 1 / (2 \pi))]}.
#'
#' @return An \R object of class '\code{\linkS4class{cyl_quadsec}}'.
#'
#' @export
#'
#' @examples
#' cop <- cyl_quadsec(a = 0.1)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot")
#'
#' @references \insertRef{Quesada-Molina1995}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
cyl_quadsec <- function(a = 1 / (2 * pi)) {

  lowbnd = -1 / (2 * pi)
  upbnd = 1 / (2 * pi)

  new(
    "cyl_quadsec",
    name = "Quad. sect. copula",
    parameters = a,
    param.names = "a",
    param.lowbnd = lowbnd,
    param.upbnd = upbnd
  )
}



#' Generate Random Samples
#' @rdname Cylcop
# @describeIn cyl_quadsec-class Generate random samples.
#' @export
setMethod("rcylcop", signature("numeric", "cyl_quadsec"), function(n, copula) {
  a <- copula@parameters[1]
  u <- runif(n)
  w <- runif(n)
  #if a==0, u and v are independent, but you would divide by zero in the conditional quantile function
  # so we have to set it for this case here separately
  if (isTRUE(all.equal(a,0))) {
    cop_uv <- cbind(u, w)
  }
  #Calcualte the inverse of the conditional distribution of V given u, C_u(v) and
  #evaluate it at w
  else{
    mat <- matrix(ncol=2,c(u,w))
    v <- cylcop::ccylcop(mat,copula,cond_on=1, inverse=T)
    cop_uv <- cbind(u, v)
  }
  return(cop_uv)
})



#' Calculate Density
#' @rdname Cylcop
# @describeIn cyl_quadsec-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_quadsec"), function(u, copula) {

  a <- copula@parameters[1]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  alpha_prime <- function(a, u) {
    a * 2 * pi * cos(2 * pi * u)
  }
  pdf <- 1 + alpha_prime(a, u) * (1 - 2 * v)
  return(c(pdf))
})



#' Calculate Distribution
#' @rdname Cylcop
# @describeIn cyl_quadsec-class Calculate the distribution.
#' @export
setMethod("pcylcop", signature("matrix", "cyl_quadsec"), function(u, copula) {
  a <- copula@parameters[1]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  alpha <- function(a, u) {
    a * sin(2 * pi * u)
  }
  cdf <- u * v + alpha(a, u) * (v - (v * v))
  return(c(cdf))
})

#' Condtional copula
#' @rdname ccylcop
# @describeIn cyl_quadsec-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_quadsec"), function(u,
                                                        copula,
                                                        cond_on = 2,
                                                        inverse = FALSE) {
  a <- copula@parameters[1]
  u_orig <- matrix(ncol=2,u)
  length <- nrow(u)
  v <- u_orig[, 2, drop = F]
  u <- u_orig[, 1, drop = F]

  if(cond_on==2){
    alpha <- a * sin(2 * pi * u)
    if(inverse==F){
      result <- u+alpha-2*v*alpha
    }
    if(inverse==T){
      result <-  numerical_inv_conditional_cop(u_orig, copula, cond_on=2)
    }
  }
  else if(cond_on==1){
    alpha_prime <- a * 2 * pi * cos(2 * pi * u)
    if(inverse==F){
      result<- v+alpha_prime*(v-v^2)
    }
    if(inverse==T){
      result<- (alpha_prime + 1 - sqrt((alpha_prime + 1) ^ 2 - 4 * alpha_prime * v)) /
        (2 * alpha_prime)
    }
  }
  else stop("cond_on must be either 1 or 2")
  return(result%>%as.numeric())
})


#-----Change attributes of existing cyl_quadsec object.-------------------------------------------
#
#' @rdname setCopParam
# @describeIn cyl_quadsec-class Change attributes of existing object.
#' @export
setMethod("setCopParam", "cyl_quadsec", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
