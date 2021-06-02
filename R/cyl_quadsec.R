#' @include cyl_cop_class.R
NULL



#' An S4 class of copulae with quadratic sections
#'
#' This class contains bivariate circular-linear copulae with quadratic sections in the linear dimension.
#' They are periodic in the circular dimension, u, and symmetric with respect to u=0.5. I.e. there is
#' symmetry between positive and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by \code{cyl_quadsec()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
#'
#' @section Extends:
#' Class \code{cyl_quadsec} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass("cyl_quadsec", contains = "cyl_copula")



#' Construction of \code{cyl_quadsec} objects
#'
#' @param a The numeric value of the parameter of the copula. It must be in [- 1 / (2 * pi)), 1 / (2 * pi))]
#'
#' @export
#'
#' @examples
#' cyl_quadsec(a = 0.1)
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



#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_quadsec"), function(n, copula) {
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
    v <- cylcop::cCopula(mat,copula,cond_on=1, inverse=T)
    cop_uv <- cbind(u, v)
  }
  return(cop_uv)
})



#' Calcualte density
#' @rdname Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_quadsec"), function(u, copula) {

  a <- copula@parameters[1]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  alpha_prime <- function(a, u) {
    a * 2 * pi * cos(2 * pi * u)
  }
  pdf <- 1 + alpha_prime(a, u) * (1 - 2 * v)
  return(c(pdf))
})



#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_quadsec"), function(u, copula) {
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
#' @rdname cCopula
#' @export
setMethod("cCopula", signature("cyl_quadsec"), function(u, copula, cond_on=2, inverse=F) {
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
#' @export
setMethod("setCopParam", "cyl_quadsec", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
