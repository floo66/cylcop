#' @include cyl_cop_class.R
NULL


#'Numerically calculate the conditional copula
#'
#' @param u Matrix of numeric values in I^2, containing as first collumn
#'  the circular (periodic) and as second the linear dimension.
#' @param copula A \code{cyl_copula} or 2-dimensional \code{copula} object.
#' @param cond_on Column number of u on which the copula is conditioned. E.g if
#' \code{cond_on=2}, the function calculates for each element in the first column of u
#' the Copula conditional on the element in the second column.
#'
#' @return
#' A vector containing the values of the distribution of the copula at
#' \code{[u,-cond_on]} conditional on the values of \code{[u,cond_on]}
#' @export
#'
numerical_conditional_cop <- function(u, copula, cond_on){
  u <- matrix(ncol=2,u)
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  if(cond_on==2){
    Ccond <- map2_dbl(u, v, function(u, v) {
      integrand <- function(x) {
        cylcop::dCopula(c(x[1], v), copula)
      }
      cond <-
        stats::integrate(f = Vectorize(integrand),
                         lower = 0,
                         upper = u)
      return(cond$value)
    })
  }
  else if(cond_on==1){
    Ccond <- map2_dbl(u, v, function(u, v) {
      integrand <- function(y) {
        cylcop::dCopula(c(u, y[1]), copula)
      }
      cond <-
        stats::integrate(f = Vectorize(integrand),
                         lower = 0,
                         upper = v)
      return(cond$value)
    })
  }
  else stop("cond_on for conditional distribution must be 1 or 2.")
  return(Ccond)
}

#'Numerically calculate the inverse of the conditional copula
#'
#' @param u Matrix of numeric values in I^2, containing as first collumn
#'  the circular (periodic) and as second the linear dimension.
#' @param copula A \code{cyl_copula} or 2-dimensional \code{copula} object.
#' @param cond_on Column number of u on which the copula is conditioned. E.g if
#' \code{cond_on=2}, the function calculates for each element in the first column of u
#' the inverse of the Copula conditional on the element in the second column.
#'
#' @return
#' A vector containing the values of the inverse distribution of the copula at
#' \code{[u,-cond_on]} conditional on the values of \code{[u,cond_on]}
#' @export
#'
numerical_inv_conditional_cop <- function(u, copula, cond_on){
  u <- matrix(ncol=2,u)
  length <- nrow(u)
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  Ccond<-rep(1,length)


  if(cond_on==2){
    for (i in 1:length) {
      cond_func_v <- function(u){
        integrand <- function(x) {
          cylcop::dCopula(c(x[1], v[i,]), copula)
        }
        cond <-
          stats::integrate(f = Vectorize(integrand),
                           lower = 0,
                           upper = u)
        return(cond$value)
      }

      cdf_cond_inv <- GoFKernel::inverse(f = cond_func_v,
                                         lower = 0,
                                         upper = 1)
      Ccond[i] <- cdf_cond_inv(u[i,])
    }
  }
  else if(cond_on==1){
    for (i in 1:length) {
      cond_func_u <- function(v){
        integrand <- function(y) {
          cylcop::dCopula(c(u[i,], y[1]), copula)
        }
        cond <-
          stats::integrate(f = Vectorize(integrand),
                           lower = 0,
                           upper = v)
        return(cond$value)
      }

      cdf_cond_inv <- GoFKernel::inverse(f = cond_func_u,
                                         lower = 0,
                                         upper = 1)
      Ccond[i] <- cdf_cond_inv(v[i,])
    }
  }
  else stop("cond_on for conditional distribution must be 1 or 2.")
  return(Ccond)
}


#-----Calculate the conditional copula of a copula-object.-------------------------------------------
#
#' @rdname cCopula
#' @export
setMethod("cCopula", "Copula", function(u, copula, cond_on, inverse) {
  if(dim(copula)!=2) stop("cylcop::cCopula() works only for copulas of dimension 2. For other copula-objects, try copula::cCopula().")

  if(cond_on==2){
    u_swap <- matrix(ncol=2,c(u[,2],u[,1]))
    tryCatch({
      copula::cCopula(u_swap,copula,indices =2, inverse = inverse) %>% as.numeric()
    }, error = function(e) {
      if(inverse==F){
      numerical_conditional_cop(u,copula,cond_on = 2)
      }
      else{
        numerical_inv_conditional_cop(u,copula,cond_on = 2)
      }
    })
  }
  else if(cond_on==1){
    tryCatch({
      copula::cCopula(u,copula,indices =2, inverse = inverse) %>% as.numeric()
    }, error = function(e) {
      if(inverse==F){
        numerical_conditional_cop(u,copula,cond_on = 1)
      }
      else{
        numerical_inv_conditional_cop(u,copula,cond_on = 1)
      }
    })
  }
  else stop("cond_on must be either 1 or 2")
})
