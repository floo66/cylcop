#' @include cyl_cop_class.R
NULL


#' Numerically Calculate the Conditional Copula
#'
#' @param u \link[base]{matrix} or \link[base]{vector} of \link[base]{numeric}
#' values in \eqn{I^2}, containing as first column
#'  the circular (periodic) and as second the linear dimension.
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param cond_on column number of \code{u} on which the copula is conditioned. E.g. if
#' \code{cond_on = 2}, the function calculates for each element in the first column of u
#' the copula conditional on the element in the second column.
#'
#' @return A vector containing the values of the distribution of the copula at
#' \code{u[,-cond_on]} conditional on the values of \code{u[,cond_on]}.
#'
#' @examples cop <- cyl_quadsec(0.1)
#' u <- cbind(c(0.3, 0.1), c(0.7, 0.3))
#' numerical_conditional_cop(u = u, cop = cop, cond_on = 1)
#'
#' @references \insertRef{Hodelmethod}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' @seealso \code{\link{ccylcop}()}, \code{\link{numerical_inv_conditional_cop}()}.
#'
#' @export
#'
numerical_conditional_cop <- function(u, copula, cond_on){
  #validate input
  tryCatch({
    check_arg_all(list(check_argument_type(u, type="numeric",
                                           length=2,
                                           lower=0,
                                           upper=1),
                       check_argument_type(u, type="matrix",
                                           ncol=2,
                                           lower=0,
                                           upper=1))
                  ,2)
    check_arg_all(list(check_argument_type(copula, type="cyl_copula"),
                       check_argument_type(copula, type="Copula"))
                  ,2)
    check_arg_all(check_argument_type(cond_on,
                                      type="numeric",
                                      values = c(1,2))
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  u <- matrix(ncol=2,u)
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  if(cond_on==2){
    Ccond <- map2_dbl(u, v, function(u, v) {
      integrand <- function(x) {
        dcylcop(c(x[1], v), copula)
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
        dcylcop(c(u, y[1]), copula)
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
#' @param u \link[base]{matrix} or \link[base]{vector} of \link[base]{numeric}
#' values in \eqn{I^2}, containing as first column
#'  the circular (periodic) and as second the linear dimension.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param cond_on column number of \code{u} on which the copula is conditioned. E.g if
#' \code{cond_on = 2}, the function calculates for each element in the first column
#' of u the inverse of the Copula conditional on the element in the second column.
#'
#' @return
#' A vector containing the values of the inverse distribution of the copula at
#' \code{[u,-cond_on]} conditional on the values of \code{[u,cond_on]}.
#'
#' @examples cop <- cyl_quadsec(0.1)
#' u <- cbind(c(0.3, 0.1), c(0.7, 0.3))
#' numerical_inv_conditional_cop(u = u, cop = cop, cond_on = 1)
#'
#' @references \insertRef{Hodelmethod}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' @seealso \code{\link{ccylcop}()}, \code{\link{numerical_conditional_cop}()}.
#'
#' @export
#'
numerical_inv_conditional_cop <- function(u, copula, cond_on){

  tryCatch({
    check_arg_all(list(check_argument_type(u, type="numeric",
                                           length=2,
                                           lower=0,
                                           upper=1),
                       check_argument_type(u, type="matrix",
                                           ncol=2,
                                           lower=0,
                                           upper=1))
                  ,2)
    check_arg_all(list(check_argument_type(copula, type="cyl_copula"),
                       check_argument_type(copula, type="Copula"))
                  ,2)
    check_arg_all(check_argument_type(cond_on,
                                      type="numeric",
                                      values = c(1,2))
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  u <- matrix(ncol=2,u)
  length <- nrow(u)
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]
  Ccond<-rep(1,length)


  if(cond_on==2){
    for (i in 1:length) {
      cond_func_v <- function(u){
        integrand <- function(x) {
          dcylcop(c(x[1], v[i,]), copula)
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
          dcylcop(c(u[i,], y[1]), copula)
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


#-----Calculate the conditional copula of a Copula-object.-------------------------------------------
#
#' @rdname ccylcop
#' @export
setMethod("ccylcop", "Copula", function(u, copula, cond_on, inverse) {
  if(dim(copula)!=2) stop("cylcop::ccylcop() works only for copulas of dimension 2. For other copula-objects, try copula::cCopula().")

  if(cond_on==2){
    u_swap <- matrix(ncol=2,c(u[,2],u[,1]))
    tryCatch({
      result <- copula::cCopula(u_swap,copula,indices =2, inverse = inverse) %>% as.numeric()
      if(any(is.na(result))) stop()
      result
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
      result <- copula::cCopula(u,copula,indices =2, inverse = inverse) %>% as.numeric()
      if(any(is.na(result))) stop()
      result
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
