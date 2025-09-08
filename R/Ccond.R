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
                       check_argument_type(copula, type="Copula",
                                           dimension = 2))
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
#' \code{u[,-cond_on]} conditional on the values of \code{u[,cond_on]}.
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
                       check_argument_type(copula, type="Copula",
                                           dimension = 2))
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

  bound_ind <- check_boundary(cbind(u,v),copula)

  for (i in seq_len(length)) {
    if(cond_on==2){
      cond_func <- function(x){
        ccylcop(cbind(x,v[i,]),cond_on = 2,inverse = F,copula = copula)
      }
    }else if(cond_on==1){
      cond_func <- function(x){
        ccylcop(cbind(u[i,],x),cond_on = 1,inverse = F,copula = copula)
      }
    }
    cdf_cond_inv <- function(y){uniroot(f = function(x){cond_func(x)-y},lower=0,upper=1)$root}


    Ccond[i] <- tryCatch({
      if(cond_on==2){
        cdf_cond_inv(u[i,]) %>% suppressWarnings()
      }else if(cond_on==1){
        cdf_cond_inv(v[i,]) %>% suppressWarnings()
      }
    }, error = function(e) {

      if(i %in% bound_ind){
        NaN
      }else{
        stop(e)
      }
    })

  }

  if(any(is.na(Ccond))){
    warning("NaN produced, because calculation on the boundaries requested.")
  }
  return(Ccond)
}


#-----Calculate the conditional copula of a Copula-object.-------------------------------------------
#
#' @rdname ccylcop
#' @export
setMethod("ccylcop", "Copula", function(u, copula, cond_on, inverse) {
  if(dim(copula)!=2) stop("cylcop::ccylcop() works only for copulas of dimension 2. For other copula-objects, try copula::cCopula().")
  calc_cCopula <- function(u, copula, cond_on, inverse){
    if("mixCopula" %in% is(copula)){
      if(inverse){ stop("do numeric")
      }
      sum <- 0
      for (i in seq_along(copula@w)) {
        sum <- sum + copula@w[i]*calc_cCopula(u, copula@cops[[i]], cond_on, inverse=F)
      }
      return(sum)
    }

    cop_orig <- copula
    u_adapt <- u
    complement_res <- F


    if(cond_on==1){
      if("rotCopula" %in% is(copula)){
        # if(inverse){ stop("do numeric")
        # }
        cop_orig <- copula@copula
        if(all(copula@flip==c(T,T))){
          u_adapt <- 1-u
          complement_res <- T
        }else if(all(copula@flip==c(T,F))){
          u_adapt[,1] <- 1-u[,1]
        }else if(all(copula@flip==c(F,T))){
          u_adapt[,2] <- 1-u[,2]
          complement_res <- T
        }
      }
    }else if (cond_on==2){
      u_adapt[,1] <- u[,2]
      u_adapt[,2] <- u[,1]
      if("rotCopula" %in% is(copula)){
        # if(inverse){ stop("do numeric")
        # }
        cop_orig <- copula@copula
        if(all(copula@flip==c(T,T))){
          u_adapt[,1] <- 1-u[,2]
          u_adapt[,2] <- 1-u[,1]
          complement_res <- T
        }else if(all(copula@flip==c(T,F))){
          u_adapt[,1] <- u[,2]
          u_adapt[,2] <- 1-u[,1]
          complement_res <- T
        }else if(all(copula@flip==c(F,T))){
          u_adapt[,1] <- 1-u[,2]
          u_adapt[,2] <- u[,1]
        }
      }
    }
    if(!any(c("archmCopula", "ellipCopula", "indepCopula","fhCopula") %in% is(cop_orig))){
      stop("do numeric")
    }
    result <- copula::cCopula(u_adapt,cop_orig,indices =2, inverse = inverse) %>% as.numeric()
    if(any(is.na(result))){
      bound_ind <- check_boundary(u_adapt,cop_orig)
      na_ind <- which(is.na(result))
      if(length(bound_ind)>0 && all(na_ind%in%bound_ind)){
        warning("NaN produced, because calculation on the boundaries requested.")
      }else{
        stop("do numeric")
      }
    }
    if(complement_res){
      result <- 1-result
    }
    result
  }

tryCatch({
    calc_cCopula(u, copula, cond_on, inverse)
  }, error = function(e) {
    if(inverse==F){
      warning("Analytical expression not available, conditional copula was calculated numerically")
      numerical_conditional_cop(u,copula,cond_on)

    }
    else{
      warning("Analytical expression not available, inverse was calculated numerically")
      numerical_inv_conditional_cop(u,copula,cond_on)
    }
  })

})
