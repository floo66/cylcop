#' Calculate the Wasserstein distance.
#'
#' The pth Wasserstein distance, based on the Euclidean distance,
#'  between two copula PDFs on a grid, or between a copula pdf and pseudo-observations.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param copula2 (optional) \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param data \link[base]{numeric} \link[base]{matrix} with 2 columns, holding the
#' angles and step lengths.
#' @param n_sample \link[base]{integer} number of samples to use when comparing two copula PDFs
#' @param n_grid \link[base]{integer} number of gridcells.
#' @param p \link[base]{integer} power (1 or 2) to which the Euclidean distance
#' between points is taken in order to compute transportation costs.
#'
#' @return
#' \link[base]{integer}, the pth Wasserstein distance
#' @export
#'
#' @examples
wasserstein <- function(copula, copula2=NULL, data=NULL, n_sample, n_grid=2500, p=2){

  if(is.null(data)&&is.null(copula2)){
    stop("provide either data or a second copula to compare copula to")
  }

  if(!is.null(data)&&!is.null(copula2)){
    stop("provide either data or a second copula to compare copula to")
  }

  if(!is.integer(p)){
    stop("p must be an integer")
  }
  if(!(p %in%c(1,2))){
    stop("p must be either 1 or 2")
  }


  n_grid <- as.integer(sqrt(n_grid))
  grid <- expand.grid(seq(1/(2*n_grid), 1-1/(2*n_grid), 1/n_grid),
                      seq(1/(2*n_grid), 1-1/(2*n_grid), 1/n_grid))
  dens_cop1 <- matrix(dcylcop(as.matrix(grid),copula),ncol=n_grid)
  dens_pgrid_cop1 <- transport::pgrid(dens_cop1,boundary = c(0,1,0,1))

  if(!is.null(data)){
    n_sample <- nrow(data)
    pseudo_obs <- copula::pobs(data, ties.method = "average")
    sample_wpp <- transport::wpp(matrix(c(pseudo_obs[,2],1-pseudo_obs[,1]),ncol=2),
                                 mass=rep(dens_pgrid_cop1$totcontmass/n_sample,n_sample))
    out_cop<- transport::semidiscrete(a=dens_pgrid_cop1, b=sample_wpp, p=p)
    distance <- out_cop$wasserstein_dist
  }else{
    dens_cop2 <- matrix(dcylcop(as.matrix(grid),copula2),ncol=n_grid)
    dens_pgrid_cop2 <- transport::pgrid(dens_cop2,boundary = c(0,1,0,1))
    if(p==2){
      distance <- transport::wasserstein(a=dens_pgrid_cop1,b=dens_pgrid_cop2,p=p,prob=T,
                                       method="networkflow")
    }else{
      distance <- transport::wasserstein(a=dens_pgrid_cop1,b=dens_pgrid_cop2,p=p,prob=T)
    }
  }
  return(distance)
}





#' Cramer-von-Mises criterion
#'
#' Calculate the Cramer-von-Mises criterion of a copula CDF
#' with a p-value (via parametric bootstrapping).
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable).
#' @param parameters \link[base]{vector} of \link[base]{character} strings
#' holding the names of the parameters to be optimized.
#'   These can be any parameters in \code{copula@@parameters}.
#' @param n_bootstrap  \link[base]{integer} number of bootstrap replicates.
#' @param optim.method \link[base]{character} string, optimizer used in
#' \code{\link[stats]{optim}()}, can be
#'  \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"},
#'  \code{"SANN"}, or \code{"Brent"}.
#' @param optim.control \link[base]{list} of additional controls passed to
#' \code{\link[stats]{optim}()}.
#'
#' @return
#' @export
#'
#' @examples
Cramer_vonMises <-function(copula, theta, x, parameters, n_bootstrap=1000,
                           optim.method="L-BFGS-B",
                           optim.control=list(maxit = 100)){
  #Get the index in copula@parameters of the parameters to be optimized and check the names

  data <- cbind(theta, x) %>% na.omit()

  start_val <- copula@parameters[which(copula@param.names %in% parameters)]
  par_num <-
    tryCatch(
      cylcop:::param_num_checked(copula, param_val = start_val, param_name = parameters),
      error = function(e) {
        #error_sound()
        rlang::abort(conditionMessage(e))
      }
    )

  #calculate empirical


  points <- data
  pseudo_obs <- pobs(points, ties.method = "average")
  emp_cop <- C.n(pseudo_obs,X=pseudo_obs, smoothing = "none",ties.method = "average")
  CvM <- sum((emp_cop-pcylcop(pseudo_obs,copula))^2)

  CvM_vec <- rep(0, n_bootstrap)
  for (i in 1:n_bootstrap) {
    cat(i, "of ",n_bootstrap,"\n")
    curr_opt_cop=1
    while(is.numeric(curr_opt_cop)){
      points <- rcylcop(nrow(data),copula)
      curr_opt_cop <-
        tryCatch(
          suppressMessages(cylcop::optML(copula=copula,theta = points[,1],
                                         x = points[,2],parameters=parameters,
                                         start = start_val,optim.method = optim.method,
                                         optim.control = optim.control
          )),
          warning = function(w) {
            curr_opt_cop+1
          })
      if(is.numeric(curr_opt_cop) &&curr_opt_cop >10)
        stop("failed to converge 10 times in a row")
    }
    curr_opt_cop <- curr_opt_cop$copula
    pseudo_obs <- pobs(points, ties.method = "average")
    emp_cop <- C.n(points,
                   X=pseudo_obs,
                   smoothing = "none",
                   ties.method = "average")
    CvM_vec[i] <- sum((emp_cop-pcylcop(points,curr_opt_cop))^2)
  }
  p_val <- (1/(n_bootstrap+1))*length(which(CvM_vec>CvM))
  return(p_val)
}
