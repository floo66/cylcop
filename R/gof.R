#' Calculate the Wasserstein distance.
#'
#' The second Wasserstein distance, based on the Euclidean distance,
#'  between two copula PDFs on a grid, or between a copula pdf and pseudo-observations.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param copula2 (optional) \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param data \link[base]{numeric} \link[base]{matrix} with 2 columns, holding the
#' angles and step lengths
#' @param n_sample \link[base]{integer} number of samples to use when comparing two copula PDFs
#' @param n_grid number of gridcells
#'
#' @return
#' \link[base]{integer}, the second Wasserstein distance
#' @export
#'
#' @examples
wasserstein <- function(copula, copula2=NULL, data=NULL, n_sample, n_grid=2500){

  if(is.null(data)&&is.null(copula2)){
    stop("provide either data or a second copula to compare copula to")
  }

  if(!is.null(data)&&!is.null(copula2)){
    stop("provide either data or a second copula to compare copula to")
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
    out_cop<- transport::semidiscrete(a=dens_pgrid_cop1, b=sample_wpp, p=2)
    distance <- out_cop$wasserstein_dist
  }else{
    dens_cop2 <- matrix(dcylcop(as.matrix(grid),copula2),ncol=n_grid)
    dens_pgrid_cop2 <- transport::pgrid(dens_cop2,boundary = c(0,1,0,1))
    distance <- transport::wasserstein(a=dens_pgrid_cop1,b=dens_pgrid_cop2,p=2,prob=T,method="networkflow")
  }
  return(distance)
}
