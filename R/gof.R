#' Calculate the Wasserstein Distance
#'
#' The pth Wasserstein distance is calculated based on the Euclidean distance between
#' two copula PDFs on a grid, or between a copula PDF and pseudo-observations.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param copula2 \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param theta (alternatively) \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable) or "circular" component of pseudo-observations.
#' @param x (alternatively) \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable) or "linear" component of pseudo-observations.
#' @param n_grid \link[base]{integer} number of grid cells at which the PDF of the copula(s) is calculated
#' Default is 2500
#' @param p \link[base]{integer} power (1 or 2) to which the Euclidean distance
#' between points is taken in order to compute transportation costs.
#'
#' @details Note that when comparing 2 copula PDFs (i.e. \code{theta = NULL} and \code{x = NULL}),
#' the calculated Wasserstein distance will depend on the number of grid cells
#' (\code{n_grid}) used to approximate the PDFs. The distance will converge to a certain
#' value with a higher number of grid cells, but the computational time will also increase.
#' The default of 2500 seems to be a good (empirically determined) compromise.
#' The same is true when calculating the Wasserstein distance between a copula
#' PDF and pseudo-observations. There, it is also important to only compare distances
#' that use the same number of observations.
#'
#' The code is based on the functions \code{transport::\link[transport]{wasserstein}()}
#' and \code{transport::\link[transport]{semidiscrete}()}.
#'
#' @return
#' \link[base]{numeric}, the pth Wasserstein distance
#'
#' @examples
#' set.seed(1234)
#' copula1 <- cyl_quadsec(0.1)
#' copula2 <- cyl_rect_combine(copula::frankCopula(2))
#' wasserstein(copula=copula1,copula2 = copula2,p=2,n_grid=200)
#' wasserstein(copula=copula1,copula2 = copula1,p=2,n_grid=200)
#' wasserstein(copula=copula1,copula2 = copula::frankCopula(2),p=2,n_grid=200)
#'
#'  sample <- rjoint(100,
#'   copula1,
#'   marginal_1 = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_2 = list(name = "weibull", coef = list(3,4))
#' )
#'
#' wasserstein(copula=copula1, theta=sample[,1], x=sample[,2], n_grid=200)
#'
#' @export
#'
wasserstein <- function(copula, copula2=NULL, theta=NULL, x=NULL, n_grid=2500, p=2){
  #validate input
  tryCatch({
    check_arg_all(list(check_argument_type(copula,
                                           type="Copula"),
                       check_argument_type(copula,
                                           type="cyl_copula")
    )
    ,2)
    check_arg_all(list(check_argument_type(copula2,
                                           type="Copula"),
                       check_argument_type(copula2,
                                           type="cyl_copula"),
                       check_argument_type(copula2,
                                          type="NULL")
    )
    ,3)
    check_arg_all(list(check_argument_type(theta,
                                           type="numeric"),
                       check_argument_type(theta,
                                           type="NULL"))
                  ,2)
    check_arg_all(list(check_argument_type(x,
                                           type="numeric"),
                       check_argument_type(x,
                                           type="NULL"))
                  ,2)
    check_arg_all(check_argument_type(n_grid,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      length=1,
                                      integer=T,
                                      values=c(1,2))
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  if((is.null(theta)||is.null(x))&&is.null(copula2)){
    stop("provide either data (theta and x) or a second copula to compare copula to")
  }

  if((!is.null(theta)||!is.null(x))&&!is.null(copula2)){
    stop("provide either data (theta and x) or a second copula to compare copula to")
  }

  if(is.null(copula2)){
    if(is.null(theta) || is.null(x)){
      stop("if copula2 is NULL, step lengths and turn angles must be provided")
    }
  }

  p <- as.integer(p)
  n_grid <- as.integer(sqrt(n_grid))
  grid <- expand.grid(seq(1/(2*n_grid), 1-1/(2*n_grid), 1/n_grid),
                      seq(1/(2*n_grid), 1-1/(2*n_grid), 1/n_grid))
  dens_cop1 <- matrix(dcylcop(as.matrix(grid),copula),ncol=n_grid)
  dens_cop1 <-  dens_cop1/sum(dens_cop1)
  dens_pgrid_cop1 <- transport::pgrid(dens_cop1,boundary = c(0,1,0,1))

  if(is.null(copula2)){
    data <- cbind(theta, x) %>% na.omit()
    n_sample <- nrow(data)
    pseudo_obs <- copula::pobs(data, ties.method = "average")
    sample_wpp <- transport::wpp(matrix(c(pseudo_obs[,2],1-pseudo_obs[,1]),ncol=2),
                                 mass=rep(dens_pgrid_cop1$totcontmass/n_sample,n_sample))
    out_cop<- transport::semidiscrete(a=dens_pgrid_cop1, b=sample_wpp, p=p)
    distance <- out_cop$wasserstein_dist
  }else{
    dens_cop2 <- matrix(dcylcop(as.matrix(grid),copula2),ncol=n_grid)
    dens_cop2 <- dens_cop2/sum(dens_cop2)
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





#' Cramér-von-Mises criterion
#'
#' Calculate the Cramér-von-Mises criterion
#' with a p-value (via parametric bootstrapping) to assess the goodness of fit
#' of a parametric copula compared to the empirical copula of the data.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}' or
#' '\code{\linkS4class{Copula}}' (package '\pkg{copula}'.
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable) or "circular" component of pseudo-observations.
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable) or "linear" component of pseudo-observations.
#' @param n_bootstrap  \link[base]{integer} number of bootstrap replicates. If
#' \code{n_bootstrap} is smaller than 1, no p-value is calculated.
#' @param parameters \link[base]{vector} of \link[base]{character} strings
#' holding the names of the parameters to be optimized when using the bootstrap
#' procedure.
#'   These can be any parameters in \code{copula@@parameters}. Default is to
#'   optimize the first 2 parameters. \code{parameters} has no effect if \code{copula}
#'   is of class '\code{\linkS4class{Copula}}' (package '\pkg{copula}'
#' @param optim.method \link[base]{character} string, optimizer used in
#' \code{\link[stats]{optim}()}, can be
#'  \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"},
#'  \code{"SANN"}, or \code{"Brent"}.
#' @param optim.control \link[base]{list} of additional controls passed to
#' \code{\link[stats]{optim}()}.
#'
#' @details The Cramér-von Misses criterion is calculated as the sum of the squared
#' differences between the empirical copula and the parameteric copula, \code{copula},
#' evaluated at the pseudo-observations obtained from \code{theta} and \code{x}.
#' If the bootstrap procedure is used, a random sample is drawn from \code{copula}
#' and converted to pseudo-observations. A new (set of) copula parameter(s) is then
#' fit to those pseudo-observations using maximum likelihood (function
#' \code{cylcop::\link[cylcop]{optML}()}).
#'
#' @references \insertRef{Genest2008}{cylcop}
#'
#' @return
#' A list of length 2 containing the Cramér-von Mises criterion and the p-value.
#' @export
#'
#' @examples
#'
#' sample <- rcylcop(100,cyl_cubsec(0.1, 0.1))
#'
#' opt_cop <- optML(copula = cyl_quadsec(),
#'   theta = sample[,1],
#'   x = sample[,2],
#'   parameters = "a",
#'   start = 0
#' )$copula
#' Cramer_vonMises(opt_cop,
#'   theta = sample[,1],
#'   x = sample[,2],
#'   n_bootstrap=100)
#'
Cramer_vonMises <-function(copula, theta, x, n_bootstrap=1000,
                           parameters = NULL,
                           optim.method = "L-BFGS-B",
                           optim.control = list(maxit = 100)){

  #validate input
  tryCatch({
    check_arg_all(list(check_argument_type(copula,
                                           type="Copula"),
                       check_argument_type(copula,
                                           type="cyl_copula")
    )
    ,2)
    check_arg_all(check_argument_type(theta,
                                           type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(n_bootstrap,
                                      type="numeric",
                                      integer=T,
                                      length=1)
                  ,1)
    check_arg_all(list(check_argument_type(parameters,
                                           type="character"),
                       check_argument_type(parameters,
                                           type="NULL"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  #Get the index in copula@parameters of the parameters to be optimized and check the names

  data <- cbind(theta, x) %>% na.omit()

  #calculate empirical


  points <- data
  pseudo_obs <- pobs(points, ties.method = "average")
  emp_cop <- C.n(pseudo_obs,X=pseudo_obs, smoothing = "none",ties.method = "average")
  CvM <- sum((emp_cop-pcylcop(pseudo_obs,copula))^2)
  p_val <- NULL

if(n_bootstrap>1){
  if(any(is(copula)%in%"cyl_copula")){
  if(is.null(parameters)){
    if(length(copula@param.names)==1){
      parameters <- copula@param.names
    }else{
      parameters <- copula@param.names[1:2]
    }
  }

  start_val <- copula@parameters[match(parameters,copula@param.names)]
  par_num <-
    tryCatch(
      param_num_checked(copula, param_val = start_val, param_name = parameters),
      error = function(e) {
        #error_sound()
        rlang::abort(conditionMessage(e))
      }
    )
  }else{
  start_val <- copula@parameters
}

  #start timer
  ptm <- proc.time()
  time <- 0

  CvM_vec <- rep(0, n_bootstrap)
  for (i in 1:n_bootstrap) {
    curr_opt_cop = 1
    while(is.numeric(curr_opt_cop)){
      points <- rcylcop(nrow(data),copula)
      if(any(is(copula)%in%"cyl_copula")){
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
      }else{
        curr_opt_cop <-
          tryCatch(
            suppressMessages(copula::fitCopula(copula=copula,
                                               data=points,
                                               start = start_val,
                                               optim.method = optim.method,
                                               optim.control = optim.control
            )),
            warning = function(w) {
              curr_opt_cop+1
            })
      }
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



    if(cylcop.env$silent==F){
      #get time for 10 steps, if n_bootstrap is at least that large, to set up the progress bar
      if (i == 10) {
        time <- (proc.time() - ptm)[3] %>% as.double()
        time <- (time / 10 * n_bootstrap)
        if (time > 20) {
          message(
            "Estimated time to generate bootstrap sample: ",
            floor(time / 60),
            " minutes, ",
            round(time - 60 * floor(time / 60)),
            " seconds"
          )
          pb <- utils::txtProgressBar(min = 1, max = n_bootstrap)
        }
      }
      if (time > 20){
        utils::setTxtProgressBar(pb, i)
      }
    }

  }
  if(cylcop.env$silent==F){
    #close progress bar
    if (time > 20){
      close(pb)
      done_sound()
    }
  }

  p_val <- (1/(n_bootstrap+1))*length(which(CvM_vec>CvM))
}
  out <- list(CvM = CvM, p_val=p_val)
  return(out)
}
