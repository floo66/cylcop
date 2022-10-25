#' Density, Distribution, Quantiles and Random Number Generation for joint
#' distributions
#'
#' The bivariate joint distributions are described in terms of two marginal
#' distributions and a copula
#'
#' @param x \link[base]{matrix} (or \link[base]{vector}) of \link[base]{numeric}
#'  values giving the points (in 2 dimensions) where the density function is evaluated.
#' @param q \link[base]{matrix} (or \link[base]{vector}) of \link[base]{numeric}
#'  values giving the points (in 2 dimensions) where
#' the distribution function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rjoint()}.
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param marginal_1 named \link[base]{list} (for parametric estimates) or
#' a '\code{\link[stats]{density}}' object (for linear kernel density estimates)
#' or a '\code{\link[circular]{density.circular}}' object (for circular kernel density estimates).
#' The output of functions \code{\link{fit_angle}()} and \code{\link{fit_steplength}()}
#' can be used here directly.
#' @param marginal_2 This input is similar to \code{marginal_1}.
#'
#' @details If entered "by hand", the named lists describing the parametric distributions
#' (\code{marginal_1} and \code{marginal_2}) must contain 2 entries:
#' \enumerate{
#'    \item{name:
#' a \link[base]{character} string denoting the name of the distribution.
#' For a circular distribution, it can be \code{"vonmises"}, \code{"vonmisesmix"}, or
#'   \code{"wrappedcauchy"}. For a linear distribution, it must be a
#'   string denoting the name of a linear distribution in the environment, i.e. the name of its
#'    distribution function without the "p",
#'   e.g. "norm" for normal distribution}
#'    \item{coef: For a circular distribution coef is a (named) \link[base]{list} of
#' parameters of the circular
#' marginal distribution as taken by the functions
#' \code{\link[circular]{qvonmises}()}, \code{\link{qvonmisesmix}()},
#' or \code{\link{qwrappedcauchy}()}. For a linear distribution, coef is
#' a named list containing the parameters of the distribution given in \code{"name"}.}
#' }
#'
#' @return
#' \itemize{
#' \item{\code{djoint()}}{ gives a \link[base]{vector} of length \code{length(x)}
#'  containing the density at \code{x}.}
#' \item{\code{pjoint()}}{ gives a
#' \link[base]{vector} of length \code{length(q)} containing
#' the distribution function at the corresponding values of \code{q}.}
#' \item{\code{rjoint()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples.}
#'}
#'
#' @examples
#' cop <- copula::normalCopula(0.6)
#' marginal_1 <- list(name="exp",coef=list(rate=2))
#' marginal_2 <- list(name="lnorm", coef=list(0,0.1))
#'
#' sample <- rjoint(10,cop,marginal_1,marginal_2)
#' pjoint(sample,cop,marginal_1,marginal_2)
#' djoint(sample,cop,marginal_1,marginal_2)
#'
#' cop <- cyl_quadsec()
#' marginal_1 <- list(name="wrappedcauchy", coef=list(location=0,scale=0.3))
#' marginal_2 <- list(name="weibull",coef=list(shape=3))
#'
#' sample <- rjoint(10,cop,marginal_1,marginal_2)
#' marginal_1 <- fit_angle(theta=sample[,1], parametric=FALSE)
#' marginal_2 <- fit_steplength(x=sample[,2],parametric="lnorm")
#' pjoint(c(0.3*pi,4),cop,marginal_1,marginal_2)
#' djoint(c(0,2),cop,marginal_1,marginal_2)
#'
#' @name joint
#'
#'
#' @aliases rjoint pjoint djoint
#'
NULL




# Random numbers
#'
#' @rdname joint
#' @export
#'
rjoint <- function(n, copula, marginal_1, marginal_2){
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      integer=T,
                                      length=1,
                                      lower=1)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  marginals <- check_get_inp_joint(copula, marginal_1, marginal_2)


  cop_sample <- rcylcop(n, copula)
  if(marginals$name_1!="dens") {x <- do.call(marginals$funcs_1$q,
                                         c(list(p=cop_sample[,1]),
                                           marginals$param_1))}
  else{x <- do.call(marginals$funcs_1$q,
                           list(cop_sample[,1],
                                marginals$param_1))}
  if(marginals$name_2!="dens") {y <- do.call(marginals$funcs_2$q,
                                             c(list(p=cop_sample[,2]),
                                               marginals$param_2))}
  else{y <- do.call(marginals$funcs_2$q,
                           list(cop_sample[,2],
                                marginals$param_2))}
return(cbind(x=x,y=y))
}

# Density
#'
#' @rdname joint
#' @export
#'
djoint <- function(x, copula, marginal_1, marginal_2){ #validate input
  tryCatch({
    check_arg_all(list(check_argument_type(x, type="numeric",
                                           length=2),
                       check_argument_type(x, type="matrix",
                                           ncol=2))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  marginals <- check_get_inp_joint(copula, marginal_1, marginal_2)
  if(!is.matrix(x)) x <- rbind(x, deparse.level = 0L)
  y <- x[, 2, drop = T]
  x <- x[, 1,  drop = T]

  if(marginals$name_1=="vonmises"||marginals$name_1=="wrappedcauchy"){
    x <- circular::circular(x)
  }
  if(marginals$name_2=="vonmises"||marginals$name_2=="wrappedcauchy"){
    y <- circular::circular(y)
  }


  if (marginals$name_1 != "dens") {
    cdf_x <- do.call(marginals$funcs_1$p,
                     c(list(x),
                       marginals$param_1))
    pdf_x <- do.call(marginals$funcs_1$d,
                     c(list(x),
                       marginals$param_1))
  } else{
    cdf_x <- do.call(marginals$funcs_1$p,
                     list(x,
                          marginals$param_1))
    pdf_x <- do.call(marginals$funcs_1$d,
                     list(x,
                          marginals$param_1))
  }
  if (marginals$name_2 != "dens") {
    cdf_y <- do.call(marginals$funcs_2$p,
                     c(list(y),
                       marginals$param_2))
    pdf_y <- do.call(marginals$funcs_2$d,
                     c(list(y),
                       marginals$param_2))
  } else{
    cdf_y <- do.call(marginals$funcs_2$p,
                     list(y,
                          marginals$param_2))
    pdf_y <- do.call(marginals$funcs_2$d,
                     list(y,
                          marginals$param_2))
  }

  out <- dcylcop(cbind(cdf_x,cdf_y),copula = copula)*pdf_x*pdf_y

  return(out)
  }

# Distribution
#'
#' @rdname joint
#' @export
#'
pjoint <- function(q, copula, marginal_1, marginal_2){
  tryCatch({
    check_arg_all(list(check_argument_type(q, type="numeric",
                                           length=2),
                       check_argument_type(q, type="matrix",
                                           ncol=2))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  marginals <- check_get_inp_joint(copula, marginal_1, marginal_2)
  if(!is.matrix(q)) q <- rbind(q, deparse.level = 0L)
  x <- q[, 1,  drop = T]
  y <- q[, 2, drop = T]

  if(marginals$name_1=="vonmises"||marginals$name_1=="wrappedcauchy"){
    x <- circular::circular(x)
  }
  if(marginals$name_2=="vonmises"||marginals$name_2=="wrappedcauchy"){
    y <- circular::circular(y)
  }

  if (marginals$name_1 != "dens") {
    cdf_x <- do.call(marginals$funcs_1$p,
                     c(list(x),
                       marginals$param_1))
  } else{
    cdf_x <- do.call(marginals$funcs_1$p,
                     list(x,
                          marginals$param_1))
  }
  if (marginals$name_2 != "dens") {
    cdf_y <- do.call(marginals$funcs_2$p,
                     c(list(y),
                       marginals$param_2))
  } else{
    cdf_y <- do.call(marginals$funcs_2$p,
                     list(y,
                          marginals$param_2))
  }

  out <- pcylcop(cbind(cdf_x,cdf_y),copula = copula)

  return(out)
}



check_get_inp_joint <- function(copula, marginal_1, marginal_2){
tryCatch({

  check_arg_all(list(check_argument_type(copula,
                                         type="cyl_copula"),
                     check_argument_type(copula,
                                         type="Copula"))
                ,2)
  check_arg_all(list(check_argument_type(marginal_1,
                                         type="list"),
                     check_argument_type(marginal_1,
                                         type="density.circular"),
                     check_argument_type(marginal_1,
                                         type="density"))
                ,3)
  check_arg_all(list(check_argument_type(marginal_2,
                                         type="list"),
                     check_argument_type(marginal_1,
                                         type="density.circular"),
                     check_argument_type(marginal_1,
                                         type="density"))
                ,3)
},
error = function(e) {
  error_sound()
  rlang::abort(conditionMessage(e))
}
)
  check_get_margin <- function(marginal, marg_id){
  if(any(c("density.circular","density") %in% is(marginal))){
    parameter <- marginal
    name <- "dens"
  }else{
    if(!"name" %in% names(marginal)){
      stop(error_sound(),
           "marginal must be a density or density.circular object or a list
           containing the entries 'name' and 'coef'."
      )
    }
    if(!is.character(marginal$name)){
      stop(error_sound(),
           "In marginal: name must be of type character."
      )
    }
    if(!marginal$name=="unif"){
      if(!"coef" %in% names(marginal)){
        stop(error_sound(),
             "marginal must be a density or density.circular object or a list
           containing the entries 'name' and 'coef'."
        )
      }
    if(!is.list(marginal$coef)){
      stop(error_sound(),
           "In marginal: coef must be of type list."
      )
    }
    }

    if(marginal$name=="vonmises"){
      if(is.null(names(marginal$coef))){
        marginal$coef[[1]] <- circular::circular(marginal$coef[[1]])
      }else{
        if(!"mu"%in%names(marginal$coef)){
          stop(error_sound,"
               coef must contain the argument mu")
        }
        marginal$coef$mu <- circular::circular(marginal$coef$mu)
      }
    }


    parameter <- marginal$coef
    name <- marginal$name
  }

    marg_funcs <- get_marg(name)

     #Test if the precision of the approximation for the wrapped cauchy marginal distribution is sufficient
    if(name=="wrappedcauchy")
      tryCatch(do.call(marg_funcs$p, c(0, parameter)),
               warning = function(w) {
                 rlang::abort(conditionMessage(w))
               })

    if(cylcop.env$silent==F){
      #echo what functions and parameters are used
      if(name != "dens"){
        paste("\n", marg_id, "marginal:")
        message <- sprintf(
          "%-18s %-18s parameters: %s",
          paste("\n", marg_id, "marginal:"),
          name,
          print_param(marg_funcs$r, parameter)
        )
      }else{
        message <- paste("\n",marg_id," marginal: non-parametric density")
      }
    }



    return(list(marg_funcs,
                parameter,
                name,
                message))
  }


  if(cylcop.env$silent==F){
    printCop(copula)
  }

  out <- list(funcs_1=NA,param_1=NA,name_1=NA,funcs_2=NA,param_2=NA,name_2=NA)
  out <- list()
  temp <- check_get_margin(marginal_1, "first")
  out$funcs_1 <- temp[[1]]
  out$param_1 <- temp[[2]]
  out$name_1 <- temp[[3]]
  if(cylcop.env$silent==F){
    message(temp[[4]])
  }

  temp <- check_get_margin(marginal_2, "second")
  out$funcs_2 <- temp[[1]]
  out$param_2 <- temp[[2]]
  out$name_2 <- temp[[3]]
  if(cylcop.env$silent==F){
    message(temp[[4]])
  }
  return(out)
}
