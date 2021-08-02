#' @include cyl_cop_class.R
NULL


#' An S4 Class of Bivariate vonMises Copulas
#'
#' This class contains circular-linear copulas that are based on the approach by
#' \insertCite{Johnson1978;textual}{cylcop} with a von Mises periodic function.
#' They are periodic in the circular dimension, u, but not symmetric with
#' respect to \eqn{u=0.5} i.e. there is no symmetry between positive and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by \code{\link{cyl_vonmises}()}.
#'
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the parameter values.
#' @slot param.names \link[base]{character} \link[base]{vector} holding the
#' parameter names.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower
#'  bounds of the parameters.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper
#'  bounds of the parameters.
#' @slot flip \link[base]{logical} value indicating whether the copula should
#' be rotated 90 degrees to capture negative correlation.
#'
#' @section Extends:
#' Class '\code{cyl_vonmises}' extends class '\code{\linkS4class{cyl_copula}}'.
#'
#' @references \insertRef{Johnson1978}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
#'
setClass("cyl_vonmises", contains = "cyl_copula", slots = "flip")




#' Construction of '\code{cyl_vonmises}' Objects
#'
#' Constructs a circular-linear von Mises copula according to
#' \insertCite{Johnson1978;textual}{cylcop} of class
#'  '\code{\linkS4class{cyl_vonmises}}'.
#' @param mu \link[base]{numeric} value giving the mean of the vonMises
#' function used to construct the copula.
#' @param kappa \link[base]{numeric} value giving the concentration of the
#' vonMises function used to construct the copula.
#' @param flip \link[base]{logical} value indicating whether the copula
#' should be rotated 90 degrees to capture negative correlation.
#'
#' @export
#'
#' @examples
#' cop <- cyl_vonmises(mu=pi, kappa=10, flip = TRUE)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' cop <- cyl_vonmises(mu=0, kappa=8, flip = FALSE)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' @references \insertRef{Johnson1978}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
cyl_vonmises <- function(mu = 0,
                         kappa = 1,
                         flip = FALSE) {

  lowbnd = c(-Inf, 0)
  upbnd = c(Inf, Inf)

  new(
    "cyl_vonmises",
    name = "vonMises copula",
    parameters = c(mu, kappa),
    param.names = c("mu", "kappa"),
    param.lowbnd = lowbnd,
    param.upbnd = upbnd,
    flip = flip
  )
}



#' Generate random samples
#' @rdname Cylcop
# @describeIn cyl_vonmises-class Generate random samples.
#' @export
setMethod("rcylcop", signature("numeric", "cyl_vonmises"), function(n, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  u <- runif(n)
  w <- runif(n)
  #Calcualte the inverse of the conditional distribution of V given u, C_u(v) and
  #evaluate it at w

  v <- cylcop::ccylcop(matrix(ncol=2,c(u,w)),copula,cond_on=1,inverse = T)
  if (copula@flip)
    u <- 1 - u
  cop_uv <- cbind(u, v)
  return(cop_uv)
})



#' Calcualte density
#' @rdname Cylcop
# @describeIn cyl_vonmises-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_vonmises"), function(u, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  #drop=F, so it aslo works with single numbers
  v <- u[, 2, drop = F]
  if (copula@flip)
    u <- 1 - u[, 1, drop = F]
  else
    u <- u[, 1, drop = F]
  #Johnson and Wehrly's 1978 equation
  circular::dvonmises(x = circular::circular(2 * pi * (u - v)),
            mu = circular::circular(mu),
            kappa = kappa) %>%
    as.double() * 2 * pi
})



#' Calcualte distribution
#' @rdname Cylcop
# @describeIn cyl_vonmises-class Calculate the distribution.
#' @export
setMethod("pcylcop", signature("matrix", "cyl_vonmises"), function(u, copula) {
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  #Slow version with numerical integration
  # integrand <- Vectorize(function(x, v, mu, kappa) pvonmises(1.99999999*pi*v-x, mu = mu, kappa = kappa, from = -x)) #1.99999 instead of 2, because function is defined on
  # out <- integrate(f = integrand, lower = 0, upper = 2*pi*u, v = v, mu = mu, kappa = kappa) #works only if f accepts vecor input
  # (1/(2*pi))*out$value

  #Fast version with sum of Bessel functions
  ppvm <- function(x, kappa) {
    tol <- 1e-20
    flag <- TRUE
    j <- 1
    sum <- 0
    while (flag) {
      term <- (besselI(x = kappa, nu = j,expon.scaled = FALSE)
               * cos(j * x)) / (j * j)
      sum <- sum + term
      j <- j + 1
      if (all(abs(term) < tol))
        flag <- FALSE
    }
    return(sum / (pi * besselI(x = kappa, nu = 0, expon.scaled = FALSE)))
  }

  if (!copula@flip) {
    cdf <- 1 / (2 * pi) *
      (
        v * 2 * pi * u + ppvm(2 * pi * (u - v) - mu, kappa) -
          ppvm(2 * pi * u - mu, kappa) -
          ppvm(-2 * pi * v - mu, kappa) +
          ppvm(-mu, kappa)
      )
  }
  else{
    #if the copula is flipped, instead of integrating from (0,0) to (u,v) we need to integrate from  (1-u,0) to (1,v).
    #Use the C-volume of the unflipped copula for that, it recursively calls this pcylcop function.
    unflipped <- copula
    unflipped@flip <- FALSE
    cdf <- prob(unflipped, l = c((1 - u), 0), u = c(1, v))
  }
  return(c(cdf))
})



#' Condtional copula
#' @rdname ccylcop
# @describeIn cyl_vonmises-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_vonmises"), function(u, copula, cond_on=2, inverse=F) {
  u_orig <- matrix(ncol=2,u)
  mu <- copula@parameters[1]
  kappa <- copula@parameters[2]
  v <- u_orig[, 2, drop = F]
  u <- u_orig[, 1, drop = F]

  if(cond_on==2){
    if(!copula@flip){
      if(inverse==F){
        result <-
          map2_dbl(u, v, ~ circular::pvonmises(2*pi*.x,
                                               mu = circular::circular(2 * pi * .y + mu),
                                               kappa = kappa,
                                               from = circular::circular(0)))
      }
      if(inverse==T){

        #this is the inverse of pvonmises(v, mu=2*pi*u), which is the integrand in the copula cdf function,
        #(1/(2*pi))*qvonmises(pvonmises(2*pi*v, mu=2*pi*u),mu=2*pi*u)=v
        result <- (1 / (2 * pi)) *
          map2_dbl(u, v, ~ circular::qvonmises(.x,
                                               mu = circular::circular(2 * pi * .y + mu),
                                               kappa = kappa,
                                               from = circular::circular(0)))
      }
    }
    else if (copula@flip){
      #if flip, we need to integrate the cdf from 1-u to 1 instead of 0 to u
      if(inverse==F){
        result <-
          map2_dbl(u, v, ~ circular::pvonmises(2*pi,
                                               mu = circular::circular(2 * pi * .y + mu),
                                               kappa = kappa,
                                               from = circular::circular((1-.x)*2*pi)))
      }
      if(inverse==T){
        #TODO find analytical formula
        result <- numerical_inv_conditional_cop(c(u,v),copula,cond_on = 2)
      }
    }
  }
  else if(cond_on==1){
    if (copula@flip)
      #we still integrate from 0 to v, but at 1-u instead of u
      u <- 1 - u
    if(inverse==F){
      result <-
        map2_dbl(u, v, ~ circular::pvonmises(2*pi*.y,
                                             mu = circular::circular(2 * pi * .x - mu),
                                             kappa = kappa,
                                             from = circular::circular(0)))
    }
    if(inverse==T){
      result <- (1 / (2 * pi)) *
        map2_dbl(u, v, ~ circular::qvonmises(.y,
                                             mu = circular::circular(2 * pi * .x - mu),
                                             kappa = kappa,
                                             from = circular::circular(0)))
    }
  }
  else stop("cond_on must be either 1 or 2")
  return(result%>%as.numeric())
})



#-----Change attributes of existing cyl_vonmises object.-------------------------------------------
#
#' @rdname setCopParam
# @describeIn cyl_vonmises-class Change attributes of existing object.
#' @export
setMethod("setCopParam", "cyl_vonmises", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  return(copula)
})
