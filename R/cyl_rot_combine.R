#' @include cyl_cop_class.R
NULL


#' An S4 Class of Circular-Linear Copulas generated from Linear Combinations
#'  of Copulas
#'
#'This class contains bivariate circular-linear copulas, generated from
#' linear-linear bivariate '\code{\linkS4class{Copula}}' objects of the package
#' '\pkg{copula}', by taking the arithmetic mean of the original copula and
#' the 90 deg rotated copula. This results in copulas that are periodic in the
#' circular dimension, u, and symmetric with respect to \eqn{u=0.5}, i.e. positive
#' and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by    \code{\link{cyl_rot_combine}()}.
#'
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the parameter values.
#' @slot param.names \link[base]{character} \link[base]{vector} the parameter names.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower bounds of the
#'   parameters.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper bounds of the
#'   parameters.
#' @slot orig.cop linear-linear 2-dimensional '\code{\linkS4class{Copula}}'
#' object of the package '\pkg{copula}'.
#' @slot shift \link[base]{logical} value indicating whether the (u-periodic)
#' copula should be shifted by 0.5 in u direction.
#'
#' @section Extends:
#' Class '\code{cyl_rot_combine}' extends class '\code{\linkS4class{Copula}}'.
#'
#' @references \insertRef{Nelsen2006}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
#'
setClass("cyl_rot_combine",
         contains = "cyl_copula",
         slots = c("orig.cop", "shift"))



#' Construction of '\code{cyl_rot_combine}' Objects
#'
#' Constructs a circular-linear copula of class
#' '\code{\linkS4class{cyl_rot_combine}}' from linear combinations
#'  of copulas.
#'
#' @param copula linear-linear 2-dimensional '\code{\linkS4class{Copula}}'
#' object of the package '\pkg{copula}'.
#' @param shift \link[base]{logical} value indicating whether the (u-periodic)
#' copula should be shifted by 0.5 in u direction.
#'
#' @export
#'
#' @examples
#' cop <- cyl_rot_combine(copula = copula::frankCopula(param = 3), shift = TRUE)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' cop <- cyl_rot_combine(copula = copula::claytonCopula(param = 10), shift = FALSE)
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' @references \insertRef{Nelsen2006}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
cyl_rot_combine <- function(copula, shift = FALSE) {
  #If the underlying linear-linear copula is a rotated copula, the parameters are not in copula, but copula@copula
  if (any(is(copula) == "rotCopula")) {
    base_copula <- copula@copula
    name <-
      class(copula@copula)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula (rotated)")
  }
  else if (any(is(copula) == "copula")) {
    base_copula <- copula
    name <-
      class(copula)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula")
  }
  else
    stop(
      error_sound(),
      "provide a (rotated) 'copula'-object from the 'copula'-package as input"
    )

  if(cylcop.env$silent==F){
    message(name,
      " made periodic in u, by taking arithmetic mean with 90 degree rotated ",
      name,
      "\n")}

  new(
    "cyl_rot_combine",
    name = paste("cyl_rot_combine of",name),
    parameters = base_copula@parameters,
    param.names = base_copula@param.names,
    param.lowbnd = base_copula@param.lowbnd,
    param.upbnd = base_copula@param.upbnd,
    orig.cop = copula,
    shift = shift
  )
}


#' cyl_rot_combine show method
#' @rdname show-cyl_copula-method
#' @export
setMethod("show", "cyl_rot_combine", function(object) {
  cat(object@name, "\n")
  for (i in seq_along(object@parameters)) {
    cat(object@param.names[i], "=", object@parameters[i], "\n")
  }
  if (object@shift)
    cat("Copula is periodically shifted by 0.5u")
})



#' Generate random samples
#' @rdname Cylcop
# @describeIn cyl_rot_combine-class Generate random samples.
#' @export
setMethod("rcylcop", signature("numeric", "cyl_rot_combine"), function(n, copula) {
  # Generate the periodic copula. We do this here and not when instantiating the copula to make
  # it easier to change parameters during MLE optimization.

  # Take linear linear copula and rotate it 90 degrees, i.e "flip" in u-direction
  copula90 <- rotCopula(copula@orig.cop, flip = c(TRUE, FALSE))
  # now take the average (i.e. a convex sum) to get a copula periodic in u
  period_cop <-
    mixCopula(list(copula@orig.cop, copula90), w = c(0.5, 0.5))

  # period_cop is a Copula object from the 'copula'-package and we can use the corresponding methods
  sample <- rcylcop(n, period_cop)
  if (n == 1L)
    sample <- cbind(sample[1], sample[2])
  if (copula@shift)
    sample[, 1] <- (sample[, 1] + 0.5) %% 1
  colnames(sample) <- c("u", "v")
  return(sample)
})



#' Calcualte density
#' @rdname Cylcop
# @describeIn cyl_rot_combine-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_rot_combine"), function(u, copula) {
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  # Generate the periodic copula. We do this here and not when instantiating the copula to make
  # it easier to change parameters during MLE optimization.

  # Take linear linear copula and rotate it 90 degrees, i.e "flip" in u-direction
  copula90 <- rotCopula(copula@orig.cop, flip = c(TRUE, FALSE))
  # now take the average (i.e. a convex sum) to get a copula periodic in u
  period_cop <-
    mixCopula(list(copula@orig.cop, copula90), w = c(0.5, 0.5))

  if (copula@shift)
    u <- (u + 0.5) %% 1
  # period_cop is a Copula object from the 'copula'-package and we can use the corresponding methods
  dcylcop(cbind(u, v), period_cop)
})



#' Calcualte distribution
#' @rdname Cylcop
# @describeIn cyl_rot_combine-class Calculate the distribution.
#' @export
setMethod("pcylcop", signature("matrix", "cyl_rot_combine"), function(u, copula) {
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  # Generate the periodic copula. We do this here and not when instantiating the copula to make
  # it easier to change parameters during MLE optimization.

  # Take linear linear copula and rotate it 90 degrees, i.e "flip" in u-direction
  copula90 <- rotCopula(copula@orig.cop, flip = c(TRUE, FALSE))
  # now take the average (i.e. a convex sum) to get a copula periodic in u
  period_cop <-
    mixCopula(list(copula@orig.cop, copula90), w = c(0.5, 0.5))

  # For explanations of these calculations for the shifted copula see text and graphics
  if (copula@shift) {
    cdf <- map2_dbl(u, v, function(u, v) {
      if (u >= 0.5) {
        cdf <- pcylcop(c(1, v), period_cop) -
          pcylcop(c(0.5, v), period_cop) +
          pcylcop(c(0.5 - (1 - u), v), period_cop)
      }
      else if (u < 0.5) {
        cdf <- pcylcop(c(1 - (0.5 - u), v), period_cop) -
          pcylcop(c(0.5, v), period_cop)
      }
    })
  }
  else
    cdf <- pcylcop(cbind(u, v), period_cop)
  return(cdf)
})



#' Conditional copula
#' @rdname ccylcop
# @describeIn cyl_rot_combine-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_rot_combine"), function(u, copula, cond_on=2, inverse=F) {
  u_orig <- matrix(ncol=2,u)
  length <- nrow(u)
  v <- u_orig[, 2, drop = F]
  u <- u_orig[, 1, drop = F]

  # Take linear linear copula and rotate it 90 degrees, i.e "flip" in u-direction
  copula90 <- rotCopula(copula@orig.cop, flip = c(TRUE, FALSE))
  # now take the average (i.e. a convex sum) to get a copula periodic in u
  period_cop <-
    mixCopula(list(copula@orig.cop, copula90), w = c(0.5, 0.5))


  if (!copula@shift) {
    result <- cylcop::ccylcop(u_orig,period_cop, cond_on, inverse)
  }
  else{
    result <- matrix(ncol=2,rep(1,length))
    result[u>=0.5,] <- cylcop::ccylcop(matrix(ncol=2,c(rep(1,nrow(v[u>=0.5])),v[u>=0.5])),period_cop, cond_on, inverse)+
      cylcop::ccylcop(matrix(ncol=2,c(rep(0.5,nrow(v[u>=0.5])),v[u>=0.5])),period_cop, cond_on, inverse)+
      cylcop::ccylcop(matrix(ncol=2,c(rep(u-0.5,nrow(v[u>=0.5])),v[u>=0.5])),period_cop, cond_on, inverse)
    result[u<0.5,] <- cylcop::ccylcop(matrix(ncol=2,c(rep(0.5,nrow(v[u<0.5])),v[u<0.5])),period_cop, cond_on, inverse)+
      cylcop::ccylcop(matrix(ncol=2,c(rep(u+0.5,nrow(v[u<0.5])),v[u<0.5])),period_cop, cond_on, inverse)
  }
  return(result)
})



#-----Change attributes of existing cyl_rot_combine object.-------------------------------------------
#
#' @rdname setCopParam
# @describeIn cyl_rot_combine-class Change attributes of existing object.
#' @export
setMethod("setCopParam", "cyl_rot_combine", function(copula, param_val, param_name) {
  if(is.null(param_name)) param_name<-copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val
  if (any(is(copula@orig.cop) == "rotCopula")) {
    copula@orig.cop@copula@parameters[param_num] <- param_val
  }
  else{
    copula@orig.cop@parameters[param_num] <- param_val
  }
  return(copula)
})
