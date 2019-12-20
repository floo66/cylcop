#' @include cyl_cop_class.R
NULL

# A class of copulas generated from linear-linear bivariate copulas of the 'copula'-package by
# taking the arithmetic mean of the original copula and the 90 deg rotated copula.
# This results in copulas that are periodic in u-direction (and not in v-direction) and therefore circular-linear
# They are symmetric between right and left turns, i.e. positive and negative angles.


#' An S4 class of circular-linear copulas from linear combinations of linear-linear copulas
#'
#' This class contains bivariate circular-linear copulas generated from
#' linear-linear bivariate copulas of the \code{copula}-package by taking the
#' arithmetic mean of the original copula and the 90 deg rotated copula. This
#' results in copulas that are periodic in the circular dimension, u, and
#' symmetric with respect to u=0.5, i.e. positive and negative angles.
#'
#' @section Objects from the Class:
#' Objects are created by \code{cyl_rot_combine()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
#' @slot orig.cop A linear-linear \code{Copula} object from the \code{copula}-package.
#' @slot shift A logical indicating whether the (u-periodic) copula shoud be shifted by 0.5 in u direction.
#'
#' @section Extends:
#' Class \code{cyl_rot_combine} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass("cyl_rot_combine",
         contains = "cyl_copula",
         slots = c("orig.cop", "shift"))



#' Construction of \code{cyl_quadsec} objects
#'
#' @param copula A linear-linear \code{Copula} object from the \code{copula}-package.
#' @param shift A logical indicating whether the (u-periodic) copula shoud be shifted by 0.5 in u direction.
#'
#' @export
#'
#' @examples
#' cyl_rot_combine(copula = copula::frankCopula(param = 3), shift = TRUE)
#' cyl_rot_combine(copula = copula::claytonCopula(param = 10), shift = FALSE)
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
      cylcop::error_sound(),
      "provide a (rotated) 'copula'-object from the 'copula'-package as input"
    )

  cat(name,
      "made periodic in u, by taking arithmetic mean with 90 degree rotated",
      name,
      "\n")

  new(
    "cyl_rot_combine",
    name = paste("'periodized'",name),
    parameters = base_copula@parameters,
    param.names = base_copula@param.names,
    param.lowbnd = base_copula@param.lowbnd,
    param.upbnd = base_copula@param.upbnd,
    orig.cop = copula,
    shift = shift
  )
}



#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_rot_combine"), function(n, copula) {
  # Generate the periodic copula. We do this here and not when instantiating the copula to make
  # it easier to change parameters during MLE optimization.

  # Take linear linear copula and rotate it 90 degrees, i.e "flip" in u-direction
  copula90 <- rotCopula(copula@orig.cop, flip = c(TRUE, FALSE))
  # now take the average (i.e. a convex sum) to get a copula periodic in u
  period_cop <-
    mixCopula(list(copula@orig.cop, copula90), w = c(0.5, 0.5))

  # period_cop is a Copula object from the 'copula'-package and we can use the corresponding methods
  sample <- rCopula(n, period_cop)
  if (n == 1L)
    sample <- cbind(sample[1], sample[2])
  if (copula@shift)
    sample[, 1] <- (sample[, 1] + 0.5) %% 1
  colnames(sample) <- c("u", "v")
  return(sample)
})



#' Calcualte density
#' @rdname Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_rot_combine"), function(u, copula) {
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
  dCopula(cbind(u, v), period_cop)
})



#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_rot_combine"), function(u, copula) {
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
        cdf <- pCopula(c(1, v), period_cop) -
          pCopula(c(0.5, v), period_cop) +
          pCopula(c(0.5 - (1 - u), v), period_cop)
      }
      else if (u < 0.5) {
        cdf <- pCopula(c(1 - (0.5 - u), v), period_cop) -
          pCopula(c(0.5, v), period_cop)
      }
    })
  }
  else
    cdf <- pCopula(cbind(u, v), period_cop)
  return(cdf)
})



#-----Change attributes of existing cyl_rot_combine object.-------------------------------------------
#
#' @rdname setCopParam
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
