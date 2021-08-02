#' An S4 Class of Bivariate Copulas on the Cylinder
#'
#'The class '\code{cyl_copula}' follows somewhat the structure of the class
#' '\code{\linkS4class{Copula}}' of the package '\pkg{copula}'. It contains
#' circular-linear copulas.
#'
#' @section Extended by:
#' '\code{cyl_copula}' is extended by the following classes:
#' * '\code{\linkS4class{cyl_vonmises}}': von Mises copulas.
#' * '\code{\linkS4class{cyl_quadsec}}': Copulas with quadratic sections.
#' * '\code{\linkS4class{cyl_cubsec}}': Copulas with cubic sections.
#' * '\code{\linkS4class{cyl_rot_combine}}': Linear combinations of copulas and their
#' 180 degree rotations.
#' * '\code{\linkS4class{cyl_rect_combine}}': Rectangular patchwork copulas.
#'
#' @section Objects from the Class:
#' Objects are created by the functions \code{\link{cyl_vonmises}()},
#' \code{\link{cyl_quadsec}()}, \code{\link{cyl_cubsec}()}, \code{\link{cyl_rot_combine}()},
#' and \code{\link{cyl_rect_combine}()}.
#'
#' @md
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the parameter values.
#' @slot param.names \link[base]{character} \link[base]{vector} holding the parameter names.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower bounds of the parameters.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper bounds of the parameters.
#'
#' @examples
#' cop <- cyl_quadsec(0.1)
#' is(cop)
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
setClass(
  "cyl_copula",
  slots = c(
    name = "character",
    parameters = "numeric",
    param.names = "character",
    param.lowbnd = "numeric",
    param.upbnd = "numeric"
  ),
  prototype = list(
    name = NA_character_,
    parameters = NA_real_,
    param.names = NA_character_,
    param.lowbnd = NA_real_,
    param.upbnd = NA_real_
  ),
  validity =
    function(object) {
      param <- object@parameters
      upper <- object@param.upbnd
      lower <- object@param.lowbnd
      lp <- length(param)
      if (!(lp == length(upper)))
        return("Parameter and upper bound have non-equal length")
      if (!(lp == length(lower)))
        return("Parameter and lower bound have non-equal length")
      intervChk <-
        function(par)
          all(is.na(param) | (lower <= param & param <= upper))
      if (!intervChk(param))
        return("Parameter value(s) out of bound")
      TRUE
    }
)


#' Print Information of '\code{cyl_copula}' Objects
#'
#' Methods for function \code{\link[methods]{show}()} in package \pkg{cylcop}
#'
#' @param object \R object of class '\code{\linkS4class{cyl_copula}}'.
#'
#' @export
#'
setMethod("show", "cyl_copula", function(object) {
  cat(object@name, "\n")
  for (i in seq_along(object@parameters)) {
    cat(object@param.names[i],"=",object@parameters[i],"\n")
  }
})


#' Plot '\code{cyl_copula}' Objects
#'
#' Methods for \code{\link[base]{plot}()} to draw a scatter plot of a
#' random sample from bivariate distributions from package \pkg{cylcop}.
#' @param x \R object of class '\code{\linkS4class{cyl_copula}}'.
#' @param n sample size of the random sample drawn from \code{x}.
#' @param ... additional arguments passed to \code{\link[base]{plot}()}.
#'
#' @examples set.seed(123)
#'
#' plot(cyl_quadsec(0.1))
#' plot(cyl_vonmises(0,2), n=100)
#'
#' @export
#'
setMethod("plot", c("cyl_copula", "missing"), function(x, n=1000, ...){
    sample <- rcylcop(n, copula = x)
    plot(sample[,2], sample[,1], xlim=0:1, ylim=0:1, xlab="v", ylab="u", main=x@name, ...)
}
)


#Print Copula Parameters
#
#Used e.g. for printing a short summary of the copula in \code{make_traj()} and other
#functions.
#
# @param copula A \code{cyl_copula} or \code{Copula} object.
#
# @return NULL, prints the first 3 parameters of \code{copula} and their values
#  in the console.
#
#
setGeneric("printCop",
           function(copula)
             standardGeneric("printCop"),
           signature = "copula")


# Print  '\code{cyl_copula}' Parameters
# @describeIn printCop Print \code{cyl_copula} parameters
setMethod("printCop", "cyl_copula", function(copula) {
  parameter_summary = paste(copula@param.names[1], " = ", round(copula@parameters[1],4))
  if (length(copula@parameters) > 1) {
    for (i in 2:length(copula@parameters)) {
      parameter_summary <-
        paste(parameter_summary,
              ", ",
              copula@param.names[i],
              " = ",
              round(copula@parameters[i],4))
      if(i>2) break
    }
  }
  #if more than 3 parameters, print only first 3
  if (length(copula@parameters)>3){
    parameter_summary <-
      paste(parameter_summary, ", ...")
  }
  cat(
    sprintf(
      "%-7s %-20s parameters: %s",
      "Copula:",
      copula@name,
      parameter_summary,
      "\n"
    )
  )
})

# Print Copula parameterss
# @describeIn printCop Print \code{Copula} parameters
setMethod("printCop", "Copula", function(copula) {
  parameter_summary = paste(copula@param.names[1], " = ", round(copula@parameters[1],4))
  if (length(copula@parameters) > 1) {
    for (i in 2:length(copula@parameters)) {
      parameter_summary <-
        paste(parameter_summary,
              ", ",
              copula@param.names[i],
              " = ",
              round(copula@parameters[i],4))
      if(i>2) break
    }
  }
  #if more than 3 parameters, print only first 3
  if (length(copula@parameters)>3){
    parameter_summary <-
      paste(parameter_summary, ", ...")
  }
  #get the name of the copula
  if (any(methods::is(copula) == "rotCopula")) {
    name <- class(copula@copula)[1] %>%
      stringr::str_remove("Copula") %>%
      stringr::str_to_sentence(locale = "en") %>%
      paste("copula (rotated)")
  }
  else {
    name <-class(copula)[1] %>%
      stringr::str_remove("Copula") %>%
      stringr::str_to_sentence(locale = "en") %>%
      paste("copula")
  }

  #print description
  cat(
    sprintf(
      "%-7s %-20s parameters: %s",
      "Copula:",
      name,
      parameter_summary,
      "\n"
    )
  )
})



#' Calculate the C-Volume of a '\code{cyl_copula}' Copula
#'
#' This is a method corresponding to the generic \code{\link[copula]{prob}()} in the
#' '\pkg{copula}' package.
#'
#' @param x \R object of class '\code{\linkS4class{cyl_copula}}'.
#' @param l \link[base]{numeric} \link[base]{vector} of length 2 holding the coordinates of the
#' lower left corner in \eqn{[0,1]^2}.
#' @param u \link[base]{numeric} \link[base]{vector} of length 2 holding the coordinates of the
#' upper right corner in \eqn{[0,1]^2}.
#'
#' @return A \link[base]{numeric} in \eqn{[0,1]}, the probability that a draw from the
#' 2-dimensional copula \code{x} falls in the rectangle defined by \code{l} and
#' \code{u}.
#'
#' @aliases prob
#'
#' @examples cop <- cyl_quadsec(0.1)
#' prob(cop, l = c(0.1, 0.3), u = c(0.3, 0.9))
#'
#' @seealso \code{copula::\link[copula]{prob}}
#'
#' @export
setMethod("prob", "cyl_copula", function(x, l, u) {
  stopifnot(is.numeric(l), is.numeric(u),
            0 <= l, l <= u, u <= 1)
  pcylcop(c(l[1], l[2]), x) - pcylcop(c(l[1], u[2]), x) - pcylcop(c(u[1], l[2]), x) + pcylcop(c(u[1], u[2]), x)
})



#' Change Attributes of '\code{cyl_copula}' Objects
#'
#' These methods can be used, e.g. in other functions, to give users limited access
#' to the parameters of a copula.
#'
#' @details Note that for a rectangular patchwork copula
#' ('\code{\linkS4class{cyl_rect_combine}}')
#' the attribute \code{rectangles_symmetric} cannot be changed by \code{setCopParam()},
#' since rectangular patchwork copulas with symmetric rectangles are treated as
#' distinct from rectangular patchwork copulas with potentially asymmetric rectangles.
#' Therefore, when changing one of the bounds of the lower rectangle of such a copula,
#' the corresponding bound of the upper rectangle is automatically changed as well
#' (see examples).
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' @param param_val \link[base]{numeric} \link[base]{vector} holding the values to which the
#' parameters given in \code{copula@@parameters} should be changed.
#' @param param_name \link[base]{vector} of \link[base]{character} strings holding the
#'  names of the parameters to be changed.
#' @param ... additional arguments.
#'
#' @return A '\code{\linkS4class{cyl_copula}}' object with the changed parameters.
#' @name setCopParam
#'
#' @examples
#' cop <- cyl_rect_combine(copula::normalCopula(0.2),low_rect = c(0.1,0.4), up_rect="symmetric")
#' cop
#' cop <- setCopParam(cop, param_val = c(0.1, 0.3), param_name = c("rho.1", "low_rect2"))
#' cop <- cyl_rect_combine(copula::normalCopula(0.2),low_rect = c(0.1,0.4), up_rect=c(0.6,0.9))
#' cop
#' cop <- setCopParam(cop, param_val = 0.3, param_name = "low_rect2")
#' cop
#' @export
#'
setGeneric("setCopParam",
           function(copula, param_val, param_name=NULL, ...)
             standardGeneric("setCopParam"),
           signature = "copula")


#' Conditional Distributions of Circular-Linear Copulas
#'
#' Calculates the conditional distributions and their inverses of circular-linear
#' copulas and 2-dimensional linear-linear copulas.
#'
#' @param u \link[base]{matrix} (or \link[base]{vector}) of \link[base]{numeric}
#' values in \eqn{[0,1]^2}, containing as first column
#'  the circular (periodic) and as second the linear dimension.
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param cond_on column number of \code{u} on which the copula is conditioned. E.g if
#' \code{cond_on = 2}, the function calculates for each element in the first
#' column of \code{u} the copula conditional on the corresponding element in the
#' second column.
#' @param inverse \link[base]{logical} indicating whether the inverse of the conditional copula is
#' calculated.
#' @param ... additional arguments.
#'
#' @return A vector containing the values of the distribution of the copula at
#' \code{[u,-cond_on]} conditional on the values of \code{[u,cond_on]}.
#'
#' @details This is a generic that calls the function \code{copula::\link[copula]{cCopula}()}
#' for 2-dimensional '\code{\linkS4class{Copula}}' objects from the '\pkg{copula}'
#' package for which \code{copula::\link[copula]{cCopula}()} is available. If
#' \code{copula::\link[copula]{cCopula}()} is not available, the conditional copula is
#' calculated numerically. For '\code{\linkS4class{cyl_copula}}' objects,
#' the conditional copula is calculated analytically or numerically
#' (depending on the copula and  the values of \code{u}).
#' Note that the input arguments and the
#' output of \code{cylcop::ccylcop()} differ from those of
#' \code{copula::\link[copula]{cCopula}()}.
#'
#' @examples cop <- cyl_quadsec(0.1)
#' #calculate C_u(v) with u = 0.1 and v = 0.5
#' cylcop::ccylcop(u = c(0.1, 0.5), copula = cop, cond_on = 1, inverse = FALSE)
#' #calculate C^-1_v(u) with u = 0.1 and v = 0.5 and with u = 0.4 and v = 0.2
#' cylcop::ccylcop(u = rbind(c(0.1, 0.5), c(0.4, 0.2)), copula = cop, cond_on = 2, inverse = TRUE)
#'
#' @references \insertRef{Nelsen2006}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{copula::\link[copula]{cCopula}()}
#'
#' @name ccylcop
#'
setGeneric("ccylcop",
           function(u, copula, cond_on=2, inverse=F, ...){
             if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
             standardGeneric("ccylcop")
             },
           signature = "copula")



# Check new '\code{cyl_copula}' parameters
#
# @param copula \R object of class '\code{\linkS4class{cyl_copula}}'
# @param param_val A \link[base]{numeric} vector holding the values to which you want to
#   change certain parameters stored in \code{copula@@parameters}.
# @param param_name A vector vector of \link[base]{character} strings holding the names of the parameters to be
#   changed.
#
# @return The function returns a \link[base]{numeric} vector holding the indices of
#   \code{param_val} in \code{copula@@parameters}.
#
#
param_num_checked <- function(copula, param_val, param_name){
  #use abort instead of stop, because I want to discern error types when catching them in MLE
  if(length(param_name)!=length(param_val)){
    rlang::abort(message=paste0("number of parameters different from the number of specified parameter values"), "wrong_param_num")
  }

  param_num <- match(param_name,copula@param.names)

  if(any(is.na(param_num))){
    rlang::abort(message=paste0("the specified copula does not have a parameter ", param_name[which(is.na(param_num))]), "inexistent_param")
  }
  if(any(param_val < copula@param.lowbnd[param_num])){
    wrong<-which(param_val < copula@param.lowbnd[param_num])
    rlang::abort(message=paste0("parameter ", param_name[wrong], " is smaller than its lower bound"), "OOB_too_small")
  }
  if(any(param_val > copula@param.upbnd[param_num])){
    wrong<-which(param_val > copula@param.upbnd[param_num])
    rlang::abort(message=paste0("parameter ", param_name[wrong], " is larger than its upper bound"), "OOB_too_large")
  }
  return(param_num)
}


#' Distribution, Density, and Random Number Generation for Circular-Linear Copulas'
#'
#' Calculate the distribution (\code{pcylcop()}), the density (\code{dcylcop()}),
#' and generate random
#' samples (\code{rcylcop()}) of a '\code{\linkS4class{cyl_copula}}' object or a
#' '\code{\linkS4class{Copula}}' object (package '\pkg{copula}', only 2-dimensional).
#' For '\code{\linkS4class{Copula}}' objects \code{pcylcop()} and \code{rcylcop()}
#' just call the functions of the '\pkg{copula}' package
#' \code{\link[copula]{pCopula}()} and \code{\link[copula]{rCopula}()}, respectively.
#' The density is, however, calculated differently in \code{dcylcop()} and
#' \code{\link[copula]{dCopula}()}. The difference is
#'  that \code{copula::\link[copula]{dCopula}()}
#'  will return a density of 0 for points on the boundary of the unit square,
#'  whereas \code{dcylcop()} will return the correct density on the boundaries
#'  for both '\code{\linkS4class{cyl_copula}}' and '\code{\linkS4class{Copula}}' objects.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' or '\code{\linkS4class{Copula}}' (package '\pkg{copula}', only 2-dimensional).
#' @param u \link[base]{matrix} (or \link[base]{vector})  of \link[base]{numeric}
#'  values in \eqn{[0,1]^2}, containing as first column the circular (periodic) and
#'  as second the linear dimension
#' @param n number of random samples to be generated with \code{rcylcop()}.
#' @param log \link[base]{logical} indicating if the logarithm of the density
#'  should be returned.
#' @param ... additional arguments
#'
#' @returns The functions \code{pcylcop()} and \code{dcylcop()} give a \link[base]{vector} of
#' length \code{nrow(u)} containing the distribution and the density, respectively,
#'  at the corresponding values of \code{u}. Te function \code{rcylcop()} generates a
#'  \link[base]{matrix} with 2 columns and \code{n} rows containing
#' the random samples.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' rcylcop(5, cop)
#' pcylcop(c(0.3, 0.1), cop)
#' pcylcop(rbind(c(0.3, 0.1), c(0.2, 1)), cop)
#'
#' cop <- cyl_rot_combine(copula::frankCopula(2), shift = TRUE)
#' dcylcop(u = rbind(c(0.1, 0.4), c(1.0, 0.2)), copula = cop)
#' dcylcop(c(0.1, 0.3), cyl_quadsec(0.1), log = TRUE)
#'
#' cop <- copula::normalCopula(0.3)
#' copula::dCopula(c(.Machine$double.eps,0.2),cop)
#' copula::dCopula(c(0,0.2),cop)
#' dcylcop(c(.Machine$double.eps,0.2),cop)
#' dcylcop(c(0,0.2),cop)
#'
#' @references \insertRef{Nelsen2006}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{copula::\link[copula]{dCopula}()},
#' \code{copula::\link[copula]{pCopula}()},
#' \code{copula::\link[copula]{rCopula}()}.
#'
#' @name Cylcop
#' @aliases rcylcop pcylcop dcylcop
#
NULL

#' Calcualte distribution
#' @rdname Cylcop
#' @export
setGeneric("pcylcop", function(u, copula, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  ## here as well, 'outside' and 'on-boundary' are equivalent:
  u[] <- pmax(0, pmin(1, u))
  standardGeneric("pcylcop")
})

#' Random numbers
#' @rdname Cylcop
#' @export
setGeneric("rcylcop", function(n, copula, ...) standardGeneric("rcylcop"))

#' Copula Density
#'
#' @rdname Cylcop
#' @export
setGeneric("dcylcop", function(u, copula, log=FALSE, ...) {

  #Code is directly taken from copula::dCopula() of package copula.
  #The only difference is marked with !!!
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  stopifnot(dim(copula) == ncol(u))

  outside.01 <- function(u, strictly=TRUE) {
    if(strictly)
      apply(u, 1, function(x) any(x <	 0, 1 <	 x))
    else
      apply(u, 1, function(x) any(x <= 0, 1 <= x))
  }

  u.is.out <- outside.01(u, strictly=TRUE)## !!! copula package uses here strictly=FALSE
  if(any.out <- any(u.is.out, na.rm=TRUE))
    u[] <- pmax(0, pmin(1, u)) # <- "needed", as some methods give error
  r <- standardGeneric("dcylcop") # the result of calling  dcylcop-method(u, copula, ..)
  if(any.out) ##  outside cube  ==> zero mass, (but not when on cube boundary, in contrast to copula-package) :
    r[u.is.out & !is.na(u.is.out)] <- if(log) -Inf else 0.
  if(log==TRUE) r <-  log(r)
  r
},
package = "cylcop")

#' Calcualte density
#' @rdname Cylcop
#' @export
setMethod("dcylcop", signature("matrix", "Copula"), function(u, copula) {
  #workaround to get the correct density at the boundaries
  u[u < .Machine$double.eps] <- .Machine$double.eps
  u[u > 1-.Machine$double.eps] <- 1-.Machine$double.eps
  copula::dCopula(u,copula)
})

#' Random number generation
#' @rdname Cylcop
#' @export
setMethod("rcylcop", signature("numeric", "Copula"), function(n, copula) {
  copula::rCopula(n,copula)
})

#' Calcualte distribution
#' @rdname Cylcop
#' @export
setMethod("pcylcop", signature("matrix", "Copula"), function(u, copula) {
  copula::pCopula(u,copula)
})
