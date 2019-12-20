#' An S4 class of bivariate copulas on the cylinder
#'
#' The class \code{cyl_copula} follows somewhat the structure of the class
#' copula of the \code{copula}-package.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the parameters.
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



#' cyl_copula show method
#' @param object A \code{cyl_copula}-object
#' @describeIn cyl_copula-class What is printed
#' @export
setMethod("show", "cyl_copula", function(object) {
  cat(object@name, "\n")
  for (i in seq_along(object@parameters)) {
    cat(object@param.names[i],"=",object@parameters[i],"\n")
  }
})



#'Print copula parameters
#'
#'Used e.g. for printing a short summary of the copula in \code{make_traj()} and other
#'functions.
#'
#'@param copula A \code{cyl_copula} or \code{Copula} object.
#'
#'@return NULL, prints the first 3 parameters of \code{copula} and their values
#'  in the console.
#'@export
#'
setGeneric("printCop",
           function(copula)
             standardGeneric("printCop"),
           signature = "copula")


#' Print cyl_copula parameters
#' @describeIn printCop Print \code{cyl_copula} parameters
#' @export
setMethod("printCop", "cyl_copula", function(copula) {
  parameter_summary = paste(copula@param.names[1], " = ", copula@parameters[1])
  if (length(copula@parameters) > 1) {
    for (i in 2:length(copula@parameters)) {
      parameter_summary <-
        paste(parameter_summary,
              ", ",
              copula@param.names[i],
              " = ",
              copula@parameters[i])
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
      "%-18s %-20s parameters: %-60s",
      "Copula:",
      copula@name,
      parameter_summary,
      "\n"
    )
  )
})

#' Print Copula parameterss
#' @describeIn printCop Print \code{Copula} parameters
#' @export
setMethod("printCop", "Copula", function(copula) {
  parameter_summary = paste(copula@param.names[1], " = ", copula@parameters[1])
  if (length(copula@parameters) > 1) {
    for (i in 2:length(copula@parameters)) {
      parameter_summary <-
        paste(parameter_summary,
              ", ",
              copula@param.names[i],
              " = ",
              copula@parameters[i])
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
      "%-18s %-20s parameters: %-60s",
      "Copula:",
      name,
      parameter_summary,
      "\n"
    )
  )
})



#' Calculate the C-volume of a bivariate copula
#'
#' This is a method corresponding to the generic \code{prob} in the
#' \code{copula}-package
#'
#' @param x A \code{cyl_copula} object
#' @param l A numeric vector of length 2 holding the coordinates of the lower left corner in I^2
#' @param u A numeric vector of length 2 holding the coordinates of the upper right corner in I^2
#'
#' @export
setMethod("prob", "cyl_copula", function(x, l, u) {
  stopifnot(is.numeric(l), is.numeric(u),
            0 <= l, l <= u, u <= 1)
  pCopula(c(l[1], l[2]), x) - pCopula(c(l[1], u[2]), x) - pCopula(c(u[1], l[2]), x) + pCopula(c(u[1], u[2]), x)
})



#' Change attributes of \code{cyl_copula} objects
#'
#' @param copula A \code{cyl_copula} object.
#' @param param_val A numeric vector holding the values to which you want to
#'   change certain parameters stored in \code{copula@@parameters}.
#' @param param_name A vector vector of character strings holding the names of the parameters to be
#'   changed.
#' @param ... Additional arguments.
#'
#' @return A \code{cyl_copula}-object with the changed parameters.
#' @name setCopParam
#' @export
#'
setGeneric("setCopParam",
           function(copula, param_val, param_name=NULL,...)
             standardGeneric("setCopParam"),
           signature = "copula")



#' Check new \code{cyl_copula} parameters
#'
#' @param copula A \code{cyl_copula}-object
#' @param param_val A numeric vector holding the values to which you want to
#'   change certain parameters stored in \code{copula@@parameters}.
#' @param param_name A vector vector of character strings holding the names of the parameters to be
#'   changed.
#'
#' @return The function returns a numeric vector holding the indices of
#'   \code{param_val} in \code{copula@@parameters}.
#' @export
#'
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


#' Density, distribution, and random number generation for \code{cyl_copula} objects.
#'
#' These methods belong to the corresponding generics of the \code{copula} package.
#'
#' @param copula An object of class \code{cyl_copula} or \code{Copula}.
#' @param u Matrix of numeric values in I^2, containing as first collumn the circular (periodic) and as second the linear dimension
#' @param n Number of random samples to be generated with \code{rCopula}.
#' @param log Logical indicating whether the logarithm of the density should be returned. For \code{cyl_copula} always \code{log = FALSE}.
#' @param ... Further arguments
#' @usage dCopula(u, copula, log=FALSE, ...)
#' pCopula(u, copula, ...)
#' rCopula(n, copula, ...)
#' @name Copula
#' @aliases rCopula dCopula pCopula
#
NULL
