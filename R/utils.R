#' Compass Bearing of a Line Between 2 Points
#'
#' The angle between a line between 2 points in Euclidean 2-D space and the line from
#' (0,0) to (0,1) is calculated. In other words, the compass bearing of a line
#' between 2 points where north is 0. Angles increase in clockwise direction.
#' @param point1 \link[base]{numeric} \link[base]{vector} holding the
#'  x and y coordinates of the first point.
#' @param point2 \link[base]{numeric} \link[base]{vector} holding the
#' x and y coordinates of the second point.
#' @param fullcirc \link[base]{logical} value indicating whether the output
#' should be an angle on \eqn{[0, 2\pi)} or \eqn{[-\pi, \pi)}.
#'
#' @return If \code{fullcirc = FALSE}, the function returns a \link[base]{numeric}
#'   value (angle) from the interval \eqn{[-\pi, \pi)}. \cr
#'   If \code{fullcirc = TRUE}, the function returns a numeric value \link[base]{numeric} from the interval
#'   \eqn{[0, 2\pi)}.
#' @export
#'
#' @examples
#' bearing(c(3,5), c(1,4))
#' bearing(c(3,5), c(1,4), fullcirc = TRUE)
bearing <- function(point1, point2, fullcirc = FALSE) {

  #validate input
  tryCatch({
    check_arg_all(check_argument_type(point1,
                                      type="numeric",
                                      length = 2)
                  ,1)
    check_arg_all(check_argument_type(point2,
                                      type="numeric",
                                      length = 2)
                  ,1)
    check_arg_all(check_argument_type(fullcirc,
                                      type="logical",
                                      length=1)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  if (all(point1 == point2)) {
    warning(
      "Distance between points is 0, cannot calculate bearing. Arbitratily assigned it an angle of 0"
    )
    return(0)
  }
  delx <- point2[1] - point1[1]
  dely <- point2[2] - point1[2]

  b <- sign(delx)
  b[b == 0] <- 1  #corrects for the fact that sign(0) == 0
  angle = b * (dely < 0) * pi + atan(delx / dely)

  if(fullcirc==TRUE){
    angle <- half2full_circ(angle)
  }
  return(angle)
}


#' Calculate the Next Position in a Trajectory from a Turn Angle and a Step Length
#'
#' The xy-coordinates of a position in 2-D space is calculated from the angle
#' between that position and the 2 previous ones in the trajectory and the
#' distance between that position and the previous one.
#'
#' @param angle \link[base]{numeric} value of the turn angle or a
#'  \code{\link[circular]{circular}} object, either in \eqn{[0, 2\pi)} or in \eqn{[-\pi, \pi)}
#' @param steplength \link[base]{numeric} value giving the distance between the
#' position and the previous one.
#' @param prevp1 \link[base]{numeric} \link[base]{vector} holding the x and y
#' coordinates of the previous position.
#' @param prevp2 \link[base]{numeric} \link[base]{vector} holding the x and y
#' coordinates of the position before the previous one.
#'
#' @return The function returns a \link[base]{numeric} \link[base]{vector}
#' holding the x and y coordinates of the position
#' @export
#'
#' @examples
#'  angstep2xy(1.5*pi, 2, prevp1 = c(1, 4), prevp2 = c(2, 7.5))
#'  angstep2xy(-0.5*pi, 2, c(1, 4), c(2, 7.5))
angstep2xy <- function(angle, steplength, prevp1, prevp2) {
  tryCatch({
    check_arg_all(list(check_argument_type(angle,
                                      type="numeric",
                                      length = 1),
                       check_argument_type(angle,
                                           type="circular"))
                  ,2)
    check_arg_all(check_argument_type(steplength,
                                      type="numeric",
                                      length = 1)
                  ,1)
    check_arg_all(check_argument_type(prevp1,
                                      type="numeric",
                                      length = 2)
                  ,1)
    check_arg_all(check_argument_type(prevp2,
                                      type="numeric",
                                      length = 2)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  #convert angle to \eqn{[0, 2\pi)} if necessary
  angle <-
    tryCatch(
      half2full_circ(angle),
      error = function(e)
        as.double(angle)
    )

  prevbear <- bearing(prevp2, prevp1, fullcirc = TRUE)
  newbear <- (prevbear + angle) %% (2 * pi)

  ##calculate changes in x and y coordinates from position-1 to position
  delx <- sin(as.double(newbear)) * steplength
  dely <- cos(as.double(newbear)) * steplength
  pos = c(prevp1[1] + delx, prevp1[2] + dely)
  return(pos)
}



#' Convert Angle from Half Circle to Full Circle
#'
#' Converts an angle from the half circle (i.e. in the interval \eqn{[-\pi, \pi)})
#' to an angle on the full circle  (i.e. in the interval \eqn{[0, 2\pi)}).
#'
#' @param angle \link[base]{numeric} value of an angle or a
#' \code{\link{circular}}-objekt in \eqn{[-\pi, \pi)}.
#'
#' @return The \link[base]{numeric} value of the angle in \eqn{[0, 2\pi)}.
#' @export
#'
#' @examples
#' half2full_circ(-1 * pi) / pi
#' half2full_circ(-0.5 * pi) / pi
#' half2full_circ(-0 * pi) / pi
#' half2full_circ(0.5 * pi) / pi
#'
half2full_circ <- function(angle) {
  tryCatch({
    check_arg_all(list(check_argument_type(angle,
                                           type="numeric"),
                       check_argument_type(angle,
                                           type="circular"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  (as.double(angle) + (2 * pi)) %% (2 * pi)
}

#' Convert Angle from Full Circle to Half Circle
#'
#' Converts an angle from the full circle (i.e. in the interval \eqn{[0, 2\pi)})
#' to an angle on the half circle  (i.e. in the interval \eqn{[-\pi, \pi)}).
#'
#' @param angle \link[base]{numeric} value of an angle or a
#' \code{\link{circular}}-objekt in \eqn{[0, 2\pi)}.
#'
#' @return The \link[base]{numeric} value of the angle in \eqn{[-\pi, \pi)}.
#' @export
#'
#' @examples
#' full2half_circ(0 * pi) / pi
#' full2half_circ(0.5 * pi) / pi
#' full2half_circ(1 * pi) / pi
#' full2half_circ(1.5 * pi) / pi
#' full2half_circ(2 * pi) / pi
#'
#'
full2half_circ <- function(angle) {
  tryCatch({
    check_arg_all(list(check_argument_type(angle,
                                           type="numeric"),
                       check_argument_type(angle,
                                           type="circular"))
                  ,2)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  angle <- as.double(angle %% (2 * pi)) %>%
    modify_if(~ .x > pi, ~ .x - (2 * pi))
  return(angle)
}


# Print arguments of a function call
#
#
#
# @param fun The function.
# @param params A list of characters holding the arguments passed to that function (named or unnamed).
# @param start An integer value, indicating with which function argument printing starts.
# @param param_print An integer value, indicating how many arguments are printed.
#
# @return A string consisting of parameter names and values with wich \code{fun} was called
# @export
#
# @examples
# print_param(dnorm, list(sd=1),start=1)
print_param <-
  function(fun,
           params = list(),
           start = 2,
           param_print = seq(start, length(formals(fun)))) {

    #Replace arguments of fun with elements of params, starting with the "start" argument.
    #Give back a string with the argument names and their values of argument numbers
    #specified in param_print
    args <- formals(fun)
    if (is.null(args)) {
      return("none")
    }
    if (length(param_print) > length(args)) {
      stop("requested to print more parameters than function has arguments",
           error_sound())
    }
    if (length(params) > length(args)) {
      stop("specified more parameters than the function has arguments",
           error_sound())
    }


    # If the parameters in params are not named, replace the default parameters consecutively

    if (length(params) != 0 && is.null(names(params))) {
      for (i in 1:length(params)) {
        if(!is.numeric(params[[i]]) && !is.character(params[[i]])){
          params[[i]]<-"object"
        }
        args[[i + start - 1]] <- params[[i]]
      }
    }


    #If the parameters in params are named, replace the default parameters by name

    else if (length(params) != 0 &&
             !"" %in% names(params) && !is.null(names(params))) {
      for (i in 1:length(params)) {
        if (!names(params)[i] %in% names(args)) {
          stop(error_sound(),
               names(params)[i],
               " is not a valid argument name")
        }
        if(!is.numeric(params[[i]]) && !is.character(params[[i]])){
          params[[i]]<-"object"
        }
        args[[names(params)[i]]] <- params[[i]]
      }
    }
    else if (length(params) != 0 && "" %in% names(params)) {
      stop(error_sound(),
           "either have all function arguments specified by name or none")
    }


    # Print

    temp <-
      paste(names(args)[param_print[1]], "=", deparse(args[[param_print[1]]]))
    if (length(param_print) > 1) {
      for (i in 2:length(param_print)) {
        temp <-
          paste(temp, ", ", names(args)[param_print[i]], " = ", deparse(args[[param_print[i]]]), sep = "")
      }
    }
    temp
  }

# Get univariate distributions
#
#
# @param marg_type A character string corresponding to the type of univariate distribution, e.g. "norm"
#
# @return Give A list of functions of a univariate distribution (random gneration, density, distribution, quantiles)
# @export
#
# @examples
# get_marg("norm")
get_marg <- function(marg_type) {
  if(marg_type=="exponential") marg_type <- "exp"
  if(marg_type=="gauss" || marg_type=="normal") marg_type <- "norm"

  tryCatch({
    list(
      "d" = (paste0("d", marg_type) %>% get()),
      "p" = (paste0("p", marg_type) %>% get()),
      "q" = (paste0("q", marg_type) %>% get()),
      "r" = (paste0("r", marg_type) %>% get())
    )
  }, error = function(e) {
    stop(
      error_sound(),
      "*",
      marg_type,
      " e.g. d",
      marg_type,
      " is not a distribution that could be found in any curently loaded package"
    )
  })
}


#' Set Package Options
#'
#' Currently the only option is to toggle verbosity on or off.
#' @param silent \link[base]{logical}, suppress all sounds and messages.
#'
#' @return No output, only side effects.
#' @examples cylcop_set_option(silent = FALSE)
#' @seealso \code{\link{cylcop_get_option}()}
#' @aliases verbose, silent
#' @export
#'
cylcop_set_option <- function(silent=FALSE){
  tryCatch({
    check_arg_all(check_argument_type(silent,
                                           type="logical",
                                           length = 1)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  assign("silent", silent, envir=cylcop.env)
}



#' Get Package Options
#'
#' Currently the only option (\code{"silent"}) is to toggle verbosity on or off.
#'
#' @param option \link[base]{character} string, the name of the option.
#'
#' @return
#' The \link[base]{numeric} value of option. If no argument is provided, a list
#' of all options is printed.
#' @examples cylcop_get_option("silent")
#' cylcop_get_option()
#' @seealso \code{\link{cylcop_set_option}()}
#' @export
#'
#'
cylcop_get_option <- function(option=NULL){
  if(is.null(option)){
    result <- ls(cylcop.env)
  }
  else{
 result<-get(option,envir=cylcop.env)
  }
  return(result)
}



#--------------Set sounds----------------
#
# Sound functions are commented out and sound-files not provided at the moment
#



# Error sound
#
# @include aaaglobal.R
# @return Play the error sound
# @export
#
error_sound <- function() {
  # if(cylcop.env$silent==F){
  #   decide <- runif(1)
  #   if(decide>0.5){
  #   sound::play(sound::loadSample(system.file("extdata", "error.wav", package = "cylcop")))
  #   }
  #   else{
  #     sound::play(sound::loadSample(system.file("extdata", "confusion.wav", package = "cylcop")))
  #
  #   }
  # }
  # else{}
}

# Warning sound
#
# @include aaaglobal.R
# @return Play the warning sound
# @export
#
warning_sound <- function() {
  # if(cylcop.env$silent==F){
  #   sound::play(sound::loadSample(system.file("extdata", "warning.wav", package = "cylcop")))
  # }
  # else{}
}

# Waiting sound
#
# @include aaaglobal.R
# @return Play the waiting sound
# @export
#
waiting_sound <- function() {
  # if(cylcop.env$silent==F){
  #   sound::play(sound::loadSample(system.file("extdata", "waiting.wav", package = "cylcop")))
  # }
  # else{}
}

# Done sound
#
# @include aaaglobal.R
# @return Play the done sound
# @export
#
done_sound <- function() {
  # if(cylcop.env$silent==F){
  #   sound::play(sound::loadSample(system.file("extdata", "done.wav", package = "cylcop")))
  # }
  # else{}
}

check_argument_type <- function(argument,
                                type = c("numeric",
                                         "logical",
                                         "character",
                                         "Copula",
                                         "cyl_copula",
                                         "list",
                                         "matrix",
                                         "NULL",
                                         "density",
                                         "density.circular",
                                         "circular",
                                         "data.frame"),
                                integer =F,
                                length = NULL,
                                ncol = NULL,
                                values = NULL,
                                lower = NULL,
                                upper = NULL) {
  cond_lst <- list(
    arg_name = deparse(substitute(argument)),
    type_cond = T,
    length_cond = T,
    ncol_cond = T,
    value_cond = T,
    lower_cond = T,
    upper_cond = T,
    NULL_cond = T,
    integer_cond = T,
    all_good=F
  )
  if(rlang::is_missing(argument)){
    stop(        paste0("Argument ",
      deparse(substitute(argument)),
      " is missing without any default value."
    ))
  }

  cond1 <- all(type %in% is(argument))
  possible_types <- c("numeric",
                      "logical",
                      "character",
                      "Copula",
                      "cyl_copula",
                      "list",
                      "matrix",
                      "NULL",
                      "density",
                      "density.circular",
                      "circular",
                      "data.frame")

  unwanted_types <- possible_types[-match(type,possible_types)]
  cond2 <- !any(unwanted_types %in% is(argument))


  if (!cond1 || !cond2) {
    cond_lst$type_cond <- paste(type,collapse=", ",sep="")
    return(cond_lst)
  }

  if (!is.null(length)) {
    if (!length == length(argument)) {
      cond_lst$length_cond <- length
      return(cond_lst)
    }
  }

  if (!is.null(ncol)) {
    if (!ncol == ncol(argument)) {
      cond_lst$ncol_cond <- ncol
      return(cond_lst)
    }
  }

  if(integer && any(type %in% c("numeric", "matrix"))){
    cond_lst$integer_cond <- is.logical(all.equal(c(argument), as.integer(argument)))
  }

  if (!is.null(values)) {
    if (!all(argument %in% values)) {
      cond_lst$value_cond <- values
      return(cond_lst)
    }
  }


  lower_cond <- T
  if (!is.null(lower)) {
    if (is.matrix(argument) && length(lower) == ncol(argument)) {
      lower_cond <- rep(T, length(lower))
      for (i in 1:length(lower)) {
        lower_cond[i] <- all(lower[i] <= argument[, i])
      }
    } else{
      lower_cond <- (lower <= argument)
    }

    if (!all(lower_cond)) {
      cond_lst$lower_cond <- lower
      return(cond_lst)
    }


  }

  upper_cond <- T
  if (!is.null(upper)) {
    if (is.matrix(argument) && length(upper) == ncol(argument)) {
      upper_cond <- rep(T, length(upper))
      for (i in 1:length(upper)) {
        upper_cond[i] <- all(upper[i] >= argument[, i])
      }
    } else{
      upper_cond <- (upper >= argument)
    }
    if (!all(upper_cond)) {
      cond_lst$upper_cond <- upper
      return(cond_lst)
    }
  }

  cond_lst$all_good <- T
  return(cond_lst)

}


check_arg_all <- function(arg_type, type_num) {
  if (type_num < 2) {
    if(arg_type$all_good) return()
    if (!isTRUE(arg_type$type_cond)) {
      stop(

        paste0(
          arg_type$arg_name,
          " must be of type ",
          arg_type$type_cond,
          "."
        )
      )
    }
    cond_lst <-  arg_type
  } else{
    type_ind <- F
    for (i in 1:type_num) {
      if (isTRUE(arg_type[[i]]$type_cond)) {
        type_ind <- i
      }
    }
    if (isFALSE(type_ind)) {
      type_message <- arg_type[[1]]$type_cond
      for (i in 2:type_num) {
        type_message <- paste(type_message, "\nor", arg_type[[i]]$type_cond)
      }
      stop(
           paste0(arg_type[[1]]$arg_name, " must be of type ", type_message, "."))
    }
    cond_lst <-  arg_type[[type_ind]]
  }


  if(cond_lst$all_good) return()
  if (!isTRUE(cond_lst$length_cond)) {
    stop(

      paste0(
        cond_lst$arg_name,
        " must be of length ",
        cond_lst$length_cond,
        "."
      )
    )
  }
  if (!isTRUE(cond_lst$ncol_cond)) {
    stop(

      paste0(
        cond_lst$arg_name,
        " must have ",
        cond_lst$ncol_cond,
        " columns."
      )
    )
  }
  if (!isTRUE(cond_lst$integer_cond)) {
    stop(
      paste0(
        cond_lst$arg_name,
        " must be an integer."
      )
    )
  }

  if (!isTRUE(cond_lst$value_cond)) {
    stop(
      paste0(
        c(
          cond_lst$arg_name,
          " can only take the following values: ",
          paste(cond_lst$value_cond,collapse=", ",sep=""),
          "."
        ),
        collapse = ""
      ))
  }

  if (!isTRUE(cond_lst$lower_cond)) {
    if (length(cond_lst$lower_cond) < 2) {
      stop(

        paste0(
          cond_lst$arg_name,
          " must be larger than ",
          cond_lst$lower_cond,
          "."
        )
      )
    } else{
      lower_message <- cond_lst$lower_cond[1]
      for (i in 2:length(cond_lst$lower_cond)) {
        lower_message <-  paste0(lower_message, ", ", cond_lst$lower_cond[i])
      }
      stop(

        paste0(
          cond_lst$arg_name,
          " must be larger than ",
          lower_message,
          "."
        )
      )
    }
  }

  if (!isTRUE(cond_lst$upper_cond)) {
    if (length(cond_lst$upper_cond) < 2) {
      stop(

        paste0(
          cond_lst$arg_name,
          " must be smaller than ",
          cond_lst$upper_cond,
          "."
        )
      )
    } else{
      upper_message <- cond_lst$upper_cond[1]
      for (i in 2:length(cond_lst$upper_cond)) {
        upper_message <- paste0(upper_message, ", ", cond_lst$upper_cond[i])
      }
      stop(

        paste0(
          cond_lst$arg_name,
          " must be smaler than ",
          upper_message,
          "."
        )
      )
    }
  }
}
