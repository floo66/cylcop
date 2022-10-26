#' Estimate Copula Parameters from Correlation Measures
#'
#' This function implements a simple search of the parameter space of a
#' '\code{\linkS4class{cyl_copula}}' object to find the
#' parameter values that lead to a correlation that is closest to the correlation
#' in the data (\code{theta} and \code{x}). In some special cases of
#' '\code{\linkS4class{cyl_rect_combine}}' copulas, the parameter can be
#' obtained analytically from Kendall's tau of the data.
#'
#' @param copula \R object of class '\code{\linkS4class{cyl_copula}}'.
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable).
#' @param acc \link[base]{numeric} value, the interval of the copula parameter
#' at which to evaluate the correlation.
#' @param n \link[base]{numeric} value, the number of sample points at each
#' optimization step.
#' @param method \link[base]{character} string describing what correlation metric
#'  to use. Either a rank-based circular-linear correlation coefficient (\code{"cor_cyl"}),
#' mutual information (\code{"mi_cyl"}), or Kendall's tau (\code{"tau"}).
#'
#' @param ... Additional parameters (see individual methods).
#'
#' @details The code assumes that the correlation captured by the copula increases
#' monotonously with the copula parameter values. It starts with a parameter value close
#' to the minimum for that copula and calculates the correlation for a sample of size \code{n}
#' from that copula. Next, the parameter is doubled and again the correlation for a sample
#' of size \code{n} calculated. After this exponential search pattern, a binary search
#' is implemented similarly between the bounds found with the exponential search. For this
#' binary search, the interval between those bounds is split into small intervals of length
#' \code{acc}. Thus, smaller values of \code{acc} lead to higher accuracy.
#'
#' If a '\code{\linkS4class{cyl_rect_combine}}' copula has rectangles spanning
#' the entire unit square and as background the independence copula, Kendall's tau can be used
#' to analytically calculate the parameter value leading to the correlation of the data.
#' No search is necessary in this case. This makes it the recommended method to use
#' for those '\code{\linkS4class{cyl_rect_combine}}' copulas.
#'
#' See also individual methods (below) for more detailed explanations.
#'
#' @return \link[base]{numeric} \link[base]{vector} containing the estimated
#' parameter value(s).
#'
#' @examples set.seed(123)
#'
#' sample <- rcylcop(100, cyl_rect_combine(copula::frankCopula(2)))
#' optCor(cyl_rect_combine(copula::frankCopula()),
#'   theta = sample[,1],
#'   x = sample[,2],
#'   method = "tau"
#' )
#'
#' optCor(cyl_rect_combine(copula::frankCopula()),
#'   theta = sample[,1],
#'   x = sample[,2],
#'   method = "mi_cyl",
#'   n = 100
#' )
#'
#' optCor(cyl_rect_combine(copula::claytonCopula()),
#'   theta = sample[,1],
#'   x = sample[,2],
#'   method = "tau"
#' )
#'
#' optCor(cyl_quadsec(), theta = sample[,1], x = sample[,2], method = "mi_cyl")
#' optCor(cyl_quadsec(), theta = sample[,1], x = sample[,2], method = "cor_cyl")
#' optCor(cyl_quadsec(),
#'   theta = sample[,1],
#'   x = sample[,2],
#'   method = "cor_cyl",
#'   n = 100,
#'   acc = 0.001
#' )
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{mi_cyl}()}, \code{\link{cor_cyl}()}, \code{\link{optML}()},
#' \code{\link{opt_auto}()}, \code{copula::\link[copula]{fitCopula}()}.
#'
#' @export
#'
setGeneric("optCor",
           function(copula,
                    theta,
                    x,
                    acc = NULL,
                    n = 10000,
                    method,
                    ...){
             missing_flag <- FALSE
             if(rlang::is_missing(method)){
               method <- "tau"
               missing_flag <- TRUE
                 }
             #validate input
             tryCatch({
               check_arg_all(check_argument_type(copula,
                                                 type="cyl_copula")
               ,1)
               check_arg_all(check_argument_type(theta,
                                                 type="numeric")
                             ,1)
               check_arg_all(check_argument_type(x,
                                                 type="numeric")
                             ,1)

               check_arg_all(list(check_argument_type(acc,
                                                      type="NULL"),
                                  check_argument_type(acc,
                                                      type="numeric",
                                                      length=1,
                                                      lower=0))
               ,2)
               check_arg_all(check_argument_type(n,
                                                 type="numeric",
                                                 length = 1,
                                                 integer=T,
                                                 lower=1)
                             ,1)
               check_arg_all(check_argument_type(method,
                                                 type="character",
                                                 values=c("cor_cyl", "mi_cyl", "tau"),
                                                 length=1)
                             ,1)
             },
             error = function(e) {
               error_sound()
               rlang::abort(conditionMessage(e))
             }
             )
             if(length(theta)!=length(x)){
               stop(
                 error_sound(),
                 "theta and x must have the same length."
               )
             }

             if(missing_flag) method <- "missing"

             standardGeneric("optCor")},
           signature = "copula")



# Roughly estimate cyl_vonmises parameter from some correlation metric
#
# It only makes sense to optimize kappa, mu does not influence correlation
#
#' @describeIn  optCor only parameter \code{"kappa"} can be optimized, since parameter
#' \code{"mu"} does not influence the correlation.
#
# @param copula A \code{cyl_vonmises} object
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
# @param acc A numeric value, the interval of kappa at which to evaluate the
#   correlation. DEFAULT: 0.5
# @param n The number of sample points at each optimization step (kappa value). DEFAULT: 10000
# @param method A string of characters, describing what correlation metric to use.
# Either a rank-based circular-linear coefficient ("cor_cyl") or
# mutual information ("mi_cyl").
#
# @return The estimated value of kappa.
#' @export
#'
setMethod("optCor", "cyl_vonmises", function(copula,
                                             theta,
                                             x,
                                             acc,
                                             n,
                                             method = "cor_cyl") {

  if(method=="missing") method <- "cor_cyl"

  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi_cyl"
  )
    stop(error_sound(),
         "only methods 'cor_cyl' and 'mi_cyl' are implemented for this copula")
  if (is.null(acc))
    acc <- 0.5
  min <- 0
  max <- 1000
  param_name <- "kappa"
  if(cylcop.env$silent==F){
    message(
    " find parameter kappa\nmethod: ",
    method,
    "\naccuracy: ",
    acc,
    "\nnumber of samples, n, in each step: ",
    n,
    "\n"
  )}
  kappa <-
    search_cor(copula, data, method, acc, n, min, max, param_name)
  if(cylcop.env$silent==F){
    message("kappa approx. ", kappa, "\n")}
  return(kappa)
})



#' Roughly estimate cyl_quadsec parameter from some correlation metric
#'
#' @describeIn  optCor the absolute value of the parameter is optimized, positive
#' and negative values give the same correlation.
#'
# @param copula A \code{cyl_quadsec} object
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
# @param acc A numeric value, the interval of a at which to evaluate the
#   correlation. DEFAULT: 0.01
# @param n The number of sample points at each optimization step (value of a). DEFAULT: 10000
# @param method A string of characters, describing what correlation metric to use.
# Either a rank-based circular-linear coefficient ("cor_cyl") or
# mutual information ("mi_cyl").
#
# @return The estimated value of a.
#' @export
#'
setMethod("optCor", "cyl_quadsec", function(copula,
                                            theta,
                                            x,
                                            acc,
                                            n,
                                            method = "cor_cyl") {
  if(method=="missing") method <- "cor_cyl"

  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi_cyl"
      )
    stop(error_sound(),
         "only methods 'cor_cyl' and 'mi_cyl' are impelemented for this copula")
  if (is.null(acc))
    acc <- 0.01
  min <- 0
  max <- 1 / (2 * pi)
  param_name <- "a"
  if (cylcop.env$silent==FALSE){
    message(
    "find parameter ",
    param_name,
    "\nmethod: ",
    method,
    "\naccuracy: ",
    acc,
    "\nnumber of samples, n, in each step: ",
    n,
    "\n"
  )}
  a <-
    search_cor(copula, data, method, acc, n, min, max, param_name)
  if(cylcop.env$silent==F){
    message("a approx. ", a, " or ", -a, "\n")}
  return(a)
})



# Roughly estimate cyl_cubsec parameters from some correlation metric
#
#
# @param copula A \code{cyl_cubsec} object
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
# @param acc A numeric value, the interval of a and/or b at which to evaluate
#   the correlation. DEFAULT: 0.01
# @param n The number of sample points at each optimization step (value of a
#   and/or b). DEFAULT: 10000
# @param method A string of characters, describing what correlation metric to
#   use. Either a rank-based circular-linear coefficient ("cor_cyl") or mutual
#   information ("mi_cyl").
#' @describeIn  optCor optimization of parameters, \code{"a"} and \code{"b"},
#' can be done separately or simultaneously.
#' @param parameter For '\code{\linkS4class{cyl_cubsec}}' copulas: A character
#' string specifying which parameter of the copula to optimize,
#'   \code{"a"}, \code{"b"}, or \code{"both"}
#'
# @return The estimated value of \code("a") or \code("b"), or a \link[base]{numeric} \link[base]{vector} holding the
#   estimated values of a and b.
#'
#'
#' @export
#'
setMethod("optCor", "cyl_cubsec", function(copula,
                                           theta,
                                           x,
                                           acc,
                                           n,
                                           method = "cor_cyl",
                                           parameter = "both") {

  if(method=="missing") method <- "cor_cyl"

  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi_cyl"
  )
    stop(error_sound(),
         "only methods 'cor_cyl' and 'mi_cyl' are impelemented for this copula")
  if (is.null(acc))
    acc <- 0.01
  min <- -1 / (2 * pi)
  max <- 1 / (2 * pi)

  if (parameter == "a" || parameter =="b") {
    param_name <- parameter
    if(cylcop.env$silent==F){message(
      "find parameter ",
      param_name,
      "\nmethod: ",
      method,
      "\naccuracy: ",
      acc,
      "\nnumber of samples, n, in each step: ",
      n,
      "\n"
    )}
    param <-
      search_cor(copula,
                 data,
                 method,
                 acc,
                 n,
                 min,
                 max,
                 param_name
      )
    if(cylcop.env$silent==F){message(param_name ," approx. ", param, "\n")}
    return(param)
  }

  else if (parameter == "both") {
    if(cylcop.env$silent==F){
    message(
      "find parameters a and b\nmethod: ",
      method,
      "\naccuracy: ",
      acc,
      "\nnumber of samples, n, in each step: ",
      n,
      "\n"
    )}

    if (method == "cor_cyl") {
      fun <- cor_cyl
    }
    else if (method == "mi_cyl" ) {
      fun <- mi_cyl
    }

    #calculate correlation values for all possible combinations of a and b, find the minimum
    cor <- fun(theta = theta, x = x)
    a <- seq(min, max, by = acc)
    b <- seq(min, max, by = acc)
    input <- expand.grid(a, b)
    diff <- Inf
    if(cylcop.env$silent==F){
      message("total steps: ", nrow(input))}
    for (i in 1:nrow(input)) {
      if (i %% 10 == 0)
        if(cylcop.env$silent==F){
          message("\ncurrent step: ", i)}
      copula <- setCopParam(copula, param_val=as.matrix(input[i,]), param_name=c("a","b"))
      sample <- rcylcop(n, copula)
      test <- fun(c(sample[, 1]), c(sample[, 2]))
      if (abs(cor - test) < diff) {
        diff <- abs(cor - test)
        final_a <- input[i, 1]
        final_b <- input[i, 2]
      }
    }
    if(cylcop.env$silent==F){
      message("\nthe optimal parameter set is: a= ",
        round(final_a, 2),
        ", b= ",
        round(final_b, 2),
        "\n")}
    return(c(final_a,final_b))
  }
  else
    stop(error_sound(),
         "argument 'parameter' must be from c(\"a\",\"b\",\"both\")")
})



# Roughly estimate cyl_rot_combine parameter from some correlation metric
#
#
# @param copula A \code{cyl_rot_combine} object
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
# @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#   correlation. DEFAULT: 0.5
# @param n The number of sample points at each optimization step (parameter value). DEFAULT: 10000
# @param method A string of characters, describing what correlation metric to use.
#   Only mutual information ("mi_cyl") makes sense
#
# @return The estimated value of the copula parameter.
#' @describeIn  optCor the circular-linear correlation coefficient will give a
#' value close to 0 for any parameter value. It therefore only makes sense to
#' use \code{method = "mi_cyl"} for the optimization.
#' @export
#'
setMethod("optCor", "cyl_rot_combine", function(copula, theta, x, acc, n, method ="mi_cyl") {

  if(method=="missing") method <- "mi_cyl"

  data <- data.frame(theta, x)
  if (method == "cor_cyl")
    stop(
      error_sound(),
      "correlation will always be approximately 0 for these types of copulas.\n",
      "Try optimization with method=mi_cyl"
    )
  if (method != "cor_cyl" &&
      method != "mi_cyl"
  )
    stop(error_sound(),
         "only method 'mi_cyl' is impelemented for this copula")

  #adhoc method of only testing positive parameter values.
  #in the future also implement negative ones.
  min <- max(0,copula@param.lowbnd)
  max <- copula@param.upbnd
  if (is.null(acc))
    acc <- 0.5
  param_name <- copula@param.names[1]
  if(cylcop.env$silent==F){
    message(
    "find parameter ",
    param_name,
    "\nmethod: ",
    method,
    "\naccuracy: ",
    acc,
    "\nnumber of samples, n, in each step: ",
    n,
    "\n"
  )}
  a <- search_cor(copula, data, "mi_cyl", acc, n, min, max, param_name)
  if(cylcop.env$silent==F){
    message(param_name, "approx. ", a, "\n")}
  return(a)
})



# Roughly estimate cyl_rect_combine parameters from some correlation metric
#
#
# @param copula A \code{cyl_rect_combine} object
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
# @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#   correlation. DEFAULT: 0.5
# @param n The number of sample points at each optimization step (parameter value). DEFAULT: 10000
# @param method A string of characters, describing what correlation metric to
#   use. Either Kendall's tau ("tau", recommended), a rank-based circular-linear coefficient
#   ("cor_cyl"), or mutual information ("mi_cyl").
#' @param background For '\code{\linkS4class{cyl_rect_combine}}' copulas :
#' A \link[base]{logical} value describing whether to optimize
#' the parameter of the background copula, (\code{background = TRUE}) or
#' the one of the copula in the rectangles (\code{background = FALSE}).
#'
#' @describeIn  optCor if the rectangles span the entire unit square and the background
#' is the independence copula, it is recommended to use \code{method = "tau"}, since this
#' calculates the copula parameter analytically. If there is a background copula,
#' other than the independence copula, its parameter can be optimized by setting
#' \code{background=TRUE}.
#
# @return The estimated value of the copula parameter.
#' @export
#'
setMethod("optCor", "cyl_rect_combine", function(copula,
                                               theta,
                                               x,
                                               acc,
                                               n,
                                               method = "tau",
                                               background = FALSE) {

  if(method=="missing") method <- "tau"

  data <- data.frame(theta, x)

  low_rect <- copula@parameters[match(c("low_rect1", "low_rect2"),copula@param.names)]
  up_rect <- copula@parameters[match(c("up_rect1", "up_rect2"),copula@param.names)]
  if ((isTRUE(all.equal(low_rect,c(0, 0.5))) &&
      isTRUE(all.equal(up_rect, c(0.5, 1))) )||
      any(is(copula@background.cop) == "indepCopula")) {
    if(background)
    stop(
      error_sound(),
      "background = TRUE only makes sense when there is a background parameter to optimize."
      )
  }

  if(method=="tau"){
    if(any(is(copula@sym.cop)=="cyl_copula")) stop(
      error_sound(),
      "method = tau is only implemented for linear-linear copulae in the rectangles"
    )

    a <- optTau(copula,theta,x)
  }
  else{
    ind <- 1
    if(any(is(copula@sym.cop)=="cyl_vonmises")) ind <- 2
  ind <- intersect(which(stringr::str_starts(copula@param.names,"bg_",negate = TRUE)),
              which(stringr::str_ends(copula@param.names,"rect.",negate = TRUE)))[ind]
  min <- max(0,copula@param.lowbnd[ind])
  if (is.null(acc))
    acc <- 0.5
  if (!background) {
    max <- copula@param.upbnd[ind]
    param_name <- copula@param.names[ind]
    if(cylcop.env$silent==F){
      message(
      "find parameter ",
      param_name,
      "\nmethod: ",
      method,
      "\naccuracy: ",
      acc,
      "\nnumber of samples, n, in each step: ",
      n,
      "\n"
    )}
  }
  else {
    max <- copula@param.upbnd[which(stringr::str_starts(copula@param.names,"bg_"))[1]]
    param_name <- copula@param.names[which(stringr::str_starts(copula@param.names,"bg_"))[1]]
    if(cylcop.env$silent==F){
      message(
      "find background parameter ",
      param_name,
      "\nmethod: ",
      method,
      "\naccuracy: ",
      acc,
      "\nnumber of samples, n, in each step: ",
      n,
      "\n"
    )}
  }

  a <-
    search_cor(copula, data, method, acc, n, min, max, param_name)
  if(cylcop.env$silent==FALSE){
    message(param_name, "approx. ", a,"\n")}
  }

  return(a)
})



# Finds the value of a specified parameter which gives correlation closest to empirical copula
#
# This function assumes that correlation and MI increase monotonously with the copula parameter
#
# @param copula A \code{cyl_copula} object.
# @param data A data.frame containing circular and linear measurements.
# @param method A string of characters, describing what correlation metric to
#   use. Either a rank-based circular-linear coefficient ("cor_cyl") or mutual
#   information ("mi_cyl").
# @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#   correlation.
# @param n The number of sample points at each optimization step (parameter value).
# @param min The minimum value of parameter to test.
# @param max The maximum value of parameter to test.
# @param param_name A character string, the name of the parameter to be varied.
#
# @return The value of parameter \code{param_name} that gives correlation closest to the empirical copula
# of the data in \code{data}.
# @export
#
search_cor <-
  function(copula,
           data,
           method = c("cor_cyl", "mi_cyl"),
           acc,
           n,
           min,
           max,
           param_name) {

    if (method == "cor_cyl") {
      fun <- cor_cyl
    }
    else if (method == "mi_cyl") {
      fun <- mi_cyl
    }
    else
      stop(error_sound(), "provide method, 'cor_cyl' or 'mi_cyl'")

    #get empirical correlation or MI of the data
    cor <- fun(theta = c(data[, 1]), x = c(data[, 2]))

    bound <- min + acc
    if(bound >= max){
      stop(error_sound(),
           paste0("The interval of the copula parameter at which the correlation is calculated\n",
                 "(acc = ",acc,") is larger than the search range of the parameter (",
                 min, " - ",max,". Reduce acc.")
           )
    }
    copula <- setCopParam(copula,  param_val=bound, param_name=param_name)
    sample <- rcylcop(n, copula)
    current <- fun(c(sample[, 1]), c(sample[, 2]))


#exponential search for upper and lower bounds

    if(cylcop.env$silent==F){
      cat("Starting exponential search: \n")}
    while ((current < cor) && (bound < max)) {
      if(cylcop.env$silent==F){
        cat(param_name, ">", bound, "\n")}
      bound <- bound * 2
      if(bound>=max) break
      copula <- setCopParam(copula, param_val=bound, param_name=param_name)
      sample <- rcylcop(n, copula)
      current <- fun(c(sample[, 1]), c(sample[, 2]))
    }


#binary search between bounds, the possible parameter values in the search "grid" vals will have a spacing of acc

    vals <- seq(from = bound / 2,
                to = min(bound + 1, max),
                by = acc)
    l <- 1
    r <- length(vals)
    m = floor((l + r) / 2)
    if(cylcop.env$silent==F){
      cat("Starting binary search :\n")}
    while (m != l & m != r) {
      if(cylcop.env$silent==F){
        cat(vals[l], " < ", param_name, " < ",  vals[r], "\n")}
      copula <- setCopParam(copula, param_val=vals[m], param_name=param_name)
      sample <- rcylcop(n, copula)
      current <- fun(c(sample[, 1]), c(sample[, 2]))
      if (current < cor)
        l <- m
      else if (current > cor)
        r <- m
      m = floor((l + r) / 2)
    }
    return(vals[m])
  }



# Estimate copula parameters based on Kendall's tau of the empirical copula
#
# @param copula A \code{Copula} object (from the copula package) or a
#   \code{cyl_rect_combine} object.
# @param theta A numeric vector of angles (measurements of a circular
#   variable).
# @param x A numeric vector of steplengths (measurements of a linear variable).
#
# @return A vector holding the parameter estimates.
# @export
#
setGeneric("optTau",
           function(copula, theta, x)
             standardGeneric("optTau"),
           signature = "copula")



# @describeIn optTau Works only for rectangles spanning the entire unit square, see text.
#
setMethod("optTau", "cyl_rect_combine", function(copula, theta, x) {
  low_rect <- copula@parameters[match(c("low_rect1", "low_rect2"),copula@param.names)]
  up_rect <- copula@parameters[match(c("up_rect1", "up_rect2"),copula@param.names)]
  if (!isTRUE(all.equal(low_rect,c(0, 0.5))) ||
      !isTRUE(all.equal(up_rect, c(0.5, 1))) ||
      !any(is(copula@background.cop) == "indepCopula")) {
    stop(
      error_sound(),
      "method = tau is only implemented for rectangles spanning the entire unit square and independent background copula"    )
  }

  data <- data.frame(theta, x) %>% na.omit()
  #calculate empirical copula
  emp_cop <- as.data.frame(pobs(data, ties.method = "average"))
  #flip  upper or lower part of empirical copula , recalculate ranks and optimize
  if (copula@flip_up)
    emp_cop$theta <- modify_if(emp_cop$theta, ~ .x > 0.5, ~ 1 - .x)
  else
    emp_cop$theta <- modify_if(emp_cop$theta, ~ .x < 0.5, ~ 1 - .x)
  emp_cop <- pobs(emp_cop, ties.method = "average")
  result <- fitCopula(copula@sym.cop, emp_cop, method = "itau")
  return(result@estimate)
})



# @describeIn optTau Only here for consistency, just calls function from copula package.
#
setMethod("optTau", "Copula", function(copula, theta, x) {
  data <- data.frame(theta, x) %>% na.omit()
  #calculate empirical copula
  emp_cop <- pobs(data, ties.method = "average")
  result <- fitCopula(copula, emp_cop, method = "itau")
  return(result@estimate)
})
