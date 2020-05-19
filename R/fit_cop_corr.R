
#' Roughly estimate cyl_copula parameters from some correlation metric
#'
#' See individual methods for more detailed explanations.
#'
#' @param copula A \code{cyl_copula} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of the copula parameter at which to evaluate the correlation.
#' @param n The number of sample points at each optimization step (parameter value).
#' @param method A string of characters, describing what correlation metric to use.
#' Either a rank-based circular-linear coefficient ("cor_cyl") or
#' mutual information ("mi").
#' @param ... Adittional parameters.
#'
#' @return A numeric value or vector containing the estimated parameter value(s).
#' @export
#'
setGeneric("optCor",
           function(copula,
                    theta,
                    x,
                    acc = NULL,
                    n = 10000,
                    method = c("cor_cyl", "mi"),
                    ...)
             standardGeneric("optCor"),
           signature = "copula")



#' Roughly estimate cyl_vonmises parameter from some correlation metric
#'
#' It only makes sense to optimize kappa, mu does not influence correlation
#'
#' @param copula A \code{cyl_vonmises} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of kappa at which to evaluate the
#'   correlation. DEFAULT: 0.5
#' @param n The number of sample points at each optimization step (kappa value). DEFAULT: 10000
#' @param method A string of characters, describing what correlation metric to use.
#' Either a rank-based circular-linear coefficient ("cor_cyl") or
#' mutual information ("mi").
#'
#' @return The estimated value of kappa.
#' @export
#'
setMethod("optCor", "cyl_vonmises", function(copula,
                                             theta,
                                             x,
                                             acc,
                                             n,
                                             method = c("cor_cyl", "mi")) {


  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi")
    stop(cylcop::error_sound(),
         "only methods 'cor_cyl' and 'mi' are impelemented")
  if (is.null(acc))
    acc <- 0.5
  min <- 0
  max <- Inf
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
#'
#' @param copula A \code{cyl_quadsec} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of a at which to evaluate the
#'   correlation. DEFAULT: 0.01
#' @param n The number of sample points at each optimization step (value of a). DEFAULT: 10000
#' @param method A string of characters, describing what correlation metric to use.
#' Either a rank-based circular-linear coefficient ("cor_cyl") or
#' mutual information ("mi").
#'
#' @return The estimated value of a.
#' @export
#'
setMethod("optCor", "cyl_quadsec", function(copula,
                                            theta,
                                            x,
                                            acc,
                                            n,
                                            method = c("cor_cyl", "mi")) {
  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi")
    stop(cylcop::error_sound(),
         "only methods 'cor_cyl' and 'mi' are impelemented")
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



#' Roughly estimate cyl_cubsec parameters from some correlation metric
#'
#'
#' @param copula A \code{cyl_cubsec} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of a and/or b at which to evaluate
#'   the correlation. DEFAULT: 0.01
#' @param n The number of sample points at each optimization step (value of a
#'   and/or b). DEFAULT: 10000
#' @param method A string of characters, describing what correlation metric to
#'   use. Either a rank-based circular-linear coefficient ("cor_cyl") or mutual
#'   information ("mi").
#' @param parameter A character string specifying which parameter to optimize,
#'   "a", "b", or "both"
#'
#' @return The estimated value of a, or b, or a numeric vector holding the
#'   estimated values of a and b.
#' @export
#'
setMethod("optCor", "cyl_cubsec", function(copula,
                                           theta,
                                           x,
                                           acc,
                                           n,
                                           method = c("cor_cyl", "mi"),
                                           parameter = "both") {
  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi")
    stop(cylcop::error_sound(),
         "only methods 'cor_cyl' and 'mi' are impelemented")
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
      fun <- function(theta, x)
        cylcop::cor_cyl(theta, x, plot = F)
    }
    else if (method == "mi") {
      fun <- cylcop::mi_binned
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
      sample <- rCopula(n, copula)
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
    stop(cylcop::error_sound(),
         "argument 'parameter' must be from c(\"a\",\"b\",\"both\")")
})



#' Roughly estimate cyl_rot_combine parameter from some correlation metric
#'
#'
#' @param copula A \code{cyl_rot_combine} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#'   correlation. DEFAULT: 0.5
#' @param n The number of sample points at each optimization step (parameter value). DEFAULT: 10000
#' @param method A string of characters, describing what correlation metric to use.
#'   Only mutual information ("mi") makes sense
#'
#' @return The estimated value of the copula parameter.
#' @export
#'
setMethod("optCor", "cyl_rot_combine", function(copula, theta, x, acc, n, method ="mi") {
  data <- data.frame(theta, x)
  if (method == "cor_cyl")
    stop(
      cylcop::error_sound(),
      "correlation will always be approximately 0 for these types of copulas.\n",
      "Try optimization with method=mi"
    )
  if (method != "cor_cyl" &&
      method != "mi")
    stop(cylcop::error_sound(),
         "only methods 'cor_cyl' and 'mi' are impelemented")

  #adhoc method of only testing positive parameter values.
  #in the future also implement negative ones.
  min <- 0
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
  a <- search_cor(copula, data, "mi", acc, n, min, max, param_name)
  if(cylcop.env$silent==F){
    message(param_name, "approx. ", a, "\n")}
  return(a)
})



#' Roughly estimate cyl_rect_combine parameters from some correlation metric
#'
#'
#' @param copula A \code{cyl_rect_combine} object
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#'   correlation. DEFAULT: 0.5
#' @param n The number of sample points at each optimization step (parameter value). DEFAULT: 10000
#' @param method A string of characters, describing what correlation metric to
#'   use. Either a rank-based circular-linear coefficient ("cor_cyl") or mutual
#'   information ("mi").
#' @param background A logical value describing whether to optimize the parameter of the background copula, (TRUE) or
#' the one of the copula in the rectangles(FALSE)
#'
#' @return The estimated value of the copula parameter.
#' @export
#'
setMethod("optCor", "cyl_rect_combine", function(copula,
                                               theta,
                                               x,
                                               acc,
                                               n,
                                               method = c("cor_cyl", "mi"),
                                               background = FALSE) {
  data <- data.frame(theta, x)
  if (method != "cor_cyl" &&
      method != "mi")
    stop(cylcop::error_sound(),
         "only methods 'cor_cyl' and 'mi' are impelemented")
  min <- 0
  if (is.null(acc))
    acc <- 0.5
  if (!background) {
    max <- copula@param.upbnd[1]
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
  if(cylcop.env$silent==F){
    message(param_name, "approx. ", a,"\n")}
  return(a)
})



#' Finds the value of a specified parameter which gives correlation closest to empirical copula
#'
#' This function assumes that correlation and MI increase monotonously with the copula parameter
#'
#' @param copula A \code{cyl_copula} object.
#' @param data A data.frame containing circular and linear measurements.
#' @param method A string of characters, describing what correlation metric to
#'   use. Either a rank-based circular-linear coefficient ("cor_cyl") or mutual
#'   information ("mi").
#' @param acc A numeric value, the interval of the copula parameter at which to evaluate the
#'   correlation.
#' @param n The number of sample points at each optimization step (parameter value).
#' @param min The minimum value of parameter to test.
#' @param max The maximum value of parameter to test.
#' @param param_name A character string, the name of the parameter to be varied.
#'
#' @return The value of parameter \code{param_name} that gives correlation closest to the empirical copula
#' of the data in \code{data}.
#' @export
#'
search_cor <-
  function(copula,
           data,
           method = c("cor_cyl", "mi"),
           acc,
           n,
           min,
           max,
           param_name) {

    if (method == "cor_cyl") {
      fun <- function(theta, x)
        cor_cyl(theta, x, plot = F)
    }
    else if (method == "mi") {
      fun <- mi_binned
    }
    else
      stop(cylcop::error_sound(), "provide method, 'cor_cyl' or 'mi'")

    #get empirical correlation or MI of the data
    cor <- fun(theta = c(data[, 1]), x = c(data[, 2]))

    bound <- min + acc
    copula <- setCopParam(copula,  param_val=bound, param_name=param_name)
    sample <- rCopula(n, copula)
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
      sample <- rCopula(n, copula)
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
      sample <- rCopula(n, copula)
      current <- fun(c(sample[, 1]), c(sample[, 2]))
      if (current < cor)
        l <- m
      else if (current > cor)
        r <- m
      m = floor((l + r) / 2)
    }
    return(vals[m])
  }



#' Estimate copula parameters based on Kendall's tau of the empirical copula
#'
#' @param copula A \code{Copula} object (from the copula package) or a
#'   \code{cyl_rect_combine} object.
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#'
#' @return A vector holding the parameter estimates.
#' @export
#'
setGeneric("optTau",
           function(copula, theta, x)
             standardGeneric("optTau"),
           signature = "copula")



#' @describeIn optTau Works only for rectangles spanning the entire unit square, see text.
#'
setMethod("optTau", "cyl_rect_combine", function(copula, theta, x) {
  low_rect <- copula@parameters[match(c("low_rect1", "low_rect2"),copula@param.names)]
  up_rect <- copula@parameters[match(c("up_rect1", "up_rect2"),copula@param.names)]
  if (!isTRUE(all.equal(low_rect,c(0, 0.5))) ||
      !isTRUE(all.equal(up_rect, c(0.5, 1))) ||
      !any(is(copula@background.cop) == "indepCopula")) {
    stop(
      cylcop::error_sound(),
      "Only implemented for rectangles spanning the entire unit square and independent background copula"
    )
  }

  data <- data.frame(theta, x) %>% na.omit()
  #calculate empirical copula
  emp_cop <- as.data.frame(pobs(data, ties.method = "average"))
  #flip  upper or lower part of empirical copula , recalcualte ranks and optimize
  if (copula@flip_up)
    emp_cop$theta <- modify_if(emp_cop$theta, ~ .x > 0.5, ~ 1 - .x)
  else
    emp_cop$theta <- modify_if(emp_cop$theta, ~ .x < 0.5, ~ 1 - .x)
  emp_cop <- pobs(emp_cop, ties.method = "average")
  fitCopula(copula@sym.cop, emp_cop, method = "itau")
})



#' @describeIn optTau Only here for consistency, just calls function from copula package.
#'
setMethod("optTau", "Copula", function(copula, theta, x) {
  data <- data.frame(theta, x) %>% na.omit()
  #calculate empirical copula
  emp_cop <- pobs(data, ties.method = "average")
  fitCopula(copula, emp_cop, method = "itau")
})
