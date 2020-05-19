

#' Estimate cyl_copula parameters according to maximum likelihood
#'
#' The code of this function is based on \code{copula::fitCopula()}.
#'
#' @param copula A \code{cyl_copula} object.
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param parameters A vector of character strings holding the names of the parameters to be optimized.
#'   These can be any parameters in \code{copula@@parameters}.
#' @param start A vector of staring values of the parameters.
#' @param lower OPTIONAL: A vector of lower bounds of the parameters.
#' @param upper OPTIONAL: A vector of upper bounds of the parameters.
#' @param optim.method optimizer used in \code{optim()}, can be
#'   "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", or "Brent"
#' @param optim.control A list of additional controls passed to \code{optim()}.
#' @param estimate.variance A logical valu, denoting whether to include an
#'   estimate of the variance (NOT YET IMPLEMENTED).
#' @param traceOpt A logical value, denoting whether to print information regarding
#'   convergence, current values, etc. during the optimization process.
#'
#' @return A list containing the same type of \code{cyl_copula} object as \code{copula},
#'  but with optimized parameters, the log likelihood and the AIC.
#' @export
#'
optML <-  function(copula,
                   theta,
                   x,
                   parameters,
                   start,
                   lower = NULL,
                   upper = NULL,
                   optim.method = "L-BFGS-B",
                   optim.control = list(maxit = 100),
                   estimate.variance = F,
                   traceOpt = F) {
# Check input and prepare data
  checked <-
    prep_n_check_ML(copula, theta, x, parameters, start, lower, upper)
  emp_cop <- checked[[1]]
  lower <- checked[[2]]
  upper <- checked[[3]]

# Find ML parameter estimates
  opt_cop <-
    fit_LL(
      copula,
      emp_cop,
      parameters,
      start,
      lower,
      upper,
      optim.method,
      optim.control,
      estimate.variance,
      traceOpt
    )
  return(opt_cop)
}



#' Check input and prepare data for MLE
#'
#' @param copula A \code{cyl_copula} object.
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param param_name A vector of character strings holding the names of the parameters to be optimized.
#'   These can be any parameters in \code{copula@@parameters}.
#' @param start A vector of staring values of the parameters.
#' @param lower OPTIONAL: A vector of lower bounds of the parameters.
#' @param upper OPTIONAL: A vector of upper bounds of the parameters.
#'
#' @return The function returns a list containing the empirical copula of the data and
#' reasonable upper and lower bounds (if not specified in the input). It also checks the input.
#' @export
#'
prep_n_check_ML <-
  function(copula,
           theta,
           x,
           param_name,
           start,
           lower,
           upper) {

    data <- cbind(theta, x) %>% na.omit()

    #calculate empirical copula
    emp_cop <- pobs(data, ties.method = "average")


# Get the index in copula@parameters of the parameters to be optimized and check the names

    par_num <-
      tryCatch(
        param_num_checked(copula, param_val = start, param_name = param_name),
        error = function(e) {
          cylcop::error_sound()
          cat("among the 'start' values: ")
          rlang::abort(conditionMessage(e))
        }
      )


# Check the lower bounds of the parameters

    if (is.null(lower)) {
      lower <- copula@param.lowbnd[par_num]
    }
    else{
      par_num <-
        tryCatch(
          param_num_checked(copula, param_val = lower, param_name = param_name),
          error = function(e) {
            cylcop::error_sound()
            cat("among the 'lower' values: ")
            rlang::abort(conditionMessage(e))
          }
        )
      param_num_checked(copula, param_val = lower, param_name = param_name)
    }


# Check the upper bounds of the parameters

    if (is.null(upper)) {
      upper <- copula@param.upbnd[par_num]
    }
    else{
      par_num <-
        tryCatch(
          param_num_checked(copula, param_val = upper, param_name = param_name),
          error = function(e) {
            cylcop::error_sound()
            cat("among the 'upper' values: ")
            rlang::abort(conditionMessage(e))
          }
        )
    }
    return_list <- list(emp_cop, lower, upper)
    return(return_list)
  }



#' Carry out MLE of copula parameters
#'
#' @param copula A \code{cyl_copula} object.
#' @param emp_cop A numeric matrix with 2 columns holding the empirical copula values.
#' @param param_name A vector of character strings holding the names of the parameters to be optimized.
#'   These can be any parameters in \code{copula@@parameters}.
#' @param start A vector of staring values of the parameters.
#' @param lower A vector of lower bounds of the parameters.
#' @param upper A vector of upper bounds of the parameters.
#' @param optim.method optimizer used in \code{optim()}, can be
#'   "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", or "Brent"
#' @param optim.control A list of additional controls passed to \code{optim()}.
#' @param estimate.variance A logical valu, denoting whether to include an
#'   estimate of the variance (NOT YET IMPLEMENTED).
#' @param traceOpt A logical value, denoting whether to print information regarding
#'   convergence, current values, etc. during the optimization process.
#'
#' @return A list containing the same type of \code{cyl_copula} object as \code{copula},
#'  but with optimized parameters, the log likelihood and the AIC.
#' @export
#'
fit_LL <-
  function(copula,
           emp_cop,
           param_name,
           start,
           lower,
           upper,
           optim.method,
           optim.control,
           estimate.variance,
           traceOpt) {

    # Determine optim() inputs
    control <-
      c(as.list(optim.control), fnscale = -1) # fnscale < 0 => maximization

    # turn every occurence of Inf and -Inf in a vector x into largest positive or negative number NA/NaN basically unchanged
    asFinite <- function(x) {
      if (any(nifi <- !is.finite(x)))
        x[nifi] <- sign(x[nifi]) * .Machine$double.xmax
      x
    }

    #if the optimization method is bounded
    if (optim.method %in% c("Brent", "L-BFGS-B")) {
      lower <-
        asFinite(asFinite(lower) + .Machine$double.eps ^ 0.5 * abs(asFinite(lower)))  #lower is slightly raised and made finite
      upper <-
        asFinite(asFinite(upper) - .Machine$double.eps ^ 0.5 * abs(asFinite(upper)))  # upper slightly lowered and made finite
    }

    #Calculate the log likelihood for a certain set of parameters
    logL <- function(param, param_name, emp_cop, copula) {

      #set parameters in copula object

      copula <-
        tryCatch(
          setCopParam(copula, param_val = param, param_name = param_name),
          OOB_too_large = function(e) {
            warning(
              cylcop:: warning_sound(),
              "a parameter has become larger than its upper bound,\nLogLik set to -Inf"
            )
            return(FALSE)
          },
          OOB_too_small = function(e) {
            warning(
              cylcop:: warning_sound(),
              "a parameter has become smaller than its lower bound,\nLogLik set to -Inf"
            )
            return(FALSE)
          }
        )


#calculate log likelihood

      if (is.logical(copula)) {
        LL <- -Inf
      }
      else{
        LL <- sum(log(cylcop::dCopula(emp_cop, copula = copula)))
        if (is.na(LL))
          LL <- -Inf
      }
      #we want finite likelihood for these optimizers
      if (optim.method %in% c("Brent", "L-BFGS-B", "BFGS")) {
        LL <- asFinite(LL)
      }
      if (traceOpt) {
        if (length(param) <= 1)
          cat(sprintf("param=%14.9g => logL=%12.8g\n", param, LL))

        else
          cat(sprintf(
            "param= %s => logL=%12.8g\n",
            paste(format(param, digits = 9), collapse = ", "),
            LL
          ))
      }
      return(LL)
    }


# Maximize the likelihood

    fit <- stats::optim(
      par = start,
      fn = logL,
      lower = lower,
      upper = upper,
      method = optim.method,
      control = control,
      copula = copula,
      emp_cop = emp_cop,
      param_name = param_name
    )

    setCopParam(copula, param_val = fit$par, param_name = param_name)
    loglik <- fit$val


#Check convergence, give warnings

    warn_tips <-
      paste(
        " - Change the starting values. Find appropriate ones using optCor() or optTau().\n",
        "- Change the optimization method 'optim.method='.\n",
        "- Worst case, try a global optimization using annealing (optim.method = \"SANN\").\n",
        "\tFor this to work properly, you need to play around with the starting temperature\n",
        "\tand the number of function evaluations at each temperature,\n",
        "\tparameters 'temp' and 'tmax' (both default 10) in 'optim.control='.\n",
        "- If no convergence with SANN, observe optimization process (traceOpt = T)\n",
        "\tto see in what direction it is going and get an idea for starting values."
      )

    if (fit[["convergence"]] != 0) {
      warning(
        cylcop:: warning_sound(),
        "Possible convergence problem: optim() gave code=",
        fit$convergence,
        "\n",
        warn_tips
      )
    }

    if(cylcop.env$silent==F){
      message("finished, optimized parameters are: ")
    for (i in seq_along(param_name)) {
      message(param_name[i], " = ", fit$par[i])
    }
    message("logL = ", loglik, "\nAIC = ", 2*length(param_name)-2*loglik)
    }
    if (loglik > 10 ^ 10 || loglik < -10^10) {
      warning(
        cylcop:: warning_sound(),
        "Optimization seems to have gone wrong. logL is very large. Try one of these:\n",
        warn_tips
      )
    }

    opt_cop <-
      setCopParam(copula, param_val = fit$par, param_name = param_name)
    return(list(copula=opt_cop,logL=loglik, AIC=(2*length(param_name)-2*loglik)))
  }
