#' Find the optimal bandwidth for circular kernel density estimate
#'
#' The ad-hoc method of finding the bandwidth parameter might give very bad results, especially for multimodal population distributions.
#' In these cases it might help to set \code{kappa.est=trigmoments}.
#'
#' @param theta A numeric vector of angles in [-pi, pi).
#' @param loss A string of characters describing the loss function, either "KullbackLeibler" (leading to an MLE estimate),
#'     or "adhoc" leading to a rule-of-thumb estimate.
#' @param kappa.est Character string describing how the spread is estimated. Either maximum likelihood "ML"
#'     or trigonometric moment "trigmoments"
#'
#' @return The numeric value of the optimized bandwidth.
#' @export
#'
opt_circ_bw <- function(theta,
                        loss = c("KullbackLeibler", "adhoc"),
                        kappa.est="ML") {
  if (loss == "KullbackLeibler") {
    bw <- suppressWarnings(circular::bw.cv.ml.circular(
      na.omit(theta),
      lower = NULL,
      upper = NULL,
      tol = 1e-4,
      kernel = "vonmises",
      K = NULL,
      min.k = 10
    ))
  }
  else if (loss == "adhoc") {
    warning(
      "ad-hoc method of finding bandwidth parameter might give very bad results, especially for multimodal population distributions\n",
      "you can use(kappa.est=\"trigmoments\" or an mle method (loss=\"KullbackLeibler\")"
    )
    bw <-
      suppressWarnings(circular::bw.nrd.circular(
          na.omit(theta),
          lower = NULL,
          upper = NULL,
          kappa.est = kappa.est,
          kappa.bias = FALSE,
          P = 3
        ))
  }
  else
    stop(
      "method of estimating bandwidth must be from c(\"Kullback-Leibler\",\"adhoc\")"
    )
  return(bw)
}





#' Find the optimal bandwidth for linear kernel density estimate
#'
#' @param x A numeric vector of linear measurements.
#' @param loss A string of characters describing the loss function, either "KullbackLeibler" (leading to an MLE estimate),
#'     or "adhoc" leading to a rule-of-thumb estimate.
#'
#' @return The numeric value of the optimized bandwidth.
#' @export
#'
opt_lin_bw <- function(x,
                       loss = c("KullbackLeibler", "adhoc")) {
  if (loss == "KullbackLeibler") {
    bw <-stats::bw.ucv(x,nb=4000)
  }
  else if (loss == "adhoc") {
    warning(
      cylcop::warning_sound(),
      "ad-hoc method of finding bandwidth parameter might give very bad results, especially for multimodal population distributions\n"
    )
    bw <-stats::bw.nrd(x)
  }
  else
    stop(
      cylcop::error_sound(),
      "method of estimating bandwidth must be from c(\"KullbackLeibler\",\"adhoc\")"
    )
  return(bw)
}



#' Fit circular univariate distribution
#'
#' This function finds parameter estimates of the marginal circular
#' distribution, or gives a kernel density estimate.
#'
#' @param theta A numeric vector of angles in [-pi, pi).
#' @param parametric Either a character string describing what distribution
#'   should be fitted "vonmises" or "wrappedcauchy" or FALSE if a non-parametric
#'   estimation (kernel density) should be made.
#' @param bandwidth The numeric value of for the kernel density bandwidth.
#'  Default is  \code{cylcop::opt_circ_bw(theta, "adhoc")}.
#'
#' @return If a parametric estimate is made, a list containing the estimated parameters, their
#' standard errors and the loglikelihood is returned.
#' If a nonparametric estimate is made, the output comes from
#' the funtion \code{circular::density.circular()}.
#' @export
#'
fit_angle <-
  function(theta,
           parametric = c("vonmises", "wrappedcauchy", FALSE),
           bandwidth = NULL) {
    theta <- na.omit(theta)

    if (!parametric == FALSE) {

      if (parametric == "vonmises") {
        distr <- suppressWarnings(circular::mle.vonmises(
          theta,
          mu = NULL,
          kappa = NULL,
          bias = FALSE
        ))
        logL <-
          suppressWarnings(circular::dvonmises(theta, distr$mu, distr$kappa)) %>% log() %>% sum()
        out <-
          list(
            mu = distr$mu %>% as.double() %>% round(3),
            kappa = distr$kappa %>% round(3),
            se.mu = distr$se.mu %>% round(3),
            se.kappa = distr$se.kappa %>% round(3),
            logL = logL %>% round(5)
          )
        msg <-
          paste0(
            "mu:\t",
            distr$mu %>% round(3),
            " (",
            distr$se.mu %>% round(3),
            ")\n",
            "kappa:\t",
            distr$kappa %>% round(3),
            " (",
            distr$se.kappa %>% round(3),
            ")\n",
            "logL:\t",
            logL %>% round(5)
          )
        cat(msg)
        return(out)
      }

      else if (parametric == "wrappedcauchy") {
        distr <- suppressWarnings(circular::mle.wrappedcauchy(
          theta,
          mu = NULL,
          rho = NULL,
          tol = 1e-15,
          max.iter = 100
        ))
        # Calculating the density does not need a numerical approximation. Therefore we can use circular::dwrappedcauchy()
        # instead of dwrappedcauchy(), but have to convert the parameter rho to scale
        logL <-
          suppressWarnings(circular::dwrappedcauchy(theta, mu = distr$mu, rho = distr$rho)) %>%
          log() %>% sum()
        #parameterization used in wprappedg("cauchy"), is scale=-ln(rho)
        out <-
          list(
            location = distr$mu %>% as.double() %>% round(3),
            scale = -log(distr$rho) %>% round(3),
            logL = logL %>% round(5)
          )
        msg <- paste0(
          "location:\t",
          distr$mu %>% round(3),
          "\n",
          "scale:\t\t",
          -log(distr$rho) %>% round(3),
          "\n",
          "logL:\t\t",
          logL %>% round(5)
        )
        cat(msg)
        return(out)
      }
    }

    else{
      if (is.null(bandwidth)) {
        bandwidth <- opt_circ_bw(theta, "adhoc")
      }
      suppressWarnings(
        circular::density.circular(
          theta,
          bw = bandwidth,
          kernel = "vonmises",
          na.rm = TRUE,
          from = circular::circular(-pi),
          to = circular::circular(pi),
          n = 4096
        )
      )
    }
  }



#' Fit linear univariate distribution
#'
#' This function finds parameter estimates of the marginal linear
#' distribution, or gives a kernel density estimate.
#'
#' @param x A numeric vector of measurements of a linear random variable.
#' @param parametric Either a character string describing what (continuous)
#'   distribution should be fitted or FALSE if a non-parametric estimation
#'   (kernel density) should be made.
#' @param bandwidth The numeric value of for the kernel density bandwidth.
#'  Default is  \code{cylcop::opt_lin_bw(x, "adhoc")}.
#' @param start OPTIONAL: A numeric vector holding the starting values  for the ML parameter estimation.
#'
#' @return If a parametric estimate is made, a list containing the estimated parameters, their
#' standard errors and the loglikelihood is returned.
#' If a nonparametric estimate is made, the output comes from
#' the funtion \code{GoFKernel::density.reflected()}.
#' @export
#'
fit_steplength<- function(x, parametric = c("beta", "cauchy", "chi-squared",
                                            "exponential", "gamma",
                                            "lognormal", "logistic",
                                            "normal","t", "weibull", FALSE),
                          start = NULL, bandwidth = NULL){
  x <- na.omit(x)
  if (!is.logical(parametric))
    parametric <- match.arg(parametric)

  if (!parametric == FALSE) {
    distr <- MASS::fitdistr(x, densfun = parametric, start = start)
    show(distr)
    cat("\nlogL= ", distr$loglik, "\n")
    return(distr)
  }

  else{
    if (is.null(bandwidth)) {
      bandwidth <- opt_lin_bw(x, "adhoc")
    }

    #Use density.reflected instead of density to avoid boundary effects
    dens<-
      GoFKernel::density.reflected(
      x,
      kernel = "gaussian",
      from=0,
      lower = 0,
      upper = Inf,
      bw = bandwidth,
      n = 4096
    )
    dens[["kernel"]] <- "norm"
    return(dens)
  }
}
