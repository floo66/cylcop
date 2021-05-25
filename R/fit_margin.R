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
#'   should be fitted "vonmises", "wrappedcauchy", "mixedvonmises", or FALSE if a non-parametric
#'   estimation (kernel density) should be made.
#' @param bandwidth If \code{parametric=F} The numeric value of the kernel density bandwidth.
#'  Default is  \code{cylcop::opt_circ_bw(theta, "adhoc")}.
#' @param mu (optional) Fixed mean direction of parametric distribution
#'
#' @return If a parametric estimate is made, a list containing the estimated parameters, their
#' standard errors and the loglikelihood is returned.
#' If a nonparametric estimate is made, the output comes from
#' the funtion \code{circular::density.circular()}.
#' @export
#'
fit_angle <-
  function(theta,
           parametric = c("vonmises", "wrappedcauchy", "mixedvonmises", FALSE),
           bandwidth = NULL, mu=NULL) {
    theta <- na.omit(theta)

    if (!parametric == FALSE) {
     if(!parametric %in% c("vonmises", "wrappedcauchy", "mixedvonmises")){
       stop(cylcop::error_sound(),
            "distribution must be vonmises, wrappedcauchy or mixedvonmises")
     }
     else if (parametric == "vonmises") {
       if(!is.null(mu) && length(mu)!=1){
         stop("If mean direction is not estimated, length of mu should be 1")
       }
        distr <- suppressWarnings(circular::mle.vonmises(
          theta,
          mu = mu,
          kappa = NULL,
          bias = FALSE
        ))
        logL <-
          suppressWarnings(circular::dvonmises(theta, distr$mu, distr$kappa)) %>% log() %>% sum()
        df <- ifelse(is.null(mu),2,1)
        out <-
          list(coef=
            list(mu = distr$mu %>% as.double() ,
            kappa = distr$kappa
            ),
            se=list(mu = distr$se.mu ,
            kappa = distr$se.kappa
            ),
            logL = logL ,
            AIC = 2*df-2*logL ,
            name="vonmises"
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
            logL %>% round(5),
            "\n",
            "AIC:\t",
            2*df-2*logL %>% round(3),
            "\n"
          )
        if(cylcop.env$silent==F){
          message(msg)}
        return(out)
      }

     else if (parametric == "mixedvonmises") {
       distr <- mle.mixedvonmises(theta, mu=mu)

#for compatibility with other distributions calculate logL with circular package
        logL <-suppressWarnings(circular::dmixedvonmises(theta,
                                                         distr$mu[1],
                                                         distr$mu[2],
                                                         distr$kappa[1],
                                                         distr$kappa[2],
                                                         distr$prop) %>%
                                  log() %>% sum())
        df <- ifelse(is.null(mu),5,3)
        out <-
          list(
            coef=list(
            mu1 = distr$mu[1] ,
            mu2 = distr$mu[2] ,
            kappa1 = distr$kappa[1] ,
            kappa2 = distr$kappa[2] ,
            prop = distr$prop
            ),
            se=list(),
            logL = logL ,
            AIC = 2*df-2*logL ,
            name="mixedvonmises"
          )
        msg <-
          paste0(
            "mu1:\t",
            distr$mu[1] %>% round(3),
            "\n",
            "mu2:\t",
            distr$mu[2] %>% round(3),
            "\n",
            "kappa1:\t",
            distr$kappa[1] %>% round(3),
            "\n",
            "kappa2:\t",
            distr$kappa[2] %>% round(3),
            "\n",
            "prop:\t",
            distr$prop %>% round(3),
            "\n",
            "logL:\t",
            logL %>% round(5),
            "\n",
            "AIC:\t",
            2*df-2*logL %>% round(3),
            "\n"
          )
        if(cylcop.env$silent==F){
          message(msg)}
        return(out)
      }

      else if (parametric == "wrappedcauchy") {
        if(!is.null(mu) && length(mu)!=1){
          stop("If mean direction is not estimated, length of mu should be 1")
        }
        distr <- suppressWarnings(circular::mle.wrappedcauchy(
          theta,
          mu = mu,
          rho = NULL,
          tol = 1e-15,
          max.iter = 100
        ))
        # Calculating the density does not need a numerical approximation. Therefore we can use circular::dwrappedcauchy()
        # instead of cylcop::dwrappedcauchy(), but have to convert the parameter rho to scale
        logL <-
          suppressWarnings(circular::dwrappedcauchy(theta, mu = distr$mu, rho = distr$rho)) %>%
          log() %>% sum()
        df <- ifelse(is.null(mu),2,1)
        #parameterization used in wprappedg("cauchy"), is scale=-ln(rho)
        out <-
          list(
            coef=list(
            location = distr$mu %>% as.double() ,
            scale = -log(distr$rho)
            ),
            logL = logL ,
            AIC = 2*df-2*logL ,
            name="wrappedcauchy"
          )
        msg <- paste0(
          "location:\t",
          distr$mu %>% round(3),
          "\n",
          "scale:\t\t",
          -log(distr$rho) %>% round(3),
          "\n",
          "logL:\t\t",
          logL %>% round(5),
          "\n",
          "AIC:\t",
          2*df-2*logL %>% round(3),
          "\n"
        )
        if(cylcop.env$silent==F){
          message(msg)}
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
#' the function \code{GoFKernel::density.reflected()}.
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
    if(!parametric %in% c("beta", "cauchy", "chi-squared",
                          "exponential", "gamma",
                          "lognormal", "logistic",
                          "normal","t", "weibull")){
      stop(cylcop::error_sound(),
           "distribution not supported")
    }
    distr <- MASS::fitdistr(x, densfun = parametric, start = start)
    out <-
      list(
        coef=as.list(distr$estimate),
        se=as.list(distr$sd),
        logL=distr$loglik ,
        AIC=2*attributes(logLik(distr))$df-2*distr$loglik,
        name=parametric
      )
    if(cylcop.env$silent==F){
      message(show(distr),"\nlogL= ", distr$loglik, "\n",
            "AIC= ", 2*attributes(logLik(distr))$df-2*logLik(distr), "\n")}
    return(out)
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
