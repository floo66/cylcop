#' Find the Optimal Bandwidth for a Circular Kernel Density Estimate
#'
#' This function basically wraps \code{circular::\link[circular]{bw.cv.ml.circular}()}
#' and \code{circular::\link[circular]{bw.nrd.circular}()} of the '\pkg{circular}'
#'  package, simplifying their inputs. For more control,
#' these '\pkg{circular}' functions could be used directly.
#' The ad-hoc method of finding the bandwidth parameter might give very bad results,
#'  especially for multimodal population distributions.
#' In these cases it can help to set \code{kappa.est = trigmoments}.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles in \eqn{[-\pi, \pi)}.
#' @param loss \link[base]{character} string describing the loss function,
#' either \code{"KullbackLeibler"} (leading to a maximum likelihood estimate), or \code{"adhoc"}
#'  leading to a rule-of-thumb estimate.
#' @param kappa.est \link[base]{character} string describing how the spread is estimated.
#' Either maximum likelihood \code{"ML"}, or trigonometric moment \code{"trigmoments"}.
#'
#' @return A \link[base]{numeric} value, the optimized bandwidth.
#'
#' @examples require(circular)
#' require(graphics)
#' set.seed(123)
#' n <- 100  #n (number of samples) is set small for performance. Increase n to
#' # a value larger than 1000 to see the effects of multimodality
#'
#' angles <- rvonmisesmix(n,
#'   mu = c(0,pi),
#'   kappa = c(2,1),
#'   prop = c(0.5,0.5)
#' )
#' bw1 <- opt_circ_bw(theta = angles, loss="adhoc")
#' bw2 <- opt_circ_bw(theta = angles, loss="adhoc", kappa.est = "trigmoments")
#' bw3 <- opt_circ_bw(theta = angles, loss="KullbackLeibler")
#'
#' dens1 <- fit_angle(theta = angles, parametric = FALSE, bandwidth = bw1)
#' dens2 <- fit_angle(theta = angles, parametric = FALSE, bandwidth = bw2)
#' dens3 <- fit_angle(theta = angles, parametric = FALSE, bandwidth = bw3)
#' true_dens <- dvonmisesmix(
#'   seq(-pi,pi,0.001),
#'   mu = c(0,pi),
#'   kappa = c(2,1),
#'   prop = c(0.5,0.5)
#' )
#'
#' plot(seq(-pi, pi, 0.001), true_dens, type = "l")
#' lines(as.double(dens1$x), as.double(dens1$y), col = "red")
#' lines(as.double(dens2$x), as.double(dens2$y), col = "green")
#' lines(as.double(dens3$x), as.double(dens3$y), col = "blue")
#'
#' @seealso \code{circular::\link[circular]{bw.cv.ml.circular}()},
#' \code{circular::\link[circular]{bw.nrd.circular}()},
#' \code{\link{opt_circ_bw}()}.
#'
#' @export
#'
opt_circ_bw <- function(theta,
                        loss = c("KullbackLeibler", "adhoc"),
                        kappa.est = "ML") {
  if (loss == "KullbackLeibler") {
    bw <- suppressWarnings(
      circular::bw.cv.ml.circular(
        na.omit(theta),
        lower = NULL,
        upper = NULL,
        tol = 1e-4,
        kernel = "vonmises",
        K = NULL,
        min.k = 10
      )
    )
  }
  else if (loss == "adhoc") {
    if (kappa.est != "trigmoments") {
      warning(
        "ad-hoc method of finding bandwidth parameter might give very bad results, especially for multimodal population distributions\n",
        "you can use(kappa.est=\"trigmoments\" or an mle method (loss=\"KullbackLeibler\")"
      )
    }
    bw <-
      suppressWarnings(
        circular::bw.nrd.circular(
          na.omit(theta),
          lower = NULL,
          upper = NULL,
          kappa.est = kappa.est,
          kappa.bias = FALSE,
          P = 3
        )
      )
  }
  else
    stop("method of estimating bandwidth must be from c(\"Kullback-Leibler\",\"adhoc\")")
  return(bw)
}





#' Find the Optimal Bandwidth for a Linear Kernel Density Estimate
#'
#' This function is basically wraps \code{stats::\link[stats]{bw.ucv}()}
#' and \code{stats::\link[stats]{bw.nrd}()} of the '\pkg{stats}'
#'  package, simplifying their inputs. For more control,
#' these '\pkg{stats}' functions could be used directly.
#'
#' @param x \link[base]{numeric} \link[base]{vector} of linear measurements.
#' @param loss \link[base]{character} string describing the loss function,
#' either \code{"KullbackLeibler"} (leading to an maximum likelihood estimate),
#'     or \code{"adhoc"} leading to a rule-of-thumb estimate.
#'
#' @return A \link[base]{numeric} value, the optimized bandwidth.
#'
#' @examples require(graphics)
#' set.seed(123)
#' n <- 1000
#'
#' x <- rweibull(n, shape = 10)
#' bw1 <- opt_lin_bw(x = x, loss="adhoc")
#' bw2 <- opt_lin_bw(x = x, loss="KullbackLeibler")
#'
#' dens1 <- fit_steplength(x = x, parametric = FALSE, bandwidth = bw1)
#' dens2 <- fit_steplength(x = x, parametric = FALSE, bandwidth = bw2)
#' true_dens <- dweibull(seq(0,max(x),length.out = 200), shape = 10)
#'
#' plot(seq(0,max(x),length.out = 200), true_dens, type = "l")
#' lines(dens1$x, dens1$y, col = "red")
#' lines(dens2$x, dens2$y, col = "green")
#'
#' @seealso \code{stats::\link[stats]{bw.ucv}()},
#' \code{stats::\link[stats]{bw.nrd}()}
#' \code{\link{opt_lin_bw}()}.
#'
#' @export
#'
opt_lin_bw <- function(x,
                       loss = c("KullbackLeibler", "adhoc")) {
  if (loss == "KullbackLeibler") {
    bw <- stats::bw.ucv(x, nb = 10000)
  }
  else if (loss == "adhoc") {
    bw <- stats::bw.nrd(x)
  }
  else
    stop(
      error_sound(),
      "method of estimating bandwidth must be from c(\"KullbackLeibler\",\"adhoc\")"
    )
  return(bw)
}



#' Fit a Circular Univariate Distribution
#'
#' This function finds parameter estimates of the marginal circular
#' distribution (with potentially fixed mean), or gives a kernel density estimate using
#' a von Mises smoothing kernel.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles in \eqn{[-\pi, \pi)}.
#' @param parametric either a \link[base]{character} string describing what distribution
#'   should be fitted (\code{"vonmises"}, \code{"wrappedcauchy"}, or
#'   \code{"vonmisesmix"}), or the \link[base]{logical} \code{FALSE} if a non-parametric
#'   estimation (kernel density) should be made.
#' @param bandwidth If \code{parametric = FALSE}, the numeric value of the kernel density bandwidth.
#'  Default is \code{cylcop::\link{opt_circ_bw}(theta, "adhoc")}.
#' @param mu (optional) \link[base]{numeric} \link[base]{vector}, fixed mean direction(s) of the
#' parametric distribution.
#' @param ncomp \link[base]{integer}, number of components of the mixed vonMises distribution.
#' Only has an effect if \code{parametric="vonmisesmix"}.
#'
#' @return If a parametric estimate is made, a \link[base]{list}
#'  containing the estimated parameters, their
#' standard errors and the log-likelihood is returned.
#' If a non-parametric estimate is made, the output is a list and comes directly
#' from the function \code{circular::\link[circular]{density.circular}()} of the '\pkg{circular}'
#' package.
#'
#' @examples require(circular)
#' require(graphics)
#' set.seed(123)
#'
#' silent_curr <- cylcop_get_option("silent")
#' cylcop_set_option(silent = TRUE)
#'
#' n <- 100  #n (number of samples) is set small for performance. Increase n to
#' # a value larger than 1000 to see the effects of multimodality
#'
#' angles <- rvonmisesmix(n,
#'   mu = c(0, pi),
#'   kappa = c(2,1),
#'   prop = c(0.5, 0.5)
#' )
#' angles <- as.double(angles)
#'
#' bw <- opt_circ_bw(theta = angles,
#'   loss="adhoc",
#'   kappa.est = "trigmoments"
#' )
#' dens_non_param <- fit_angle(theta = angles,
#'   parametric = FALSE,
#'   bandwidth = bw
#' )
#'
#' param_estimate <- fit_angle(theta = angles,
#'   parametric = "vonmisesmix"
#' )
#' param_estimate_fixed_mean <- fit_angle(theta = angles,
#'   parametric = "vonmisesmix",
#'   mu = c(0, pi)
#' )
#'
#' true_dens <- dvonmisesmix(seq(-pi,pi,0.001),
#'   mu = c(0, pi),
#'   kappa = c(2,1),
#'   prop = c(0.5, 0.5)
#' )
#' dens_estimate <- dvonmisesmix(seq(-pi,pi,0.001),
#'   mu= param_estimate$coef$mu,
#'   kappa = param_estimate$coef$kappa,
#'   prop = param_estimate$coef$prop
#' )
#' dens_estimate_fixed_mean <- dvonmisesmix(seq(-pi,pi,0.001),
#'   mu= param_estimate_fixed_mean$coef$mu,
#'   kappa = param_estimate_fixed_mean$coef$kappa,
#'   prop = param_estimate_fixed_mean$coef$prop
#' )
#'
#' plot(seq(-pi, pi, 0.001), true_dens, type = "l")
#' lines(as.double(dens_non_param$x), as.double(dens_non_param$y), col = "red")
#' lines(seq(-pi, pi, 0.001), dens_estimate, , col = "green")
#' lines(seq(-pi, pi, 0.001), dens_estimate_fixed_mean, , col = "blue")
#'
#' cylcop_set_option(silent = silent_curr)
#'
#' @seealso \code{circular::\link[circular]{density.circular}()},
#' \code{\link{fit_angle}()}, \code{\link{opt_circ_bw}()}.
#'
#'
#' @export
#'
#'
fit_angle <-
  function(theta,
           parametric = c("vonmises", "wrappedcauchy", "vonmisesmix", FALSE),
           bandwidth = NULL,
           mu = NULL,
           ncomp = 2) {
    theta <- na.omit(theta)

    if (!parametric == FALSE) {
      if (!parametric %in% c("vonmises", "wrappedcauchy", "vonmisesmix")) {
        stop(error_sound(),
             "distribution must be vonmises, wrappedcauchy or vonmisesmix")
      }
      else if (parametric == "vonmises") {
        if (!is.null(mu) && length(mu) != 1) {
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
        df <- ifelse(is.null(mu), 2, 1)
        out <-
          list(
            coef =
              list(
                mu = distr$mu %>% as.double() ,
                kappa = distr$kappa
              ),
            se = list(mu = distr$se.mu ,
                      kappa = distr$se.kappa),
            logL = logL ,
            AIC = 2 * df - 2 * logL ,
            name = "vonmises"
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
            2 * df - 2 * logL %>% round(3),
            "\n"
          )
        if (cylcop.env$silent == F) {
          message(msg)
        }
        return(out)
      }

      else if (parametric == "vonmisesmix") {
        distr <- mle.vonmisesmix(theta, mu = mu, ncomp = ncomp)
        logL <- suppressWarnings(
          dvonmisesmix(
            theta,
            mu = distr$mu,
            kappa = distr$kappa,
            prop = distr$prop
          ) %>%
            log() %>% sum()
        )
        df <- ifelse(is.null(mu), 5, 3)
        out <-
          list(
            coef = list(
              mu = distr$mu ,
              kappa = distr$kappa ,
              prop = distr$prop
            ),
            se = list(),
            logL = logL ,
            AIC = 2 * df - 2 * logL ,
            name = "vonmisesmix"
          )
        msg <-
          paste0(
            "mu:\t",
            paste(distr$mu %>% round(3), collapse = ", "),
            "\n",
            "kappa:\t",
            paste(distr$kappa %>% round(3), collapse = ", "),
            "\n",
            "prop:\t",
            paste(distr$prop %>% round(3), collapse = ", "),
            "\n",
            "logL:\t",
            logL %>% round(5),
            "\n",
            "AIC:\t",
            2 * df - 2 * logL %>% round(3),
            "\n"
          )
        if (cylcop.env$silent == F) {
          message(msg)
        }
        return(out)
      }

      else if (parametric == "wrappedcauchy") {
        if (!is.null(mu) && length(mu) != 1) {
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
        df <- ifelse(is.null(mu), 2, 1)
        #parameterization used in wprappedg("cauchy"), is scale=-ln(rho)
        out <-
          list(
            coef = list(
              location = distr$mu %>% as.double() ,
              scale = -log(distr$rho)
            ),
            logL = logL ,
            AIC = 2 * df - 2 * logL ,
            name = "wrappedcauchy"
          )
        msg <- paste0(
          "location:\t",
          distr$mu %>% round(3),
          "\n",
          "scale:\t\t",-log(distr$rho) %>% round(3),
          "\n",
          "logL:\t\t",
          logL %>% round(5),
          "\n",
          "AIC:\t",
          2 * df - 2 * logL %>% round(3),
          "\n"
        )
        if (cylcop.env$silent == F) {
          message(msg)
        }
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



#' Fit a Linear Univariate Distribution
#'
#' This function finds parameter estimates of the marginal linear
#' distribution, or gives a kernel density estimate using
#' a Gaussian smoothing kernel.
#'
#' @param x \link[base]{numeric} \link[base]{vector} of measurements of a linear
#' random variable in \eqn{[0,\infty)}.
#' @param parametric either a \link[base]{character} string describing what distribution
#'  should be fitted (\code{"beta"}, \code{"cauchy"}, \code{"chi-squared"},
#'  \code{"exponential"}, \code{"gamma"}, \code{"lognormal"}, \code{"logistic"},
#'  \code{"normal"}, \code{"t"}, or \code{"weibull"}), or the \link[base]{logical}
#'  \code{FALSE} if a non-parametric estimation (kernel density) should be made.
#' @param bandwidth \link[base]{numeric} value for the kernel density bandwidth.
#'  Default is  \code{cylcop::\link{opt_lin_bw}(x, "adhoc")}.
#' @param start (optional, except when \code{parametric = "chi-squared"})
#' named \link[base]{list} containing the parameters to be optimized with inital
#' values.
#' @param ncomp \link[base]{integer}, number of components of the mixed distribution.
#' Only has an effect if \code{parametric %in% c("normalmix", "weibullmix", "gammamix",
#' "lnormmix")}.
#'
#' @return If a parametric estimate is made, a \link[base]{list} containing the
#' estimated parameters, their standard errors and the log-likelihood is returned.
#' If a non-parametric estimate is made, the output is a list, and comes directly
#' from the function \code{GoFKernel::\link[GoFKernel]{density.reflected}()} of the '\pkg{GoFKernel}'
#' package.
#'
#' @examples require(graphics)
#' set.seed(123)
#'
#' silent_curr <- cylcop_get_option("silent")
#' cylcop_set_option(silent = TRUE)
#'
#' n <- 10000
#'
#' x <- rweibull(n, shape = 10)
#'
#' dens_non_param <- fit_steplength(x = x, parametric = FALSE)
#' weibull <- fit_steplength(x = x, parametric = "weibull")
#' gamma <- fit_steplength(x = x, parametric = "gamma")
#' chisq <- fit_steplength(x = x, parametric = "chi-squared", start = list(df = 1))
#'
#' true_dens <- dweibull(seq(0, max(x), length.out = 200),
#'   shape = 10
#' )
#' dens_weibull <- dweibull(seq(0, max(x),length.out = 200),
#'   shape = weibull$coef$shape,
#'   scale = weibull$coef$scale
#' )
#' dens_gamma <- dgamma(seq(0, max(x),length.out = 200),
#'   shape = gamma$coef$shape,
#'   rate = gamma$coef$rate
#' )
#' dens_chisq <- dchisq(seq(0, max(x),length.out = 200),
#'   df = chisq$coef$df
#' )
#'
#' plot(seq(0,max(x),length.out = 200), true_dens, type = "l")
#' lines(dens_non_param$x, dens_non_param$y, col = "red")
#' lines(seq(0,max(x),length.out = 200), dens_weibull, col = "green")
#' lines(seq(0,max(x),length.out = 200), dens_gamma, col = "blue")
#' lines(seq(0,max(x),length.out = 200), dens_chisq, col = "cyan")
#'
#' cylcop_set_option(silent = silent_curr)
#'
#' @seealso \code{GoFKernel::\link[GoFKernel]{density.reflected}()},
#' \code{\link{fit_angle}()}, \code{\link{opt_lin_bw}()}.
#'
#' @export
#'
fit_steplength <-
  function(x,
           parametric = c(
             "beta",
             "cauchy",
             "chi-squared",
             "exponential",
             "gamma",
             "lognormal",
             "logistic",
             "normal",
             "t",
             "weibull",
             "normalmix",
             "weibullmix",
             "gammamix",
             "lnormmix",
             FALSE
           ),
           start = NULL,
           bandwidth = NULL,
           ncomp = 2) {
    x <- na.omit(x)
    if (!is.logical(parametric))
      parametric <- match.arg(parametric)

    if (!parametric == FALSE) {
      if (!parametric %in% c(
        "cauchy",
        "chi-squared",
        "exponential",
        "gamma",
        "lognormal",
        "logistic",
        "normal",
        "t",
        "weibull",
        "normalmix",
        "weibullmix",
        "gammamix",
        "lnormmix"
      )) {
        stop(error_sound(),
             "distribution not supported")
      }
      if (parametric == "chi-squared") {
        if (is.null(start)) {
          stop(
            error_sound(),
            "For a chi-squared distribution, you have to provide start values for df"
          )
        }
        distr <-
          MASS::fitdistr(x,
                         densfun = "chi-squared",
                         start = start,
                         method = "BFGS")
      }

      low <- switch(
        parametric,
        "cauchy" = c(0, 0),
        "gamma" = c(0, 0),
        "logistic" = c(0, 0),
        "t" = c(0, 0, 0),
        "weibull" = c(0, 0),
        NULL
      )
      if (!parametric %in% c("chi-squared",
                             "normalmix",
                             "weibullmix",
                             "gammamix",
                             "lnormmix")) {
        distr <-
          MASS::fitdistr(x,
                         densfun = parametric,
                         start = start,
                         lower = low)
      }
      if (!parametric %in% c("normalmix",
                             "weibullmix",
                             "gammamix",
                             "lnormmix")) {
        if (parametric == "lognormal")
          parametric <- "lnorm"
        out <-
          list(
            coef = as.list(distr$estimate),
            se = as.list(distr$sd),
            logL = distr$loglik ,
            AIC = 2 * attributes(logLik(distr))$df - 2 * distr$loglik,
            name = parametric
          )
        if (cylcop.env$silent == F) {
          message(
            show(distr),
            "\nlogL= ",
            distr$loglik,
            "\n",
            "AIC= ",
            2 * attributes(logLik(distr))$df - 2 * logLik(distr),
            "\n"
          )
        }
      } else{
        family <- switch(
          parametric,
          "normalmix" = "normal",
          "weibullmix" = "weibull",
          "gammamix" = "gamma",
          "lnormmix" = "lnorm"
        )
        distr <- mixR::mixfit(x, ncomp = ncomp, family = family)
        if (parametric == "normalmix") {
          out <-
            list(
              coef = list(distr[[1]], distr[[2]], distr[[3]]),
              se = list(),
              logL = distr$loglik,
              AIC = distr$aic,
              name = parametric
            )
          names(out$coef) <- c("prop", "mu", "sd")
        } else if (parametric == "gammamix") {
          out <-
            list(
              coef = list(distr[[1]], distr[[4]], 1/distr[[5]]),
              se = list(),
              logL = distr$loglik,
              AIC = distr$aic,
              name = parametric
            )
          names(out$coef) <- c("prop", "shape", "scale")
        }
        else{
          out <-
            list(
              coef = list(distr[[1]], distr[[4]], distr[[5]]),
              se = list(),
              logL = distr$loglik,
              AIC = distr$aic,
              name = parametric
            )
          names(out$coef) <- c("prop", names(distr)[4], names(distr)[5])
          if (parametric == "weibullmix") {
            names(out$coef) <- c("prop", "shape", "scale")
          }
        }

        if (cylcop.env$silent == F) {
          message(show(distr))
        }
      }
      return(out)
    }

    else{
      if (is.null(bandwidth)) {
        bandwidth <- opt_lin_bw(x, "adhoc")
      }

      #Use density.reflected instead of density to avoid boundary effects
      dens <-
        GoFKernel::density.reflected(
          x,
          kernel = "gaussian",
          from = 0,
          lower = 0,
          upper = Inf,
          bw = bandwidth,
          n = 4096
        )
      dens[["kernel"]] <- "norm"
      return(dens)
    }
  }
