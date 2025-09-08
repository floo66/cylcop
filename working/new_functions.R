#this replaces fit_steplength


fit_lin_param <-   function(x,
                            densfun = c(
                              "beta",
                              "cauchy",
                              "chi-squared",
                              "chisq",
                              "exponential",
                              "exp",
                              "gamma",
                              "lognormal",
                              "lnorm",
                              "lognorm",
                              "logistic",
                              "normal",
                              "t",
                              "weibull",
                              "normalmix",
                              "weibullmix",
                              "gammamix",
                              "lnormmix"),
                            start = NULL,
                            ncomp = 2) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(densfun,
                                           type="character",
                                           values = c(
                                             "beta",
                                             "cauchy",
                                             "chi-squared",
                                             "chisq",
                                             "exponential",
                                             "exp",
                                             "gamma",
                                             "lognormal",
                                             "lnorm",
                                             "lognorm",
                                             "logistic",
                                             "normal",
                                             "t",
                                             "weibull",
                                             "normalmix",
                                             "weibullmix",
                                             "gammamix",
                                             "lnormmix"),
                                           length=1)
                  ,1)

    check_arg_all(list(check_argument_type(start,
                                           type="NULL"),
                       check_argument_type(start,
                                           type="list"))
                  ,2)
    check_arg_all(check_argument_type(ncomp,
                                      type="numeric",
                                      length=1,
                                      integer=T,
                                      lower=1)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  x <- na.omit(x)
  if (!is.logical(densfun))
    densfun <- match.arg(densfun)
  if (densfun == "lnorm"||densfun == "lognorm")
    densfun <- "lognormal"
  if (densfun == "chisq")
    densfun <- "chi-squared"
  if (densfun == "exp")
    densfun <- "exponential"

    if(!densfun %in% c("normal",
                       "t",
                          "logistic",
                          "normalmix",
                          "cauchy")){
      if(any(x <= 0)){
        stop("With the selected distribution, all entries of x must be larger than 0.")
      }
    }
    if (densfun == "chi-squared") {
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
    if(densfun == "beta"){
      if (is.null(start)) {
        stop(
          error_sound(),
          "For a beta distribution, you have to provide start values for shape1
            and shape2"
        )}
      distr <-
        MASS::fitdistr(x,
                       densfun = "beta",
                       start = start,
                       method = "BFGS")
    }

    if (!densfun %in% c("chi-squared",
                           "beta",
                           "normalmix",
                           "weibullmix",
                           "gammamix",
                           "lnormmix")) {
      distr <-
        MASS::fitdistr(x,
                       densfun = densfun,
                       start = start)
    }
    if (!densfun %in% c("normalmix",
                           "weibullmix",
                           "gammamix",
                           "lnormmix")) {

      if (densfun == "lognormal")
        densfun <- "lnorm"
      if (densfun == "exponential")
        densfun <- "exp"
      if (densfun == "chi-squared")
        densfun <- "chisq"


      out <-
        list(
          coef = as.list(distr$estimate),
          se = as.list(distr$sd),
          logL = distr$loglik ,
          AIC = 2 * attributes(logLik(distr))$df - 2 * distr$loglik,
          name = densfun
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
        densfun,
        "normalmix" = "normal",
        "weibullmix" = "weibull",
        "gammamix" = "gamma",
        "lnormmix" = "lnorm"
      )
      distr <- mixR::mixfit(x, ncomp = ncomp, family = family)

      if (densfun == "normalmix") {
        coef <- list(prop=distr[[1]],
                     mu=distr[[2]],
                     sd=distr[[3]])
      } else if (densfun == "gammamix") {
        coef <- list(prop=distr[[1]],
                     shape=distr[[4]],
                     rate=distr[[5]])
      }
      else if(densfun == "lnormmix"){
        coef <- list(prop=distr[[1]],
                     meanlog=distr[[4]],
                     sdlog=distr[[5]])
      }
      else if(densfun == "weibullmix"){
        coef <- list(prop=distr[[1]],
                     shape=distr[[4]],
                     scale=distr[[5]])
      }
        out <-
          list(
            coef = coef,
            se = list(),
            logL = distr$loglik,
            AIC = distr$aic,
            name = densfun
          )

      if (cylcop.env$silent == F) {
        message(show(distr))
      }
    }
    return(out)
  }






fit_lin_np <-
  function(x,
           bandwidth = NULL,
           limits=c(-Inf,Inf)) {
    #validate input
    tryCatch({
      check_arg_all(check_argument_type(x,
                                        type="numeric")
                    ,1)

      check_arg_all(list(check_argument_type(bandwidth,
                                             type="NULL"),
                         check_argument_type(bandwidth,
                                             type="numeric",
                                             lower=0))
                    ,2)
      check_arg_all(check_argument_type(limits,
                                        type="numeric",
                                        length=2)
                    ,1)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )
    x <- na.omit(x)

    # Check if every value in x is within the specified limits
    if (any(x < limits[1] | x > limits[2])) {
      stop(
        paste0(
          "Error: All values in x must be between ",limits[1]," and ",limits[2],". Found values outside this range.",
          limits[1], limits[2]
        )
      )
    }


    if (is.null(bandwidth)) {
      bandwidth <- opt_lin_bw(x, "nrd")
    }

      if(limits[1] > -Inf){
        from <- limits[1]
        lower <- limits[1]
      }else{
        from <- range(x)[1]-3*bandwidth
        lower <- -Inf
      }

    if(limits[2] < Inf){
      to <- limits[2]
      upper <- limits[2]
    }else{
      to <- range(x)[2]+3*bandwidth
      upper <- Inf
    }

      #Use density.reflected instead of density to avoid boundary effects
      dens <-
        GoFKernel::density.reflected(
          x,
          kernel = "gaussian",
          from = from,
          to = to,
          lower = lower,
          upper = upper,
          bw = bandwidth,
          n = 4096
        )
      dens[["kernel"]] <- "norm"
      return(dens)
  }

fit_circ_param <-
  function(theta,
           densfun = c("vonmises", "wrappedcauchy", "vonmisesmix"),
           mu = NULL,
           ncomp = 2) {

    #validate input
    tryCatch({
      check_arg_all(check_argument_type(theta,
                                        type="numeric")
                    ,1)
      check_arg_all(check_argument_type(densfun,
                                             type="character",
                                             values = c("vonmises", "wrappedcauchy", "vonmisesmix"))
                    ,1)
      check_arg_all(check_argument_type(ncomp,
                                        type="numeric",
                                        length=1,
                                        integer=T,
                                        lower=1)
                    ,1)
      mu_length <- 1
      if(densfun=="vonmisesmix"){
        mu_length <- ncomp
      }

        check_arg_all(list(check_argument_type(mu,
                                               type="NULL"),
                           check_argument_type(mu,
                                               type="numeric",
                                               length=mu_length))
                      ,2)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )
    theta <- na.omit(theta)
    theta <- full2half_circ(theta)

      if (densfun == "vonmises") {
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

      else if (densfun == "vonmisesmix") {
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

      else if (densfun == "wrappedcauchy") {
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
        #parameterization used in dwprappedcauchy(), is scale=-ln(rho)
        distr$scale <- -log(distr$rho)
        distr$mu <- as.double(distr$mu)
        logL <-
          dwrappedcauchy(theta, location = distr$mu, scale = distr$scale) %>%
          log() %>% sum()
        df <- ifelse(is.null(mu), 2, 1)

        out <-
          list(
            coef = list(
              location = distr$mu ,
              scale = distr$scale
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


fit_circ_np <-
  function(theta,
           bandwidth = NULL) {

    #validate input
    tryCatch({
      check_arg_all(check_argument_type(theta,
                                        type="numeric")
                    ,1)
      check_arg_all(list(check_argument_type(bandwidth,
                                             type="NULL"),
                         check_argument_type(bandwidth,
                                             type="numeric",
                                             lower=0))
                    ,2)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )
    theta <- na.omit(theta)
    theta <- full2half_circ(theta)

    if (is.null(bandwidth)) {
        bandwidth <- opt_circ_bw(theta, "nrd")
      }
      suppressWarnings(
      test <-   circular::density.circular(
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

fit_steplength() <- function(...){
  mc <- list(...)
  if(!mc$paramtric){
    .Deprecated("fit_lin_np", package = "cylcop")
    fit_lin_np(x = mc$x, bandwidth = mc$bandwidth, limits = c(0,Inf))
  }else{
    .Deprecated("fit_lin_param", package = "cylcop")
    fit_lin_param(x = mc$x, densfun = mc$parametric, start = mc$start, ncomp = mc$ncomp)
  }
}

fit_angle() <- function(...){
  mc <- list(...)
  if(!mc$paramtric){
    .Deprecated("fit_circ_np", package = "cylcop")
    fit_circ_np(theta = mc$theta, bandwidth = mc$bandwidth, )
  }else{
    .Deprecated("fit_circ_param", package = "cylcop")
    fit_circ_param(theta = mc$theta, densfun = mc$parametric, mu = mc$mu, ncomp = mc$ncomp)
  }
}




new_pcylcop <- function(u,copula){
# Pre-allocate result vector
  res <- rep(NA, nrow(u))

  # Handle cases where u contains 0 or 1
  zero_ind <- unique(which(u == 0, arr.ind = TRUE)[, 1])
  one_ind <- unique(which(u == 1, arr.ind = TRUE)[, 1])

  # For rows with 1, handle two cases:
  # 1. If a row has both entries as 1, return 1.
  # 2. If a row has only one entry as 1, return the other value.
  if (length(one_ind) > 0) {
    one_value <- apply(u[one_ind, ,drop=F], 1, function(row) {
      if (all(row == 1)) {
        return(1)  # Both columns are 1, so return 1.
      } else {
        return(row[row != 1])  # Return the non-1 value.
      }
    })

    # Ensure it's a vector (not a list or matrix)
    one_value <- unlist(one_value)

    res[one_ind] <- one_value
  }

  # Set values for rows with u == 0
  res[zero_ind] <- 0

  # Handle remaining rows that are neither fully 0 nor 1
  ind_remain <- setdiff(seq_len(nrow(u)), c(zero_ind, one_ind))

  if (length(ind_remain) > 0) {
    res[ind_remain] <- copula::pCopula(u[ind_remain, ], copula)
  }

  return(res)
}
