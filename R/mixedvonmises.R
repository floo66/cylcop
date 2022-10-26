#' Density, Distribution, Quantiles and Random Number Generation for the mixed von
#' Mises Distribution
#'
#' The number of components in the mixed von Mises distribution is specified by the length
#' of the parameter vectors. The quantiles are numerically obtained from the distribution function using
#' monotone cubic splines.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} giving the angles where
#' the density or distribution function is evaluated.
#' @param p \link[base]{numeric} \link[base]{vector} giving the probabilities where
#' the quantile function is evaluated.
#' @param n \link[base]{integer} value, the number of random samples to be
#' generated with \code{rvonmisesmix()}.
#' @param mu \link[base]{numeric} \link[base]{vector} holding the mean directions.
#' @param kappa \link[base]{numeric} \link[base]{vector} holding the concentration
#' parameters.
#' @param prop \link[base]{numeric} \link[base]{vector}, holding the mixing proportions
#' of the components.
#'
#' @return
#' \itemize{
#' \item{\code{dvonmisesmix()}}{ gives a \link[base]{vector} of length \code{length(theta)}
#'  containing the density at \code{theta}.}
#' \item{\code{pvonmisesmix()}}{ gives a
#' \link[base]{vector} of length \code{length(theta)} containing
#' the distribution function at the corresponding values of \code{theta}.}
#' \item{\code{qvonmisesmix()}}{ gives a \link[base]{vector} of length \code{length(p)}
#' containing the quantiles at the corresponding values of \code{p}.}
#' \item{\code{rvonmisesmix()}}{ generates a \link[base]{vector} of length \code{n}
#' containing the random samples, i.e. angles in \eqn{[-\pi, \pi)}.}
#'}
#'
#' @examples
#'
#' rvonmisesmix(10, mu = c(0, pi, pi/2), kappa = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' dvonmisesmix(c(0, 2, pi, 1), mu = c(0, pi), kappa = c(2, 2), prop = c(0.6, 0.4))
#'
#' prob <- pvonmisesmix(c(0.1, pi), mu = c(0, pi, pi/2), kappa = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#' prob
#' qvonmisesmix(prob, mu = c(0, pi, pi/2), kappa = c(2, 2, 4), prop = c(0.6, 0.3, 0.1))
#'
#' @name vonmisesmix
#'
#'
#' @aliases rvonmisesmix pvonmisesmix dvonmisesmix qvonmisesmix
#'
NULL












# Random numbers
#'
#' @rdname vonmisesmix
#' @export
#'
rvonmisesmix <- function(n, mu, kappa, prop) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(n,
                                      type="numeric",
                                      length = 1,
                                      integer=T,
                                      lower=1)
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(kappa,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  if (!is.numeric(prop) || !is.numeric(mu) || !is.numeric(kappa)) {
    stop("prop, mu, and kappa must be numeric vectors")
  }

  if (length(prop) != length(mu) ||
      length(prop) != length(kappa) || length(mu) != length(kappa)) {
    stop("prop, mu, and kappa must have the same length")
  }
  mu <- as.vector(mu, mode = "numeric")
  kappa <- as.vector(kappa, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)

  dist_component <-
    sample(1:length(prop),
           size = n,
           replace = TRUE,
           prob = prop)
  draws_per_component <-
    tabulate(dist_component, nbins = length(prop))
  x <- rep(0, n)
  for (i in 1:length(prop)) {
    x[dist_component == i] <-
      suppressWarnings(
        circular::rvonmises(
          draws_per_component[i],
          mu = circular::circular(mu[i]),
          kappa = kappa[i]
        )
      ) %>%full2half_circ()
  }

  return(x)
}


# Density
#'
#' @rdname vonmisesmix
#' @export
#'
dvonmisesmix <- function(theta, mu, kappa, prop) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(theta,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(kappa,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  if (length(prop) != length(mu) ||
      length(prop) != length(kappa) || length(mu) != length(kappa)) {
    stop("prop, mu, and kappa must have the same length")
  }
  mu <- as.vector(mu, mode = "numeric")
  kappa <- as.vector(kappa, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)

  out <- 0
  for (i in 1:length(prop)) {
    out <-
      out + prop[i] * suppressWarnings(
        circular::dvonmises(
          x = circular::circular(theta),
          mu = circular::circular(mu[i]),
          kappa = kappa[i]
        )
      )
  }
  out
}


# Distribution
#'
#' @rdname vonmisesmix
#' @export
#'
pvonmisesmix <- function(theta, mu, kappa, prop) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(theta,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(kappa,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  if (length(prop) != length(mu) ||
      length(prop) != length(kappa) || length(mu) != length(kappa)) {
    stop("prop, mu, and kappa must have the same length")
  }
  mu <- as.vector(mu, mode = "numeric")
  kappa <- as.vector(kappa, mode = "numeric")
  prop <- as.vector(prop, mode = "numeric")
  prop <- prop / sum(prop)

  out <- 0
  for (i in 1:length(prop)) {
    out <-
      out + prop[i] * suppressWarnings(
        circular::pvonmises(
          q = circular::circular(theta),
          mu = circular::circular(mu[i]),
          kappa = kappa[i],
          from = circular::circular(-pi)
        )
      )
  }
  out
}




# Quantiles of the Mixed von Mises Distribution
#'
#' @rdname vonmisesmix
#' @export
#'
qvonmisesmix <- function(p, mu, kappa, prop) {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(p,
                                      type="numeric",
                                      lower=0,
                                      upper=1)
                  ,1)
    check_arg_all(check_argument_type(mu,
                                      type="numeric")
                  ,1)
    check_arg_all(check_argument_type(kappa,
                                      type="numeric",
                                      lower=0)
                  ,1)
    check_arg_all(check_argument_type(prop,
                                      type="numeric",
                                      lower=0)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  if (all(p > 0.999999))
    return(rep(pi, length(p)))

  #Numerically find the inverse. This can get pretty slow, because the entire calculation is repeated
  #every time qvonmisesmix() is called
  # fun <- function(q) pvonmisesmix(q, mu, kappa, prop) - p
  # result <- suppressWarnings(uniroot(fun, c(-pi, pi - 1e-8))$root)

  #alternative: Interpolation with splines. When the qvonmisesmix() is called
  #for the first time after loading the package, the cdf values are calculated and
  #stored in a  new environment. This takes some time. All subsequent calls to qvonmisesmix()
  #use the stored cdf values and are very fast.

  npts <- 1000
  pts <- seq(-pi, (pi - 0.000001), length.out = npts)

  interpolate_from_prev_calc <- function(mu, kappa, prop) {
    if (any(c(mu, kappa, prop) != cylcop_vonmisesmix.env$parameter_vals)) {
      stop("not the parameters in the environment")
    }
    result <-
      stats::spline(cylcop_vonmisesmix.env$cdf,
                    pts,
                    method = "hyman",
                    xout = p)$y
    return(result)
  }

  result <- tryCatch(
    interpolate_from_prev_calc(mu, kappa, prop),
    error = function(e) {
      #cylcop.env$cdf doesn't exist yet, i.e. first time calling the function
      #Or it corresponds to other parameter values than the ones the function is currently called with
      assign("parameter_vals", c(mu, kappa, prop), envir =
               cylcop_vonmisesmix.env)
      assign("cdf",
             pvonmisesmix(pts, mu, kappa, prop),
             envir = cylcop_vonmisesmix.env)
      result <-
        stats::spline(cylcop_vonmisesmix.env$cdf,
                      pts,
                      method = "hyman",
                      xout = p)$y
      return(result)
    }
  )
  result[which(p > 0.999999)] <- pi
  return(result)
}

cylcop_vonmisesmix.env <- new.env(parent = emptyenv())

#' Mixed von Mises Maximum Likelihood Estimates
#'
#' Computes the maximum likelihood estimates for the parameters of a mixed
#' von Mises distribution: the mean directions, the concentration parameters,
#' and the proportions of the distributions. The code is a simplified version of
#'  \code{movMF::\link[movMF]{movMF}()} with the added
#' feature of optionally fixed mean directions \insertCite{Hornik2014}{cylcop}.
#' @param theta \link[base]{numeric} \link[base]{vector} of angles.
#' @param mu (optional) \link[base]{numeric} \link[base]{vector} of length \code{ncomp}
#'  holding the mean directions (angles). If not specified
#' the mean directions are estimated.
#' @param ncomp positive \link[base]{integer} specifying the number of components
#'  of the mixture model.
#'
#' @details The function complements the '\pkg{circular}' package, which
#' provides functions to make maximum likelihood estimates of e.g. von Mises
#' (\code{circular::\link[circular]{mle.vonmises}()}), or wrapped Cauchy distributions
#' (\code{circular::\link[circular]{mle.wrappedcauchy}()})
#'
#' @return A list containing the optimized parameters \code{mu}, \code{kappa},
#' and \code{prop}.
#'
#' @examples set.seed(123)
#'
#' n <- 10
#' angles <- rvonmisesmix(n,
#'   mu = c(0, pi),
#'   kappa = c(2, 1),
#'   prop = c(0.4,0.6)
#' )
#' mle.vonmisesmix(theta = angles)
#' mle.vonmisesmix(theta = angles, mu = c(0, pi))
#' mle.vonmisesmix(theta = angles, ncomp=3)
#'
#' @references \insertRef{Hornik2014}{cylcop}.
#'
#' @seealso
#' \code{movMF::\link[movMF]{movMF}()},
#' \code{circular::\link[circular]{mle.vonmises}()},
#' \code{\link{dvonmisesmix}()},
#' \code{\link{qvonmisesmix}()}.
#'
#'
#' @export
#'
mle.vonmisesmix <-  function(theta, mu = NULL, ncomp = 2)  {
  #validate input
  tryCatch({
    check_arg_all(check_argument_type(theta,
                                      type="numeric")
                  ,1)
    check_arg_all(list(check_argument_type(mu,
                                      type="numeric"),
                       check_argument_type(mu,
                                           type="NULL"))
                  ,2)
    check_arg_all(check_argument_type(ncomp,
                                      type="numeric",
                                      lower=1,
                                      integer=T)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  if (!is.null(mu) && length(mu) != ncomp) {
    stop(
      paste0(
        "Need to either estimate the means mu or provide a vector holding ",
        ncomp,
        " numerical values"
      )
    )
  }

  #convert angles to points on unit circle
  x <- cbind(cos(half2full_circ(theta)), sin(half2full_circ(theta)))
  if (!is.null(mu)) {
    mu_circ <- cbind(cos(half2full_circ(mu)), sin(half2full_circ(mu)))
  }
  #number of components in mixture
  k <- ncomp
  n <- nrow(x)
  maxiter <- 600

  do_kappa <- function(norms, alpha) {
    solve_kappa_Newton_Fourier(Rbar = (norms / (n * alpha)), d = 2)
  }

  G <- NULL

  M <- as.matrix(x[sample.int(nrow(x), k), , drop = FALSE])
  D <- pmax(1 - tcrossprod(x, M), 0)
  P <- memberships_from_cross_dissimilarities(D, 2)


  opt_old <- list()
  L_old <- -Inf
  iter <- 0L

  while (iter < maxiter) {
    ## M step.
    #max of expected likelihood of the mixture class probabilities
    alpha <- colMeans(P)
    if (any(alpha <= 0))
      stop("Cannot handle empty components")
    #max of expected likelihood of mean mu or given mu
    M <- crossprod(P, x)
    norms <- sqrt(rowSums(M ^ 2, dims = 1))
    if (is.null(mu)) {
      M <- M / ifelse(norms > 0, norms, 1)
    }
    else{
      M <- mu_circ
    }
    #max of expected likelihood of concentration kappa
    kappa <- do_kappa(norms, alpha)

    ## E step
    #posterior probability
    G <- cadd(tcrossprod(x, kappa * M),
              log(alpha) -  lH(kappa, 0))
    g <- log_row_sums(G)
    L_new <- sum(g)
    if (L_new > L_old) {
      L_old <- L_new
      #collect parameters, output
      param_theta <- kappa * M
      for (i in 1:ncomp) {
        opt_old$mu[i] <- atan2(y = param_theta[i, 2], x = param_theta[i, 1])
        opt_old$kappa[i] <- sqrt(sum(param_theta[i, ] ^ 2))
      }
      opt_old$prop <- alpha
    }
    P <- exp(G  - g)
    iter <- iter + 1L
  }

  if (is.null(opt_old))
    stop("EM algorithm did not converge for any run")
  return(opt_old)
}





#-------utility functions copied without modification from package movMF version
#-------0.2-3, movMF.R----------------
cadd <- function(A, x) {
  A + rep.int(x, rep.int(nrow(A), ncol(A)))
}

log_row_sums <- function(x) {
  M <- x[cbind(seq_len(nrow(x)), max.col(x, "first"))]
  M + log(rowSums(exp(x - M)))
}


solve_kappa_Newton_Fourier <-
  function(Rbar,
           d,
           tol = 1e-6,
           maxiter = 100L) {
    if (any(Rbar >= 1))
      stop("Cannot handle infinite concentration parameters")
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)

    sapply(seq_along(Rbar),
           function(i) {
             r <- Rbar[i]
             D <- d[i]
             nu <- D / 2 - 1
             lower <- Rinv_lower_Amos_bound(r, nu)
             upper <- Rinv_upper_Amos_bound(r, nu)
             iter <- 1L
             while (iter <= maxiter) {
               AA <- A(lower, D)
               AAprime <- Aprime(lower, D, A = AA)
               lower <- lower - (AA - r) / AAprime
               AA <- A(upper, D)
               upper <- upper - (AA - r) / AAprime
               if ((upper - lower) < tol * (lower + upper)) {
                 if ((upper - lower) < -tol * (lower + upper))
                   stop("no convergence")
                 break
               }
               iter <- iter + 1L
             }
             ## <FIXME>
             ## What should we really return?
             (lower + upper) / 2
             ## </FIXME>
           })
  }

beta_SS <- function(nu) {
  sqrt((nu + 1 / 2) * (nu + 3 / 2))
}

G <- function(kappa, alpha, beta) {
  n <- max(length(kappa), length(alpha), length(beta))
  kappa <- rep_len(kappa, n)
  alpha <- rep_len(alpha, n)
  beta <- rep_len(beta, n)

  y <- numeric(n)

  ind <- ((alpha == 0) & (beta == 0))
  if (any(ind))
    y[ind] <- 1
  ind <- !ind
  if (any(ind)) {
    denom <- alpha[ind] + sqrt(kappa[ind] ^ 2 + beta[ind] ^ 2)
    y[ind] <- ifelse(denom == 0, NaN, kappa[ind] / denom)
  }

  y
}

Ginv <- function(rho, alpha, beta) {
  ## Perhaps recycle arguments eventually ...

  sigma <- rho ^ 2
  rho * (alpha + sqrt(alpha ^ 2 * sigma + beta ^ 2 * (1 - sigma))) /
    (1 - sigma)
}

Rinv_lower_Amos_bound <- function(rho, nu) {
  ## Perhaps recycle arguments eventually ...
  pmax(Ginv(rho, nu, nu + 2),
       Ginv(rho, nu + 1 / 2, beta_SS(nu)))
}

Rinv_upper_Amos_bound <- function(rho, nu) {
  ## Perhaps recycle arguments eventually ...
  Ginv(rho, nu + 1 / 2, nu + 3 / 2)
}

S <- function(kappa, alpha, beta) {
  n <- max(length(kappa), length(alpha), length(beta))
  kappa <- rep_len(abs(kappa), n)
  alpha <- rep_len(alpha, n)
  beta <- rep_len(abs(beta), n)

  ab <- alpha + beta
  s <- double(n)
  ind <- (ab < 0) & (kappa ^ 2 > alpha ^ 2 - beta ^ 2)
  if (any(ind))
    s[ind] <- NaN
  ind <- !ind
  if (any(ind)) {
    u <- sqrt(kappa[ind] ^ 2 + beta[ind] ^ 2)
    s[ind] <-
      u - beta[ind] - ifelse(alpha[ind] == 0, 0, alpha[ind] * log((alpha[ind] + u) / ab[ind]))
  }
  s
}


lH_asymptotic <- function(kappa, nu) {
  ## Compute a suitable asymptotic approximation to
  ## \log(H_\nu(\kappa)).

  ## Add range checking eventually.
  n <- max(length(kappa), length(nu))
  kappa <- rep_len(kappa, n)
  nu <- rep_len(nu, n)

  y <- double(n)
  ipk <- (kappa > 0)
  ipn <- (nu > 0)
  ind <- ipk & !ipn
  if (any(ind)) {
    ## For \log(H_0) = \log(I_0), use the asymptotic approximation
    ##   I_0(\kappa) \approx e^\kappa / \sqrt{2 \pi \kappa}
    ## (e.g., http://dlmf.nist.gov/10.40).
    y[ind] <- kappa[ind] - log(2 * pi * kappa[ind]) / 2
  }
  ind <- ipk & ipn
  if (any(ind)) {
    ## For \nu > 0, use the Amos-type based approximation discussed
    ## above.
    kappa <- kappa[ind]
    nu <- nu[ind]
    beta <- beta_SS(nu)
    kappa_L <- pmin(kappa, sqrt((3 * nu + 11 / 2) * (nu + 3 / 2)))
    y[ind] <- S(kappa, nu + 1 / 2, beta) +
      (S(kappa_L, nu, nu + 2) - S(kappa_L, nu + 1 / 2, beta))
  }
  y
}

lH <- function(kappa, nu) {
  ## Compute \log(H_\nu(\kappa)) (or an approximation to it).
  ## See the implementation notes for details.

  n <- max(length(kappa), length(nu))
  kappa <- rep_len(kappa, n)
  nu <- rep_len(nu, n)

  y <- lH_asymptotic(kappa, nu)
  ## If the value from the asymptotic approximation is small enough,
  ## we can use direct Taylor series summation.
  ind <- y <= 699.5
  if (any(ind))
    y[ind] <- log(H(kappa[ind], nu[ind]))
  ## For intermediate values of the asymptotic approximation, we use
  ## rescaling and direct Taylor series summation.
  ind <- !ind & (y <= 1399)
  if (any(ind)) {
    v <- y[ind] / 2
    y[ind] <- v + log(H(kappa[ind], nu[ind], exp(-v)))
  }
  ## (Otherwise, we use the asymptotic approximation.)
  y
}

#-------utility functions adapted, but modified from package movMF version 0.2-3,
#--------movMF.R----------------

## Utility functions for computing H and log(H).

H <- function(kappa, nu, v0 = 1) {
  ## Compute v0 H_\nu(\kappa) by direct Taylor series summation.

  ## Add range checking eventually.
  n <- max(length(kappa), length(nu), length(v0))
  kappa <- rep_len(kappa, n)
  nu <- rep_len(nu, n)
  v0 <- rep_len(v0, n)

  my0F1(as.integer(n),
        as.double(kappa ^ 2 / 4),
        as.double(nu + 1),
        as.double(v0),
        y = double(n))
}

## Translated funtion my0F1 in movMF/src/mysum.c from C to R
my0F1 <- function(n, x, nu, y0, y) {
  for (i in 1:n) {
    z <- x[i]
    t <- y0[i]
    u <- nu[i]
    v <- 1.0
    r <- 0.0
    s <- t

    #Iterate until terms reach maximum.
    m <- (-u - 1 + sqrt((u - 1) * (u - 1)  + 4 * z)) / 2 + 1
    while (v < m) {
      r <- s
      t <- t * (z / (u * v))
      s <- s + t
      u <- u + 1
      v <- v + 1
    }
    #Iterate until numeric convergence if possible.
    while (s > r) {
      r <- s

      t <- t * (z / (u * v))
      s <- s + t

      u <- u + 1
      v <- v + 1
    }
    y[i] <- s
  }
  return(y)
}

A <- function(kappa, d) {
  n <- max(length(kappa), length(d))
  kappa <- rep_len(kappa, n)
  d <- rep_len(d, n)
  A <- vector("numeric", length = n)

  tol <- 1e-6
  index <- kappa >= tol
  if (sum(index)) {
    ## Compute A_d via a ratio of H functions:
    ##   A_d(\kappa)
    ##     = (kappa/d) * H_{d/2}(\kappa) / H_{d/2-1}(\kappa).
    s <- d[index] / 2 - 1
    kappai <- kappa[index]
    y <- kappai / d[index]
    a <- lH_asymptotic(kappai, s + 1)
    ind <- (a <= 699.5)
    if (any(ind))
      y[ind] <- y[ind] *
      (H(kappai[ind], s + 1) /
         H(kappai[ind], s))
    ind <- !ind & (a <= 1399)
    if (any(ind)) {
      v <- exp(-a[ind] / 2)
      y[ind] <- y[ind] *
        (H(kappai[ind], s + 1, v) /
           H(kappai[ind], s, v))
    }
    ind <- (a > 1399)
    if (any(ind))
      y[ind] <- y[ind] *
      exp(a[ind] - lH_asymptotic(kappai[ind], s))
    if (any(y >= 1))
      stop("RH evaluation gave infeasible values which are not in the range [0, 1)")
    A[index] <- y
  }
  if (sum(!index)) {
    di <- d[!index]
    kappai <- kappa[!index]
    A[!index] <-
      kappai / di - kappai ^ 3 / (di ^ 2 * (di + 2)) + 2 * kappai ^ 5 / (di ^
                                                                           3 * (di + 2) * (di + 4))
  }
  A
}

Aprime <-
  function(kappa, d, A = NULL)
  {
    n <- max(length(kappa), length(d), length(A))
    kappa <- rep_len(kappa, n)
    d <- rep_len(d, n)
    if (!is.null(A))
      A <- rep_len(A, n)
    aprime <- vector("numeric", length = n)
    tol <- 1e-6
    index <- kappa >= tol
    if (sum(index)) {
      if (is.null(A))
        A <- A(kappa[index], d[index])
      else
        A <- A[index]
      aprime[index] <- 1 - A ^ 2 - A * (d[index] - 1) / kappa[index]
    }
    if (sum(!index)) {
      di <- d[!index]
      aprime[!index] <-
        1 / di - 3 / (di ^ 2 * (di + 2)) * kappa[!index] ^ 2 + 10 / (di ^ 3 * (di + 2) * (di + 4)) * kappa[!index] ^
        4
    }
    aprime
  }



#-------utility functions copied without modification from package clue version 0.3-42, membership.R----------------
memberships_from_cross_dissimilarities <- function(d, power = 2) {
  ## For a given matrix of cross-dissimilarities [d_{bj}], return a
  ## matrix [u_{bj}] such that \sum_{b,j} u_{bj}^p d_{bj}^q => min!
  ## under the constraint that u is a stochastic matrix.
  ## If only one power is given, it is taken as p, with q as 1.
  ## <NOTE>
  ## This returns a plain matrix of membership values and not a
  ## cl_membership object (so that it does not deal with possibly
  ## dropping or re-introducing unused classes).
  ## </NOTE>
  exponent <- if (length(power) == 1L)
    1 / (1 - power)
  else
    power[2L] / (1 - power[1L])
  u <- matrix(0, nrow(d), ncol(d))
  zero_incidences <- !(d > 0)
  n_of_zeroes <- rowSums(zero_incidences)
  if (any(ind <- (n_of_zeroes > 0)))
    u[ind,] <-
    zero_incidences[ind, , drop = FALSE] / n_of_zeroes[ind]
  if (any(!ind)) {
    ## Compute d_{bj}^e / \sum_k d_{bk}^e without overflow from very
    ## small d_{bj} values.
    d <- exponent * log(d[!ind, , drop = FALSE])
    d <- exp(d - d[cbind(seq_len(nrow(d)), max.col(d))])
    u[!ind,] <- d / rowSums(d)
  }
  u
}
