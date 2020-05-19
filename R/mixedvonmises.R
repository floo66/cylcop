#' Quantiles of the mixed vonMises distribution
#'
#' @param p vector of probabilities
#' @param mu1 mean of the first component
#' @param mu2 mean of the second component
#' @param kappa1 concentration parameter of the first component
#' @param kappa2 concentration parameter of the first component
#' @param prop mixing proportion
#'
#' @return quantiles of the mixed vonMises distribution.
#' @export
#'
qmixedvonmises <- function(p, mu1, mu2, kappa1, kappa2, prop) {

  if(p>0.999999) return(pi)

  #Numerically find the inverse. This can get pretty slow, because the entire calculation is repeated
  #every time qmixedvonmises() is called
  # fun <- function(q) pmixedvonmises(q, mu1, mu2, kappa1, kappa2, prop, from=-pi) - p
  # result <- suppressWarnings(uniroot(fun, c(-pi, pi - 1e-8))$root)

  #alternative: Interpolation with splines. When the qmixedvonmises() is called
  #for the first time after loading the package, the cdf values are calculated and
  #stored in a  new environment. This takes some time. All subsequent calls to qmixedvonmises()
  #use the stored cdf values and are very fast.

  npts<-1000
  pts<-seq(-pi,(pi-0.000001),length.out = npts)

  result <- tryCatch(stats::spline(cylcop_mixedvonmises.env$cdf, pts, method = "hyman", xout = p)$y,
            error=function(e){
              #cylcop.env$cdf doesn't exist yet, i.e. first time calling the function
              assign("cdf",
                     suppressWarnings(pmixedvonmises(pts, mu1, mu2, kappa1, kappa2, prop, from=-pi)),
                     envir=cylcop_mixedvonmises.env)
              result <- stats::spline(cylcop_mixedvonmises.env$cdf, pts, method = "hyman", xout = p)$y
              return(result)
            })
  return(result)
}

cylcop_mixedvonmises.env <- new.env(parent = emptyenv())

#' mixed von Mises Maximum Likelihood Estimates
#'
#'Computes the maximum likelihood estimates for the parameters of a von Mises distribution:
#'the mean directions, the concentration parameters, and the proportion of the 2 distributions.
#'
#' @details Code is a simplified version of \code{movMF::movMF()} with the added
#' feature of optionally fixed mean directions, see references.
#' @param theta vector of angles
#' @param mu (optional) vector holding the 2 mean directions (angles). If not specified
#' the mean directions are estimated.
#'
#' @return A list containing the optimized parameters \code{mu}, \code{kappa}, and \code{prop}.
#' @references Hornik, K. and Grun, B. (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions. _J. Stat. Softw., 58(10)_, 164-168.
#'
#' @export
#'
mle.mixedvonmises <-  function(theta, mu=NULL)  {
  if(!is.null(mu) && length(mu)!=2){
    stop("Need to either estimate the means mu or provide a vector holding 2 numerical values")
  }

  #convert angles to points on unit circle
  x <- cbind(cos(cylcop::half2full_circ(theta)), sin(cylcop::half2full_circ(theta)))
  if(!is.null(mu)){
    mu_circ<-matrix(nrow=2,ncol=2)
    mu_circ[1,] <- cbind(cos(cylcop::half2full_circ(mu[1])), sin(cylcop::half2full_circ(mu[1])))
    mu_circ[2,] <- cbind(cos(cylcop::half2full_circ(mu[2])), sin(cylcop::half2full_circ(mu[2])))
  }
  #number of components in mixture
  k<-2
  n <- nrow(x)
  maxiter <- 600

  do_kappa <- function(norms, alpha){
    solve_kappa_Newton_Fourier(Rbar=(norms / (n * alpha)), d=2)}

  G <- NULL

  M <- as.matrix(x[sample.int(nrow(x), k), , drop = FALSE])
  D <- pmax(1 - tcrossprod(x, M), 0)
  P <- memberships_from_cross_dissimilarities(D, 2)


  opt_old <- list()
  L_old <- -Inf
  iter <- 0L

  while(iter < maxiter) {

    ## M step.
    #max of expected likelihood of the mixture class probabilities
    alpha <- colMeans(P)
    if(any(alpha <= 0))
      stop("Cannot handle empty components")
    #max of expected likelihood of mean mu or given mu
    M <- crossprod(P, x)
    norms <- sqrt(rowSums(M^2, dims = 1))
    if(is.null(mu)){
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
    if(L_new > L_old) {
      L_old <- L_new
      #collect parameters, output
      param_theta <- kappa * M
      opt_old$mu[1] <- atan2(y = param_theta[1,2],x=param_theta[1,1])
      opt_old$mu[2] <- atan2(y = param_theta[2,2],x=param_theta[2,1])
      opt_old$kappa[1] <- sqrt(sum(param_theta[1,]^2))
      opt_old$kappa[2] <- sqrt(sum(param_theta[2,]^2))
      opt_old$prop <- alpha[1]
    }
    P <- exp(G  - g)
    iter <- iter + 1L
  }

  if(is.null(opt_old))
    stop("EM algorithm did not converge for any run")
  return(opt_old)
}





#-------utility functions copied without modification from package movMF version 0.2-3, movMF.R----------------
cadd <- function(A, x){
  A + rep.int(x, rep.int(nrow(A), ncol(A)))
}

log_row_sums <- function(x){
  M <- x[cbind(seq_len(nrow(x)), max.col(x, "first"))]
  M + log(rowSums(exp(x - M)))
}


solve_kappa_Newton_Fourier <- function(Rbar, d, tol = 1e-6, maxiter = 100L){
  if(any(Rbar >= 1))
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
           while(iter <= maxiter) {
             AA <- A(lower, D)
             AAprime <- Aprime(lower, D, A = AA)
             lower <- lower - (AA - r) / AAprime
             AA <- A(upper, D)
             upper <- upper - (AA - r) / AAprime
             if((upper - lower) < tol * (lower + upper)) {
               if ((upper - lower) < - tol * (lower + upper))
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

beta_SS <- function(nu){
  sqrt((nu + 1/2) * (nu + 3/2))
}

G <- function(kappa, alpha, beta){
  n <- max(length(kappa), length(alpha), length(beta))
  kappa <- rep_len(kappa, n)
  alpha <- rep_len(alpha, n)
  beta <- rep_len(beta, n)

  y <- numeric(n)

  ind <- ((alpha == 0) & (beta == 0))
  if(any(ind)) y[ind] <- 1
  ind <- !ind
  if(any(ind)) {
    denom <- alpha[ind] + sqrt(kappa[ind]^2 + beta[ind]^2)
    y[ind] <- ifelse(denom == 0, NaN, kappa[ind] / denom)
  }

  y
}

Ginv <- function(rho, alpha, beta){
  ## Perhaps recycle arguments eventually ...

  sigma <- rho^2
  rho * (alpha + sqrt(alpha^2 * sigma + beta^2 * (1 - sigma))) /
    (1 - sigma)
}

Rinv_lower_Amos_bound <- function(rho, nu){
  ## Perhaps recycle arguments eventually ...
  pmax(Ginv(rho, nu, nu + 2),
       Ginv(rho, nu + 1/2, beta_SS(nu)))
}

Rinv_upper_Amos_bound <- function(rho, nu){
  ## Perhaps recycle arguments eventually ...
  Ginv(rho, nu + 1/2, nu + 3/2)
}

S <- function(kappa, alpha, beta){
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


lH_asymptotic <- function(kappa, nu){
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
  if(any(ind)) {
    ## For \log(H_0) = \log(I_0), use the asymptotic approximation
    ##   I_0(\kappa) \approx e^\kappa / \sqrt{2 \pi \kappa}
    ## (e.g., http://dlmf.nist.gov/10.40).
    y[ind] <- kappa[ind] - log(2 * pi * kappa[ind]) / 2
  }
  ind <- ipk & ipn
  if(any(ind)) {
    ## For \nu > 0, use the Amos-type based approximation discussed
    ## above.
    kappa <- kappa[ind]
    nu <- nu[ind]
    beta <- beta_SS(nu)
    kappa_L <- pmin(kappa, sqrt((3 * nu + 11 / 2) * (nu + 3 / 2)))
    y[ind] <- S(kappa, nu + 1/2, beta) +
      (S(kappa_L, nu, nu + 2) - S(kappa_L, nu + 1/2, beta))
  }
  y
}

lH <-function(kappa, nu){
  ## Compute \log(H_\nu(\kappa)) (or an approximation to it).
  ## See the implementation notes for details.

  n <- max(length(kappa), length(nu))
  kappa <- rep_len(kappa, n)
  nu <- rep_len(nu, n)

  y <- lH_asymptotic(kappa, nu)
  ## If the value from the asymptotic approximation is small enough,
  ## we can use direct Taylor series summation.
  ind <- y <= 699.5
  if(any(ind))
    y[ind] <- log(H(kappa[ind], nu[ind]))
  ## For intermediate values of the asymptotic approximation, we use
  ## rescaling and direct Taylor series summation.
  ind <- !ind & (y <= 1399)
  if(any(ind)) {
    v <- y[ind] / 2
    y[ind] <- v + log(H(kappa[ind], nu[ind], exp(-v)))
  }
  ## (Otherwise, we use the asymptotic approximation.)
  y
}

#-------utility functions adapted, but modified from package movMF version 0.2-3, movMF.R----------------

## Utility functions for computing H and log(H).

H <- function(kappa, nu, v0 = 1){
  ## Compute v0 H_\nu(\kappa) by direct Taylor series summation.

  ## Add range checking eventually.
  n <- max(length(kappa), length(nu), length(v0))
  kappa <- rep_len(kappa, n)
  nu <- rep_len(nu, n)
  v0 <- rep_len(v0, n)

  my0F1(
    as.integer(n),
    as.double(kappa ^ 2 / 4),
    as.double(nu + 1),
    as.double(v0),
    y = double(n)
  )
}

## Translated funtion my0F1 in movMF/src/mysum.c from C to R
my0F1 <- function(n, x, nu, y0, y){
  for (i in 1:n) {
    z <- x[i]
    t <- y0[i]
    u <- nu[i]
    v <- 1.0
    r <- 0.0
    s <- t

    #Iterate until terms reach maximum.
    m <- (- u - 1 + sqrt((u - 1) * (u - 1)  + 4 * z)) / 2 + 1
    while(v < m) {
      r <- s
      t <- t*(z / (u * v))
      s <- s+t
      u<-u+1
      v<-v+1
    }
    #Iterate until numeric convergence if possible.
    while(s > r) {
      r <- s;
      t <- t*(z / (u * v))
      s <- s+t;
      u<-u+1
      v<-v+1
    }
    y[i] <- s
  }
  return(y)
}

A <- function(kappa, d){
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
    if(any(ind))
      y[ind] <- y[ind] *
      (H(kappai[ind], s + 1) /
         H(kappai[ind], s))
    ind <- !ind & (a <= 1399)
    if(any(ind)) {
      v <- exp(- a[ind] / 2)
      y[ind] <- y[ind] *
        (H(kappai[ind], s + 1, v) /
           H(kappai[ind], s, v))
    }
    ind <- (a > 1399)
    if(any(ind))
      y[ind] <- y[ind] *
      exp(a[ind] - lH_asymptotic(kappai[ind], s))
    if(any(y >= 1))
      stop("RH evaluation gave infeasible values which are not in the range [0, 1)")
    A[index] <- y
  }
  if (sum(!index)) {
    di <- d[!index]
    kappai <- kappa[!index]
    A[!index] <- kappai / di - kappai^3 / (di^2 * (di + 2)) + 2 * kappai^5 / (di^3 * (di + 2) * (di + 4))
  }
  A
}

Aprime <-
  function(kappa, d, A = NULL)
  {
    n <- max(length(kappa), length(d), length(A))
    kappa <- rep_len(kappa, n)
    d <- rep_len(d, n)
    if(!is.null(A))
      A <- rep_len(A, n)
    aprime <- vector("numeric", length = n)
    tol <- 1e-6
    index <- kappa >= tol
    if (sum(index)) {
      if(is.null(A))
        A <- A(kappa[index], d[index])
      else
        A <- A[index]
      aprime[index] <- 1 - A ^ 2 - A * (d[index] - 1) / kappa[index]
    }
    if (sum(!index)) {
      di <- d[!index]
      aprime[!index] <- 1 / di - 3 / (di^2 * (di + 2)) * kappa[!index]^2 + 10 / (di^3 * (di + 2) * (di + 4)) * kappa[!index]^4
    }
    aprime
  }



#-------utility functions copied without modification from package clue version 0.3-42, membership.R----------------
memberships_from_cross_dissimilarities <- function(d, power = 2){
  ## For a given matrix of cross-dissimilarities [d_{bj}], return a
  ## matrix [u_{bj}] such that \sum_{b,j} u_{bj}^p d_{bj}^q => min!
  ## under the constraint that u is a stochastic matrix.
  ## If only one power is given, it is taken as p, with q as 1.
  ## <NOTE>
  ## This returns a plain matrix of membership values and not a
  ## cl_membership object (so that it does not deal with possibly
  ## dropping or re-introducing unused classes).
  ## </NOTE>
  exponent <- if(length(power) == 1L)
    1 / (1 - power)
  else
    power[2L] / (1 - power[1L])
  u <- matrix(0, nrow(d), ncol(d))
  zero_incidences <- !(d > 0)
  n_of_zeroes <- rowSums(zero_incidences)
  if(any(ind <- (n_of_zeroes > 0)))
    u[ind, ] <-
    zero_incidences[ind, , drop = FALSE] / n_of_zeroes[ind]
  if(any(!ind)) {
    ## Compute d_{bj}^e / \sum_k d_{bk}^e without overflow from very
    ## small d_{bj} values.
    d <- exponent * log(d[!ind, , drop = FALSE])
    d <- exp(d - d[cbind(seq_len(nrow(d)), max.col(d))])
    u[!ind, ] <- d / rowSums(d)
  }
  u
}
