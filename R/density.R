#' Density, distribution, random number generation and quantiles for a kernel density estimate
#'
#' @param density A list containing information about the kernel density
#'   estimate. The structure of the list must be as returned by \code{fit_angle()} or
#'   \code{fit_steplength()}.
#' @param n Number of random samples to be generated.
#' @name dens
#' @export
#'
rdens <- function(n, density) {
  # Draw a points from the domain of the KDE
  sample(as.double(density$x),
         n,
         replace = TRUE,
         prob = density$y)

  #Now draw from the kernel distributions at these points,
  #with kappa (for vonmises) or sd (for Gaussian) equal to the bandwidth
  kernel<-get_marg(density$kernel)
  map(sample, ~kernel$r(1,sample,dens$bw))
}



# Probability density from a kernel density estimate
#
#' @rdname dens
#' @param x A numeric vector of measurements.
#'
#' @export
#'
ddens <- function(x, density) {
  dens_interpol<-stats::approxfun(density$x,density$y)
  map_dbl(x,~dens_interpol(.x))
}



# Distribution function from a kernel density estimate
#
#' @rdname dens
#'
#' @export
#'
pdens <- function(x, density) {

  # could do something like this (below) to precisely calculate the cdf, but this is quite slow,
  # better use linear interpolation
  # rowMeans(pnorm(outer(x,na.omit(traj$steplength),"-"),0,density$bw))

  density$x <- as.double(density$x)
  density$y <- as.double(density$y)
  delx <-
    (density$x[length(density$x)] - density$x[1]) / length(density$x)

  #possible speedup by saving preious values of sum and just adding a term, already fast enough, though
  #right Riemann sum, then set the first value, cdf(min(x))=0
  cdf<-map_dbl(seq(1,length(density$y)-1,1),~sum(density$y[1:.x]) * delx)
  cdf <- c(0,cdf,1)

  #add one more x-value, its cdf is 1
  density$x<-c(density$x,2*density$x[length(density$x)]-density$x[length(density$x)-1])
  cdf_interpol <- stats::approxfun(c(density$x,Inf),c(cdf,1))
  map_dbl(x,~cdf_interpol(.x))
}



# Quantiles of a kernel density estimate
#
#' @rdname dens
#' @param p A numeric vector of probabilities.
#'
#' @export
#'
qdens <- function(p, density) {
  density$x <- as.double(density$x)
  density$y <- as.double(density$y)
  delx <-
    (density$x[length(density$x)] - density$x[1]) / length(density$x)

  bin_search <- function(p) {
    l <- 0L
    r <- length(density$x)
    while (abs(r - l) > 1) {
      m = floor((l + r) / 2)
      prob <- sum(density$y[1:m]) * delx
      if (prob < p)
        l <- m
      else if (prob > p)
        r <- m
    }

    cdf<-function(index){
      if(index==0L) return(0)
      else if(index==length(density$y)) return(1)
      else return(sum(density$y[1:index]) * delx)
      }

    #add one more x-value, its cdf is 1
    density$x<-c(density$x,2*density$x[length(density$x)]-density$x[length(density$x)-1])

    #linear interpolation, faster than precise calcualtion, see pdens()
    q_interpol <- stats::approxfun(c(cdf(l), cdf(r)),c(density$x[l+1],density$x[r+1]))
    return(q_interpol(p))
  }
  map_dbl(p,  ~ bin_search(.x))
}
