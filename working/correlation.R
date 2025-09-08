#' Estimate a Rank-Based Circular-Linear Correlation Coefficient
#'
#' The code is based on \insertCite{Mardia1976;textual}{cylcop},
#' \insertCite{Solow1988;textual}{cylcop} and \insertCite{Tu2015;textual}{cylcop}.
#' The function returns a numeric value between 0 and 1, not -1 and 1, positive
#' and negative correlation cannot be discerned. Note also that the correlation
#' coefficient is independent of the marginal distributions.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable).
#'
#' @return A \link[base]{numeric} value between 0 and 1, the circular-linear
#' correlation coefficient.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#'
#' #draw samples and calculate the correlation coefficient
#' sample <- rcylcop(100, cop)
#' cor_cyl(theta = sample[,1], x = sample[,2])
#'
#' #the correlation coefficient is independent of the marginal distribution.
#' sample <- traj_sim(100,
#'   cop,
#'   marginal_circ = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_lin = list(name = "weibull", coef = list(shape = 2))
#' )
#' cor_cyl(theta = sample$angle, x = sample$steplength)
#' cor_cyl(theta = sample$cop_u, x = sample$cop_v)
#'
#' # Estimate correlation of samples drawn from circular-linear copulas with
#' # perfect correlation
#' cop <- cyl_rect_combine(copula::normalCopula(1))
#' sample <- rcylcop(100, cop)
#' cor_cyl(theta = sample[,1], x = sample[,2])
#'
#' @references \insertRef{Mardia1976}{cylcop}
#'
#' \insertRef{Solow1988}{cylcop}
#'
#' \insertRef{Tu2015}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{mi_cyl_bin}()}, \code{\link{fit_cylcop_cor}()}.
#'
#' @export
#'
cor_cyl <- function(theta, x) {
  tryCatch({
    check_arg_all(check_argument_type(theta, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x, type="numeric")
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

#Variable names as in Solow 1988 but with radians instead of degrees.

# Assigning ranks to angular and linear measurements

  data <- data.frame(theta, x) %>% na.omit() %>% dplyr::arrange(x) %>%
    mutate(r_theta = data.table::frank(.data$theta, ties.method = "average"))
  n <- nrow(data)
  data <- mutate(data, r_theta_star = .data$r_theta * 2 * pi / n) %>%
    mutate(r_x = data.table::frank(.data$x, ties.method = "average"))


# Calcualte correlation

  C <- dplyr::summarize(data, sum(.data$r_x * cos(.data$r_theta_star))) %>% as.double()
  S <- dplyr::summarize(data, sum(.data$r_x * sin
                           (.data$r_theta_star))) %>% as.double()

  if (n %% 2 == 0) {
    a <- 1 / (1 + 5 * (1 / tan(pi / n)) ^ 2 + 4 * (1 / tan(pi / n)) ^ 4)
  }
  else{
    a <- 2 * (sin(pi / n)) ^ 4 / ((1 + (cos(pi / n))) ^ 3)
  }

  D <- a * (C ^ 2 + S ^ 2)

  return(D)
}



#' Estimate the Mutual Information Between a Circular and a Linear Random
#' Variable with a Binning Approach
#'
#' The empirical copula is obtained from the data (\code{theta} and \code{x}),
#' and the mutual information of the 2 components is calculated by binning the data.
#' This gives a non-negative number that can be normalized to lie between 0 and 1.
#'
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles (measurements of a circular
#'   variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths (measurements of a linear
#'   variable).
#' @param normalize \link[base]{logical} value whether the mutual information should be
#'   normalized to lie within \eqn{[0,1]}.
#' @param symmetrize \link[base]{logical} value whether it should be assumed that right and left
#' turns are equivalent. If \code{theta} can take values in \eqn{[-\pi, \pi)},
#' this means that positive and negative angles are equivalent.
#' @param nbins \link[base]{integer} value, the number of bins to use in each dimension.
#' The default value (when \code{nbins} is \code{NULL}) is \code{length(theta)^(1/3)}).
#' @param method \link[base]{character} string, the name of the entropy estimator.
#' At the moment, \code{"ML"} and \code{"MM"} are implemented.
#'
#' @details First, the two components of the empirical copula, \eqn{u} and \eqn{v}
#' are obtained. Then the differential entropy of the joint distirbution is estimated via
#' discretizing \eqn{u} and \eqn{v} into \code{nbins} bins each.
#' The mutual information can then be calculated from the entropy. It cam also
#' normalized to lie between 0 and 1 by dividing by the entropy of the marginal distribution.
#' Since the marginal distribution is uniform, the entropy is just the log of the
#' number of bins.
#'
#' Even if \code{u} and \code{v} are perfectly correlated
#' (i.e. \code{\link{cor_cyl}} goes to 1 with large sample sizes),
#' the normalized mutual information will not be 1 if the underlying copula is periodic and
#' symmetric. E.g. while \code{normalCopula(1)} has a correlation of 1 and a density
#' that looks like a line going from \eqn{(0,0)} to \eqn{(1,1)},
#'  \code{cyl_rect_combine(normalCopula(1))}
#'  has a density that looks like "<". The mutual information will be 1 in the first case,
#'  but not in the second. Therefore, we can set \code{symmetrize = TRUE} to first
#'  convert (if necessary) theta to lie in \eqn{[-\pi, \pi)} and then multiply all angles
#'  larger than 0 with -1. The empirical copula is then calculated and the mutual information
#' is obtained from those values. It is exactly 1 in the case of
#' perfect correlation as captured by e.g.
#' \code{cyl_rect_combine(normalCopula(1))}.
#'
#' Note also that the mutual information is independent of the marginal distributions.
#' However, \code{symmetrize=TRUE} only works with angles, not with pseudo-observations.
#' When \code{x} and \code{theta} are pseudo-observations, information is lost
#' due to the ranking, and symmetrization will fail.
#'
#' Currently there are 2 methods implemented to calculate the joint entropy via binning.
#' with \code{method="ML"} the probability associated with each bin is just the maximum likelihood estimate,
#' i.e. the number of samples in that bin divided by the total number of samples.
#' \code{method="MM"} refers to the Miller–Madow estimator, which just adds a correction term to the
#' maximum likelihood estimate of the entropy to account for bias.
#'
#' @return A \link[base]{numeric} value, the mutual information between \code{theta} and \code{x}
#' in nats.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' marg1 <- list(name="vonmises",coef=list(0,4))
#' marg2 <- list(name="lnorm",coef=list(2,3))
#'
#' #draw samples and calculate the mutual information.
#' sample <- rjoint(100,cop,marg1,marg2)
#' mi_cyl_bin(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = TRUE,
#'   symmetrize = FALSE
#' )
#'
#' #the correlation coefficient is independent of the marginal distribution.
#'  sample <- traj_sim(100,
#'   cop,
#'   marginal_circ = list(name = "vonmises", coef  = list(0, 1)),
#'   marginal_lin = list(name = "weibull", coef = list(shape = 2))
#' )
#'
#' mi_cyl_bin(theta = sample$angle,
#'   x = sample$steplength,
#'   normalize = TRUE,
#'   symmetrize = FALSE)
#' mi_cyl_bin(theta = sample$cop_u,
#'   x = sample$cop_v,
#'   normalize = TRUE,
#'   symmetrize = FALSE)
#'
#' # Estimate correlation of samples drawn from circular-linear copulas
#' # with perfect correlation.
#' cop <- cyl_rect_combine(copula::normalCopula(1))
#' sample <- rjoint(100,cop,marg1,marg2)
#' # without normalization
#' mi_cyl_bin(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = FALSE,
#'   symmetrize = FALSE
#' )
#' #with normalization
#' mi_cyl_bin(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = TRUE,
#'   symmetrize = FALSE
#' )
#' #only with normalization and symmetrization do we get a value of 1
#' mi_cyl_bin(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = TRUE,
#'   symmetrize = TRUE
#' )
#'
#' @references \insertRef{Jian2011}{cylcop}
#'
#' \insertRef{Calsaverini2009}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cor_cyl}()}, \code{\link{fit_cylcop_cor}()}.
#'
#' @export
#'
mi_cyl_bin <- function(theta,
                      x,
                      normalize = TRUE,
                      symmetrize = FALSE,
                      nbins = NULL,
                      method="ML") {



  tryCatch({
    check_arg_all(check_argument_type(theta, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(normalize, type="logical")
                  ,1)
    check_arg_all(check_argument_type(symmetrize, type="logical")
                  ,1)
    check_arg_all(list(check_argument_type(nbins, type="numeric", length = 1, lower = 2, integer=T),
                  check_argument_type(nbins, type="NULL"))
                 ,2)
    check_arg_all(check_argument_type(method,
                                      type="character",
                                      values = c("ML", "MM"),
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


  if(is.null(nbins)){
    nbins <- ceiling(max(2,length(theta)^(1/3)))
  }

    data <- data.frame(theta, x) %>% na.omit()

  if (symmetrize) {
    data$theta <- full2half_circ(data$theta)
    ind <- which(data[,1]>0)
    data[ind,1] <- 0-data[ind,1]
  }

  #calculate empirical copula
  emp_cop <- pobs(data, ties.method = "average")%>%as.data.frame()



  #set up bins and count datapoint per bin
  bins <- seq(0, 1, length.out = nbins + 1)
  count <- as.data.frame(table(cut(emp_cop[, 1], breaks = bins), cut(emp_cop[, 2], breaks = bins)))
  p <-  count[,3] / nrow(emp_cop)
  ent <- p * log(p)
  empty_bins <- length(which(is.na(ent)))
  ent[is.na(ent)] <- 0
  ent <- -sum(ent)

  if(method=="MM"){
    ent <- ent+(nbins^2-1-empty_bins)/nrow(emp_cop)
  }

  #add correction because differential entropy is not limit of Shannon entropy
  mi <- -(ent + log(1 / (nbins^2)))

  if (normalize) {
    mi <- mi / log(nbins)
  }




  # discr_theta <- infotheo::discretize(emp_cop[,1],nbins=nbins)
  # discr_x <- infotheo::discretize(emp_cop[,2],nbins=nbins)
  # if (normalize) {
  #   mi <-
  #     infotheo::mutinformation(discr_theta,discr_x) / sqrt(infotheo::entropy(discr_theta) * infotheo::entropy(discr_x))
  # }
  #
  # else{
  #   mi <- infotheo::mutinformation(discr_theta,discr_x)
  # }

  return(mi)
}





#' Estimate the Mutual Information Between a Circular and a Linear Random
#' Variable with a K-Nearest-Neighbor Approach
#'
#' The empirical copula is obtained from the data (\code{theta} and \code{x}),
#' and the mutual information of the 2 components is calculated.
#'
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles (measurements of a circular
#'   variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths (measurements of a linear
#'   variable).
#' @param periodic \link[base]{logical} value whether the metric to calculate the k
#' nearest neighbors should be periodic (circular) in dimension theta. Note that
#' to set \code{symemtrize=T} and \code{periodic=T} will not make sense in most cases.
#' @param symmetrize \link[base]{logical} value whether it should be assumed that right and left
#' turns are equivalent. If \code{theta} can take values in \eqn{[-\pi, \pi)},
#' this means that positive and negative angles are equivalent.
#' @param k \link[base] {integer} value, the number of nearest neighbors to calculate.
#'
#' @details First, the two components of the empirical copula, \eqn{u} and \eqn{v}
#' are obtained. Then the joint differential entropy is estimated using a k-nearest
#' neighbor approach according to Kozachenko and Leonenko. The k-nearest neighbors
#' are found using the Euclidean norm. If \code{periodic=T}, the Euclidean distance is
#' calculated across a periodic boundary in the theta-dimension.
#'
#' For information and illustrations on the effect of \code{symmertrize}, see the
#' details examples of \code{\link{mi_cyl_bin}()}. Note, however, that the function
#' \code{mi_cyl_knn} does not allow for normalization of the mutual information.
#' (Except for a transformation to a non-linear correlation coefficient using the
#' function  \code{\link{mi2cc}()}).
#'
#' The number of nearest neighbors, \code{k}, influences the result values of 3 to 5
#' seem to work well.
#'
#' Note, however that the mutual information will be biased, with a bias of approximately
#' N^(0.5) where N is the number of samples. Due to this bias, it is possible that
#' for weak correlations, the mutual information estimate becomes positive. In those
#' cases it is probably best to assume the mutual information to be 0.
#'
#' @return A \link[base]{numeric} value, the mutual information between \code{theta} and \code{x}
#' in nats.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' marg1 <- list(name="vonmises",coef=list(0,4))
#' marg2 <- list(name="lnorm",coef=list(2,3))
#'
#' #draw samples and calculate the mutual information.
#' sample <- rjoint(100,cop,marg1,marg2)
#' mi_cyl_knn(
#'   theta = sample[,1],
#'   x = sample[,2],
#'   periodic = TRUE
#'   symmetrize = FALSE,
#'   k=3
#' )
#'
#' @references
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cor_cyl}()}, \code{\link{fit_cylcop_cor}()}.
#'
#' @export
#'
mi_cyl_knn <- function(theta,
                       x,
                       periodic = TRUE,
                       symmetrize = FALSE,
                       k = 3) {

#hardcoded buffer_width for periodic metric.
#Exact number not really important. All points falling outside the buffer width
# Will be considered in second sweep.Minor performance gain by tuning this to
# a specific problem. Maybe user input in the future.
buffer_width <- 0.1

#At the moment only Euclidean metric is implemented. Code for max-norm is commented out
#The problem is that ANN library has to be compiled fresh when the metric is switched
#So far, there are no packages on CRAN that give an interface to the ANN library
#a norm other than Euclidean. Some pacakges allow other norms, but only when you provide
#the pre-calculated distance matrix. Defeats the purpose, memory issues, performance issues.
#Maybe code something myself in teh future.
metric="euclid"


  tryCatch({
    check_arg_all(check_argument_type(theta, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(periodic, type="logical")
                  ,1)
    check_arg_all(check_argument_type(symmetrize, type="logical")
                  ,1)
    check_arg_all(check_argument_type(k, type="numeric", length = 1, lower = 1, integer=T)
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


  data <- data.frame(theta, x) %>% na.omit()

  if (symmetrize) {
    data$theta <- full2half_circ(data$theta)
    ind <- which(data[,1]>0)
    data[ind,1] <- 0-data[ind,1]
  }

  #calculate empirical copula
  emp_cop <- pobs(data, ties.method = "average")%>%as.data.frame()



  N <- nrow(emp_cop)
  ranks <- emp_cop*(N+1)
  ranks <- round(ranks[order(ranks[,1]),],0)
  if(periodic){
    buffer_width <- ceiling(buffer_width*N)
    ranks <- rbind(tail(ranks,buffer_width),ranks,head(ranks,buffer_width))
    ranks[,1] <- ranks[,1]+c(rep(-N,buffer_width),rep(0,N),rep(N,buffer_width))+buffer_width
    points_2_eval <- seq_len(N) + buffer_width
  }else{
    points_2_eval <- seq_len(N)
  }

  # if(metric=="max"){
  #   cd <- 1
  #   rotate_ranks_45 <- rotate(ranks,45)
  #   neighbors <-  RANN1::nn2(rotate_ranks_45,
  #                            rotate_ranks_45[points_2_eval,],
  #                            k=(k+1))
  #   dist <- round(neighbors$nn.dists[, (k+1)] / sqrt(1^2 + 1^2), 0)
  # }else if(metric=="euclid"){
    cd <- pi/4
    neighbors <-  dbscan::kNN(x=ranks,
                              query = ranks[points_2_eval,],
                              k=k+1,
                              bucketSize = 10,
                              sort=F,
                              search="kdtree",
                              splitRule="SL_MIDPT")
    dist <- neighbors$dist[, (k+1)]
  # }

  if(periodic){
    u_dist_nearest_border <- pmin(seq_len(N) + buffer_width - 1, N + buffer_width - seq_len(N))
    outside_buffer <- u_dist_nearest_border < dist
    if(any(outside_buffer)){
      ind_outside_buff <- which(outside_buffer)
      new_width <- ceiling(max(dist[ind_outside_buff]))
      new_ranks <- round(ranks[order(ranks[,1]),],0)
      new_ranks <- rbind(tail(new_ranks,new_width),new_ranks,head(new_ranks,new_width))
      new_ranks[,1] <- new_ranks[,1]+c(rep(-N,new_width),rep(0,N),rep(N,new_width))+new_width
      new_points_2_eval <- (seq_len(N) + new_width)[ind_outside_buff]


      # if(metric=="max"){
      #   new_rotate_ranks_45 <- as.data.frame(rotate(new_ranks,45))
      #   new_neighbors <-  RANN1::nn2(new_rotate_ranks_45,
      #                                new_rotate_ranks_45[new_points_2_eval,],
      #                                k=(k+1))
      #   dist_new <- round(new_neighbors$nn.dists[,(k+1)]/sqrt(1^2+1^2),0)
      # }else if(metric=="euclid"){
        new_neighbors <-  dbscan::kNN(x=new_ranks,
                                      query = new_ranks[new_points_2_eval,],
                                      k=k+1,
                                      bucketSize = 10,
                                      sort=F,
                                      search="kdtree",
                                      splitRule="SL_MIDPT")
        dist_new <- new_neighbors$dist[,(k+1)]
      # }
      dist[ind_outside_buff] <- dist_new
    }
  }

  mi <- -(digamma(N) - digamma(k)+log(cd) +2/N*sum(log(dist)) - 2*log(N+1))

  if(mi<0){
    warning(warning_sound(), paste0("MI is negative! With the current number of samples, the bias can be up
          to about ", N^(-0.5)))
  }

  return(mi)
}




#' Estimate the Mutual Information of a Copula.
#'
#' The mutual information is obtained from the negative
#' differential entropy of the copula PDF using a Monte Carlo approach.
#'
#'
#' @param copula '\code{\linkS4class{cyl_copula}}' or a '\code{\linkS4class{Copula}}' object
#' from the package '\pkg{copula}'.
#' @param n \link[base]{integer} value, the number of samples to use in the Monte
#' Carlo approach.
#'
#'
#' @return A \link[base]{numeric} value, the mutual information of the copula \code{cop}
#' in nats.
#'
#' @examples
#' set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' mi_cyl_mc(cop, n = 100)
#'
#' @references
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cor_cyl}()}, \code{\link{mi_cyl_knn}()}, \code{\link{mi_cyl_bin}()}.
#'
#' @export
#'
mi_cyl_mc <- function(cop,
                      n = 100000) {


  tryCatch({
    check_arg_all(list(
      check_argument_type(copula,
                          type = "Copula"),
      check_argument_type(copula,
                          type = "cyl_copula")
    )
    , 2)
    check_arg_all(check_argument_type(n, type="numeric", length = 1, lower = 1, integer=T)
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )

  samp <- rcylcop(n,cop)
  log_cop_samp <- log(dcylcop(samp,copula = cop))
  entrop <- -sum(log_cop_samp)/n
  entrop_var <- sum((log_cop_samp-entrop)^2)/(n+1)
  sqrt(entrop_var)
  entrop

  return(mi)
}

#' Calculate a non-linear correlation coefficient from MI
#'
#' From the mutual information of 2 univariate random variables, the Pearson
#' correlation coefficient of a bivariate Gaussian distribution with
#' the same mutual information is calculated.
#'
#'
#' @param mi \link[base]{numeric} value, the (unnormalized!) mutual information
#' obtained e.g. with  \code{\link{mi_cyl_knn}()}, \code{\link{mi_cyl_bin}()}, or
#' \code{\link{mi_cyl_mc}()}.
#'
#' @return A \link[base]{numeric} value, the non-linear correlation coefficient
#'
#' @examples
#' set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' mi <- mi_cyl_mc(cop, n = 100)
#' mi2cc(mi)
#'
#' @references
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{{mi_cyl_mc}()}, \code{\link{mi_cyl_knn}()}, \code{\link{mi_cyl_bin}()}.
#'
#' @export
#'
#'
mi2cc <- function(mi){
  tryCatch({
    check_arg_all(
      check_argument_type(mi,
                          type = "numeric",
                          lower = 0
                          ), 1)
    },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  })

  sqrt(1-exp(-2*mi))
}


#' Calculate a Non-Linear Correlation Coefficient Using Optimal Transport
#'
#'
#'
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles (measurements of a circular
#'   variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths (measurements of a linear
#'   variable).
#' @param nbins \link[base]{integer} value, the number of bins to use in each dimension.
#' for the calculation of the optimal transport plan.
#'
#' @details The function estimates the Wasserstein distance between the empirical copula
#' of the data and the independence copula and
#' the empirical copula and the distance to the closest Frechet-Hoeffding
#' copula. The correlation coefficient is then the ratio between the distance to the
#' independence copula and the sum of the distances between the empirical copula
#' and the independence copula and the empirical copula and the closest
#' Frechet-Hoeffding copula.
#'
#' @return A \link[base]{numeric} value, non-linear correlation between \code{theta} and \code{x}.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#' marg1 <- list(name="vonmises",coef=list(0,4))
#' marg2 <- list(name="lnorm",coef=list(2,3))
#'
#' #draw samples and calculate the mutual information.
#' sample <- rjoint(100,cop,marg1,marg2)
#' cor_wasser(
#'   theta = sample[,1],
#'   x = sample[,2],
#'   nbins = 5
#' )
#'
#' @references
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cor_cyl}()}, \code{\link{mi2cc}()}.
#'
#' @export
#'
cor_wasser <- function(theta,
                       x,
                       nbins = 50){

  tryCatch({
    check_arg_all(check_argument_type(theta, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(nbins, type="numeric", length = 1, lower = 1, integer=T)
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

  pi_d <- wasserstein(copula = copula::indepCopula(),theta=theta,x=x,emp_binned = T,
                           nbins=nbins,method="aha")
  pi_d/(pi_d+min(wasserstein(copula = copula::upfhCopula(),theta=theta,x=x,emp_binned = T,
                             nbins=nbins,method="aha"),
                 wasserstein(copula = copula::lowfhCopula(),theta=theta,x=x,emp_binned = T,
                             nbins=nbins,method="aha")))
}
