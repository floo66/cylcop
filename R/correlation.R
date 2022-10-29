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
#' @seealso \code{\link{mi_cyl}()}, \code{\link{fit_cylcop_cor}()}.
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
#' Variable
#'
#' The empirical copula is obtained from the data (\code{theta} and \code{x}),
#' and the mutual information of the 2 components is calculated. This gives a
#' non-negative number that can be normalized to lie between 0 and 1.
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
#'
#' @details First, the two components of the empirical copula, \eqn{u} and \eqn{v}
#' are obtained. Then the mutual information is calculated via discretizing \eqn{u} and \eqn{v}
#' into \code{length(theta)^(1/3)} bins. The mutual information can be
#' normalized to lie between 0 and 1 by dividing by the product of the entropies
#' of \code{u} and \code{v}. This is done using functions from the '\pkg{infotheo}'
#' package.
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
#' mi_cyl(theta = sample[,1],
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
#' mi_cyl(theta = sample$angle,
#'   x = sample$steplength,
#'   normalize = TRUE,
#'   symmetrize = FALSE)
#' mi_cyl(theta = sample$cop_u,
#'   x = sample$cop_v,
#'   normalize = TRUE,
#'   symmetrize = FALSE)
#'
#' # Estimate correlation of samples drawn from circular-linear copulas
#' # with perfect correlation.
#' cop <- cyl_rect_combine(copula::normalCopula(1))
#' sample <- rjoint(100,cop,marg1,marg2)
#' # without normalization
#' mi_cyl(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = FALSE,
#'   symmetrize = FALSE
#' )
#' #with normalization
#' mi_cyl(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = TRUE,
#'   symmetrize = FALSE
#' )
#' #only with normalization and symmetrization do we get a value of 1
#' mi_cyl(theta = sample[,1],
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
mi_cyl <- function(theta,
                      x,
                      normalize = TRUE,
                      symmetrize = FALSE) {

  tryCatch({
    check_arg_all(check_argument_type(theta, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x, type="numeric")
                  ,1)
    check_arg_all(check_argument_type(normalize, type="logical")
                  ,1)
    check_arg_all(check_argument_type(symmetrize, type="logical")
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

  nbins <-max(2,nrow(emp_cop)^(1/3))

  discr_theta <- infotheo::discretize(emp_cop[,1],nbins=nbins)
  discr_x <- infotheo::discretize(emp_cop[,2],nbins=nbins)
  if (normalize) {
    mi <-
      infotheo::mutinformation(discr_theta,discr_x) / sqrt(infotheo::entropy(discr_theta) * infotheo::entropy(discr_x))
  }

  else{
    mi <- infotheo::mutinformation(discr_theta,discr_x)
  }

  return(mi)
}
