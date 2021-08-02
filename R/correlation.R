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
#' sample <- make_traj(100,
#'   cop,
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape = 2)
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
#' @seealso \code{\link{mi_cyl}()}, \code{\link{optCor}()}.
#'
#' @export
#'
cor_cyl <- function(theta, x) {

#Variable names as in Solow 1988 but with radians instead of degrees.

# Assigning ranks to angular and linear measurements

  data <- data.frame(theta, x) %>% na.omit() %>% dplyr::arrange(x) %>%
    mutate(r_theta = frank(.data$theta, ties.method = "average"))
  n <- nrow(data)
  data <- mutate(data, r_theta_star = .data$r_theta * 2 * pi / n) %>%
    mutate(r_x = frank(.data$x, ties.method = "average"))


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
#' The mutual information can be normalized to lie between 0 ans 1
#' by dividing by the product of the entropies of \code{x} and \code{theta}.
#' Even if \code{x} and \code{theta} are perfectly correlated, the normalized
#' mutual information will not be 1 if the underlying copula is periodic and
#' symmetric. Therefore, we can set \code{symmetrize = TRUE} to set all u-values of
#' the empirical copula that are larger than \eqn{0.5} to \eqn{1-0.5}. The mutual information
#' is then calculated from those values and the is exactly 1 in the case of
#' perfect correlation as captured by e.g.
#' \code{cyl_rect_combine(normalCopula(1))}. The estimate (output of
#' \code{mi_cyl()}) will be less than one for numerical reasons. Note also that
#' the mutual information is independent of the marginal distributions.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles (measurements of a circular
#'   variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths (measurements of a linear
#'   variable).
#' @param normalize \link[base]{logical} value whether the mutual information should be
#'   normalized to lie within \eqn{[0,1]}.
#' @param symmetrize \link[base]{logical} value whether it should be assumed that positive
#'   and negative angles are equivalent.
#'
#' @return A \link[base]{numeric} value, the mutual information between \code{theta} and \code{x}.
#'
#' @examples set.seed(123)
#'
#' cop <- cyl_quadsec(0.1)
#'
#' #draw samples and calculate the mutual information.
#' sample <- rcylcop(100, cop)
#' mi_cyl(theta = sample[,1],
#'   x = sample[,2],
#'   normalize = TRUE,
#'   symmetrize = FALSE
#' )
#'
#' #the correlation coefficient is independent of the marginal distribution.
#' sample <- make_traj(100,
#'   cop,
#'   marginal_circ = "vonmises",
#'   parameter_circ = list(0, 1),
#'   marginal_lin = "weibull",
#'   parameter_lin = list(shape = 2)
#' )
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
#' sample <- rcylcop(100, cop)
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
#' #only with normaliztion and symmetrization do we get a value close to 1
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
#' \insertRef{Tenzer2016}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @seealso \code{\link{cor_cyl}()}, \code{\link{optCor}()}.
#'
#' @export
#'
mi_cyl <- function(theta,
                      x,
                      normalize = T,
                      symmetrize = F) {
  data <- data.frame(theta, x) %>% na.omit()

  #calculate empirical copula
  emp_cop <- pobs(data, ties.method = "average")%>%as.data.frame()

  if (symmetrize) {
    emp_cop$theta <- modify_if(emp_cop$theta, ~ .x > 0.5, ~ 1 - .x)
  }
  discr_theta <- infotheo::discretize(emp_cop$theta)
  discr_x <- infotheo::discretize(emp_cop$x)
  if (normalize) {
    mi <-
      infotheo::mutinformation(discr_theta,discr_x) / sqrt(infotheo::entropy(discr_theta) * infotheo::entropy(discr_x))
  }

  else{
    mi <- infotheo::mutinformation(discr_theta,discr_x)
  }

  return(mi)
}
