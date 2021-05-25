#' Estimate a rank-based circular-linear correlation-coefficient-like parameter
#'
#' The code is based on \insertCite{Mardia1976;textual}{cylcop}, \insertCite{Solow1988;textual}{cylcop}
#' and \insertCite{Tu2015;textual}{cylcop}. The function returns a numeric value
#' between 0 and 1, not -1 and 1, hence "correlation-coefficient-LIKE".
#'
#' @param theta A numeric vector of angles (measurements of a circular variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#'
#' @return A numeric value between 0 and 1
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
  else{symmetrize
    a <- 2 * (sin(pi / n)) ^ 4 / ((1 + (cos(pi / n))) ^ 3)
  }

  D <- a * (C ^ 2 + S ^ 2)

  return(D)
}



#' Estimate the mutual information between a circular and a linear random
#' variable
#'
#' The mutual information can be normalized to lie between 0 ans 1
#' by dividing by the product of the entropies of \code{x} and \code{theta}.
#' Even if \code{x} and \code{theta} are perfectly correlated, the normalized
#' mutual information will not be 1 if the underlying copula is periodic and
#' symmetric. Therefore, we can set \code{symmetrize=T} to set all u-values of
#' the empirical copula that are larger than 0.5 to 1-0.5. The mutual information
#' is then calculated from those values and is exactly 1 in the case of
#' perfect correlation as captured by e.g.
#' \code{cyl_rect_combine(normalCopula(1))}.
#'
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear
#'   variable).
#' @param normalize A logical value whether the mutual information should be
#'   normalized to lie within [0,1].
#' @param symmetrize A logical value whether it should be assumed that positive
#'   and negative angles are equivalent.
#'
#' @return A numeric value, the mutual information between \code{theta} and \code{x}.
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
