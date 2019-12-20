#' Estimate a rank-based circular-linear correlation-coefficient-like parameter
#'
#' Nomenclature afer Solow 1988 but with radians instead of degrees.
#'
#' @param theta A numeric vector of angles (measurements of a circular variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#' @param plot A logical value whether a plot of the regression line should be dislpayed.
#'
#' @return The function returns a numeric value between 0 and 1, not -1 and 1, hence "correlation-coefficient-LIKE"
#' @export
#'
cor_cyl <- function(theta, x, plot = TRUE) {

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

# Plot regression line
  if (plot) {
    lm_test <- lm(r_x ~ cos(r_theta_star) + sin(r_theta_star), data = data)
    reg_line <- data.frame(theta = seq(0, 2 * pi, 0.01))
    reg_line <-
      mutate(
        reg_line,
        x = as.double(coef(lm_test)[1]) + as.double(coef(lm_test)[2]) * cos(theta) +
          as.double(coef(lm_test)[3]) * sin(theta)
      )
    p <- ggplot() +
      geom_point(data = data, aes(x = .data$r_x, y = .data$r_theta_star)) +
      geom_point(data = reg_line, aes(x = x, y = theta), color = "red") +
      theme_bw() +
      xlab("linear rank") +
      ylab("circular rank") +
      theme(
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none"
      )
    plot(p)
  }

  return(D)
}



#' Estimate the mutual information between a circular and a linear random
#' variable
#'
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear
#'   variable).
#' @param normalize A logical value whether the mutual information should be
#'   normalized to lie within [0,1].
#' @param symmetrize A logical value whether it should be assumed that positive
#'   and negative angles are equivalent. If they are indeed exactly equivalent
#'   and \code{symmetrize=F}, the maximal normalized mutual information will be about
#'   0.5.
#'
#' @return A numeric value, the mutual information between \code{theta} and \code{x}.
#' @export
#'
mi_binned <- function(theta,
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
