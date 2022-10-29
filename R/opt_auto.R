
#' Automatically Find the Best Fitting Copula
#'
#' The parameters of 15 different circular-linear copulas are fitted to data
#' and sorted
#' according to AIC. For each copula, first, a starting value for the maximum
#' likelihood estimation (MLE) is found using \code{\link{fit_cylcop_cor}()}.
#' Then, MLE is carried out with a "reasonable" setup using \code{\link{fit_cylcop_ml}()}.
#' If MLE fails, parameters obtained with \code{\link{fit_cylcop_cor}()} are reported.
#'
#' @param theta \link[base]{numeric} \link[base]{vector} of angles
#' (measurements of a circular variable).
#' @param x \link[base]{numeric} \link[base]{vector} of step lengths
#' (measurements of a linear variable).
#'
#' @return A list containing 3 lists: Descriptions of the copulas, the
#' '\code{\linkS4class{cyl_copula}}' objects with fitted parameters, and the AIC.
#' The lists are sorted by ascending AIC.
#' If \code{\link{fit_cylcop_ml}()} has failed, the reported parameters are the ones obtained
#'  with \code{\link{fit_cylcop_cor}()} and the AIC is set to \code{NA}.
#'
#' @examples set.seed(123)
#'
#' #Optimal copula is independent of marginals.
#' data <- rcylcop(100,cyl_quadsec(0.1))
#'
#' #This takes a few seconds to run.
#' \donttest{copula_lst <- opt_auto(theta = data[,1], x = data[,2])}
#'
#' @seealso \code{\link{fit_cylcop_cor}()}, \code{\link{fit_cylcop_ml}()}
#'
#' @references \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
#' @export
#'
opt_auto <- function(theta, x) {
  #validate input
  tryCatch({

    check_arg_all(check_argument_type(theta,
                                           type="numeric")
                  ,1)
    check_arg_all(check_argument_type(x,
                                      type="numeric")
                  ,1)
  },
  error = function(e) {
    error_sound()
    rlang::abort(conditionMessage(e))
  }
  )
  original_silent <- cylcop.env$silent
  cylcop.env$silent <- T

  name_lst <- vector(mode = "list", length = 15)
  copula_lst <- vector(mode = "list", length = 15)
  AIC_lst <- vector(mode = "list", length = 15)

  if (original_silent == FALSE) {
    message("Optimizing copula 1, cyl_vonmises")
  }
  out <- calc_vonmises(theta, x, flip=FALSE)
  name_lst[[1]] <- out$name
  copula_lst[[1]] <- out$copula
  AIC_lst[[1]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 2, cyl_vonmises")
  }
  out <- calc_vonmises(theta, x, flip=TRUE)
  name_lst[[2]] <- out$name
  copula_lst[[2]] <- out$copula
  AIC_lst[[2]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 3, cyl_quadsec")
  }
  out <- calc_quadsec(theta, x)
  name_lst[[3]] <- out$name
  copula_lst[[3]] <- out$copula
  AIC_lst[[3]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 4, cyl_cubsec")
  }
  out <- calc_cubsec(theta, x)
  name_lst[[4]] <- out$name
  copula_lst[[4]] <- out$copula
  AIC_lst[[4]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 5, cyl_rot_combine(frankCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::frankCopula(2), shift = FALSE)
  name_lst[[5]] <- out$name
  copula_lst[[5]] <- out$copula
  AIC_lst[[5]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 6, cyl_rot_combine(claytonCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::claytonCopula(2), shift = FALSE)
  name_lst[[6]] <- out$name
  copula_lst[[6]] <- out$copula
  AIC_lst[[6]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 7, cyl_rot_combine(gumbelCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::gumbelCopula(4), shift = FALSE)
  name_lst[[7]] <- out$name
  copula_lst[[7]] <- out$copula
  AIC_lst[[7]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 8, cyl_rot_combine(frankCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::frankCopula(2), shift = TRUE)
  name_lst[[8]] <- out$name
  copula_lst[[8]] <- out$copula
  AIC_lst[[8]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 9, cyl_rot_combine(claytonCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::claytonCopula(2), shift = TRUE)
  name_lst[[9]] <- out$name
  copula_lst[[9]] <- out$copula
  AIC_lst[[9]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 10, cyl_rot_combine(gumbelCopula())")
  }
  out <- calc_rot_combine(theta, x, lin_cop=copula::gumbelCopula(4), shift = TRUE)
  name_lst[[10]] <- out$name
  copula_lst[[10]] <- out$copula
  AIC_lst[[10]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 11, cyl_rect_combine(frankCopula())")
  }
  out <- calc_rect_combine(theta, x, lin_cop=copula::frankCopula(2), flip_up = TRUE)
  name_lst[[11]] <- out$name
  copula_lst[[11]] <- out$copula
  AIC_lst[[11]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 12, cyl_rect_combine(claytonCopula())")
  }
  out <- calc_rect_combine(theta, x, lin_cop=copula::claytonCopula(2), flip_up = TRUE)
  name_lst[[12]] <- out$name
  copula_lst[[12]] <- out$copula
  AIC_lst[[12]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 13, cyl_rect_combine(gumbelCopula())")
  }
  out <- calc_rect_combine(theta, x, lin_cop=copula::gumbelCopula(2), flip_up = TRUE)
  name_lst[[13]] <- out$name
  copula_lst[[13]] <- out$copula
  AIC_lst[[13]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 14, cyl_rect_combine(claytonCopula())")
  }
  out <- calc_rect_combine(theta, x, lin_cop=copula::claytonCopula(2), flip_up = FALSE)
  name_lst[[14]] <- out$name
  copula_lst[[14]] <- out$copula
  AIC_lst[[14]] <- out$AIC

  if (original_silent == FALSE) {
    message("Optimizing copula 15, cyl_rect_combine(gumbelCopula())")
  }
  out <- calc_rect_combine(theta, x, lin_cop=copula::gumbelCopula(2), flip_up = FALSE)
  name_lst[[15]] <- out$name
  copula_lst[[15]] <- out$copula
  AIC_lst[[15]] <- out$AIC


  cylcop.env$silent <- original_silent

  if (any(is.infinite(unlist(AIC_lst)))) {
    warning(warning_sound(),
            "For at least one copula, MLE did not converge!")
  }

  ordered_ind <- order(unlist(AIC_lst))
  name_lst_ordered <- name_lst
  copula_lst_ordered <- copula_lst
  AIC_lst_ordered <- AIC_lst


  for (i in 1:length(name_lst)) {
    name_lst_ordered[[i]] <- name_lst[[ordered_ind[i]]]
    copula_lst_ordered[[i]] <- copula_lst[[ordered_ind[i]]]
    if (is.infinite(AIC_lst[[ordered_ind[i]]])) {
      AIC_lst_ordered[[i]] <- "MLE NOT CONVERGED"
    } else{
      AIC_lst_ordered[[i]] <- AIC_lst[[ordered_ind[i]]]
    }
  }
  return(list(
    name = name_lst_ordered,
    copula = copula_lst_ordered,
    AIC = AIC_lst_ordered
  ))
}


calc_vonmises <- function(theta, x, flip){
  out_lst <- list()
  start_val <-
    tryCatch(
      fit_cylcop_cor(
        copula = cyl_vonmises(flip = flip),
        theta = theta,
        x = x,
        method = "cor_cyl",
        acc = 1,
        n = 1000
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      fit_cylcop_ml(
        copula = cyl_vonmises(flip = flip),
        theta = theta,
        x = x,
        parameters = c("kappa", "mu"),
        start = c(start_val, 0),
        lower = c(0, -pi),
        upper = c(1000, pi)
      ),
      error = function(e) {
        return(0)
      },
      warning = function(e) {
        return(0)
      }
    )

  if (is.numeric(copula)) {
    if (isFALSE(start_val)) {
      out_lst$copula <- "ERROR"
    } else{
      out_lst$copula <- cyl_vonmises(mu = 0,
                                     kappa = start_val,
                                     flip = flip)
    }
    out_lst$AIC <- Inf
  } else{
    out_lst$copula <- copula$copula
    out_lst$AIC <- copula$AIC
  }
  if(!flip){
    out_lst$name <- "cyl_vonmises, positive correlation"
  } else{
    out_lst$name <- "cyl_vonmises, negative correlation"
  }

  return(out_lst)
}


calc_quadsec <- function(theta, x){
  out_lst <- list()
start_val <- tryCatch(
  fit_cylcop_cor(copula = cyl_quadsec(0),
         theta, x, method = "cor_cyl"),
  error = function(e) {
    return(FALSE)
  }
)
copula1 <- tryCatch(
  fit_cylcop_ml(
    copula = cyl_quadsec(),
    theta = theta,
    x = x,
    parameters = "a",
    start = abs(start_val),
    lower = 0
  ),
  error = function(e) {
    return(0)
  },
  warning = function(e) {
    return(0)
  }
)
copula2 <- tryCatch(
  fit_cylcop_ml(
    copula = cyl_quadsec(),
    theta = theta,
    x = x,
    parameters = "a",
    start = -abs(start_val),
    upper = 0
  ),
  error = function(e) {
    return(0)
  },
  warning = function(e) {
    return(0)
  }
)

if (is.numeric(copula1) && is.numeric(copula2)) {
  if (isFALSE(start_val)) {
    out_lst$copula <- "ERROR"
  } else{
    out_lst$copula <- cyl_quadsec(start_val)
  }
  out_lst$AIC <- Inf
} else if (is.numeric(copula1) && is.list(copula2)) {
  copula <- copula2
} else if (is.numeric(copula2) && is.list(copula1)) {
  copula <- copula1
} else{
  if (copula1$AIC < copula2$AIC) {
    copula <- copula1
  } else{
    copula <- copula2
  }
  out_lst$copula <- copula$copula
  out_lst$AIC <- copula$AIC
}
out_lst$name <- "cyl_quadsec"
return(out_lst)
}

calc_cubsec <- function(theta, x){
  out_lst <- list()
  start_val <-
    tryCatch(
      fit_cylcop_cor(
        copula = cyl_cubsec(),
        theta = theta,
        x = x,
        method = "cor_cyl",
        acc = 0.05,
        n = 2000
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      fit_cylcop_ml(
        copula = cyl_cubsec(),
        theta = theta,
        x = x,
        parameters = c("a", "b"),
        start = start_val
      ),
      error = function(e) {
        return(0)
      },
      warning = function(e) {
        return(0)
      }
    )
  if (is.numeric(copula)) {
    if (isFALSE(start_val)) {
      out_lst$copula <- "ERROR"
    } else{
      out_lst$copula <- cyl_cubsec(a = start_val[1], b = start_val[2])
    }
    out_lst$AIC <- Inf
  } else{
    out_lst$copula <- copula$copula
    out_lst$AIC <- copula$AIC
  }
  out_lst$name <- "cyl_cubsec"
  return(out_lst)
}


calc_rot_combine <- function(theta, x, lin_cop, shift){
  out_lst <- list()
  start_val <-
    tryCatch(
      fit_cylcop_cor(
        copula = cyl_rot_combine(lin_cop, shift = shift),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  name_text <- "cyl_rot_combine,"
  if("frankCopula" %in% is (lin_cop)){
    lower <- NULL
    upper <- NULL
    name_text <- paste(name_text, "Frank copula,")
  }
  if("claytonCopula" %in% is (lin_cop)){
    lower <- 0.0001
    upper <- NULL
    name_text <- paste(name_text, "Clayton copula,")
  }
  if("gumbelCopula" %in% is (lin_cop)){
    lower <- 1
    upper <- start_val + 30
    name_text <- paste(name_text, "Gumbel copula,")
  }
  copula <-
    tryCatch(
      fit_cylcop_ml(
        copula = cyl_rot_combine(lin_cop, shift = shift),
        theta = theta,
        x = x,
        parameters = c("alpha"),
        start = start_val,
        lower = lower,
        upper = upper
      ),
      error = function(e) {
        return(0)
      },
      warning = function(e) {
        return(0)
      }
    )

  if (is.numeric(copula)) {
    if (isFALSE(start_val)) {
      out_lst$copula <- "ERROR"
    } else{
      lin_cop@parameters <- start_val
      out_lst$copula <-
        cyl_rot_combine(lin_cop, shift = shift)
    }
    out_lst$AIC <- Inf
  } else{
    if (("claytonCopula" %in% is (lin_cop)) && (copula$copula@parameters[1] <= 0.0005)) {
      copula$copula@parameters[1] <- 0
      out_lst$copula <- copula$copula
      out_lst$AIC <- 2
    } else{
      out_lst$copula <- copula$copula
      out_lst$AIC <- copula$AIC
    }
  }

  if(shift){
    out_lst$name <- paste(name_text, "'diamond'-shape")
  } else{
    out_lst$name <- paste(name_text, "'X'-shape")
  }
  return(out_lst)
}


calc_rect_combine <- function(theta, x, lin_cop, flip_up){
  out_lst <- list()
  start_val <-
    tryCatch(
      suppressWarnings(optTau(
        copula = cyl_rect_combine(
          copula = lin_cop,
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = flip_up
        ),
        theta = theta,
        x = x
      )),
      error = function(e) {
        return(FALSE)
      },
      warning = function(e) {
        return(FALSE)
      }
    )
  name_text <- "cyl_rect_combine,"
  if("frankCopula" %in% is (lin_cop)){
    lower <- NULL
    upper <- NULL
    name_text <- paste(name_text, "Frank copula,")
  }
  if("claytonCopula" %in% is (lin_cop)){
    lower <- 0.0001
    upper <- NULL
    name_text <- paste(name_text, "Clayton copula,")
  }
  if("gumbelCopula" %in% is (lin_cop)){
    lower <- 1
    upper <- start_val + 30
    name_text <- paste(name_text, "Gumbel copula,")
  }
  copula <-
    tryCatch(
      fit_cylcop_ml(
        copula = cyl_rect_combine(
          copula = lin_cop,
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = flip_up
        ),
        theta = theta,
        x = x,
        parameters = "alpha",
        start = start_val,
        upper = upper,
        lower = lower
      ),
      error = function(e) {
        return(0)
      },
      warning = function(e) {
        return(0)
      }
    )


  if (is.numeric(copula)) {
    if (isFALSE(start_val)) {
      out_lst$copula <- "ERROR"
    } else{
      lin_cop@parameters <- start_val
      out_lst$copula <-
        cyl_rect_combine(
          copula = lin_cop,
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = flip_up
        )
    }
    out_lst$AIC <- Inf
  } else{
    if (("claytonCopula" %in% is (lin_cop)) && (copula$copula@parameters[1] <= 0.0005)) {
      copula$copula@parameters[1] <- 0
      out_lst$copula <- copula$copula
      out_lst$AIC <- 2
    } else{
      out_lst$copula <- copula$copula
      out_lst$AIC <- copula$AIC
    }
  }

  if(flip_up){
    out_lst$name <- paste(name_text, "'>'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]")
  } else{
    out_lst$name <- paste(name_text,  "'<'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]")
  }
  return(out_lst)
}









