
#' Black-box to find best fitting copula
#'
#' The parameters of 15 different circular-linear copulas are fitted and sorted
#' according to AIC. For each copula, first, a starting value for the maximum
#' likelihood estimation (MLE) is found using \code{cylcop::optCor()} or
#' \code{cylcop::optTau()}. Then MLE with a "reasonable" setup is carried out.
#' If MLE fails, parameters obtained with \code{cylcop::optCor()} or
#' \code{cylcop::optTau()} are reported.
#'
#' @param theta A numeric vector of angles (measurements of a circular
#'   variable).
#' @param x A numeric vector of steplengths (measurements of a linear variable).
#'
#' @return A list containing 3 lists: Descriptions of the copulae, the
#' \code{cyl_copula}-objects with fitted parameters, the AIC.
#' @export
#'
opt_auto <- function(theta, x) {
  original_silent <- cylcop.env$silent
  cylcop.env$silent <- T
  original_silent <- F

  name_lst <- vector(mode = "list", length = 15)
  copula_lst <- vector(mode = "list", length = 15)
  AIC_lst <- vector(mode = "list", length = 15)

  if (original_silent == F) {
    message("Optimizing copula 1, cyl_vonmises")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_vonmises(flip = F),
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
      optML(
        copula = cyl_vonmises(flip = F),
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
      copula_lst[[1]] <- "ERROR"
    } else{
      copula_lst[[1]] <- cyl_vonmises(mu = 0,
                                      kappa = start_val,
                                      flip = F)
    }
    AIC_lst[[1]] <- Inf
  } else{
    copula_lst[[1]] <- copula$copula
    AIC_lst[[1]] <- copula$AIC
  }
  name_lst[[1]] <- "cyl_vonmises, positive correlation"

  if (original_silent == F) {
    message("Optimizing copula 2, cyl_vonmises")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_vonmises(flip = T),
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
      optML(
        copula = cyl_vonmises(flip = T),
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
      copula_lst[[2]] <- "ERROR"
    } else{
      copula_lst[[2]] <- cyl_vonmises(mu = 0,
                                      kappa = start_val,
                                      flip = T)
    }
    AIC_lst[[2]] <- Inf
  } else{
    copula_lst[[2]] <- copula$copula
    AIC_lst[[2]] <- copula$AIC
  }
  name_lst[[2]] <- "cyl_vonmises, negative correaltion"

  if (original_silent == F) {
    message("Optimizing copula 3, cyl_quadsec")
  }
  start_val <- tryCatch(
    cylcop::optCor(copula = cyl_quadsec(0),
                   theta, x, method = "cor_cyl"),
    error = function(e) {
      return(FALSE)
    }
  )
  copula1 <- tryCatch(
    cylcop::optML(
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
    cylcop::optML(
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
      copula_lst[[3]] <- "ERROR"
    } else{
      copula_lst[[3]] <- cyl_quadsec(start_val)
    }
    AIC_lst[[3]] <- Inf
  } else if (is.numeric(copula1) && is.list(copula2)) {
    copula == copula2
  } else if (is.numeric(copula2) && is.list(copula1)) {
    copula == copula1
  } else{
    if (copula1$AIC < copula2$AIC) {
      copula <- copula1
    } else{
      copula <- copula2
    }
    copula_lst[[3]] <- copula$copula
    AIC_lst[[3]] <- copula$AIC
  }
  name_lst[[3]] <- "cyl_quadsec"


  if (original_silent == F) {
    message("Optimizing copula 4, cyl_cubsec")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_cubsec(),
        theta = theta,
        x = x,
        method = "cor_cyl",
        acc = 0.05,
        n = 5000
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
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
      copula_lst[[4]] <- "ERROR"
    } else{
      copula_lst[[4]] <- cyl_cubsec(a = start_val[1], b = start_val[2])
    }
    AIC_lst[[4]] <- Inf
  } else{
    copula_lst[[4]] <- copula$copula
    AIC_lst[[4]] <- copula$AIC
  }
  name_lst[[4]] <- "cyl_cubsec"


  if (original_silent == F) {
    message("Optimizing copula 5, cyl_rot_combine(frankCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(frankCopula(2), shift = F),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(frankCopula(2), shift = F),
        theta = theta,
        x = x,
        parameters = c("alpha"),
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
      copula_lst[[5]] <- "ERROR"
    } else{
      copula_lst[[5]] <-
        cyl_rot_combine(frankCopula(start_val), shift = F)
    }
    AIC_lst[[5]] <- Inf
  } else{
    copula_lst[[5]] <- copula$copula
    AIC_lst[[5]] <- copula$AIC
  }
  name_lst[[5]] <- "cyl_rot_combine, Frank copula, 'X'-shape"


  if (original_silent == F) {
    message("Optimizing copula 6, cyl_rot_combine(claytonCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(claytonCopula(2), shift = F),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(claytonCopula(2), shift = F),
        theta = theta,
        x = x,
        parameters = c("alpha"),
        start = start_val,
        lower = 0.0001
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
      copula_lst[[6]] <- "ERROR"
    } else{
      copula_lst[[6]] <-
        cyl_rot_combine(claytonCopula(start_val), shift = F)
    }
    AIC_lst[[6]] <- Inf
  } else{
    if (copula$copula@parameters[1] <= 0.0005) {
      copula$copula@parameters[1] <- 0
      copula_lst[[6]] <- copula$copula
      AIC_lst[[6]] <- 2
    } else{
      copula_lst[[6]] <- copula$copula
      AIC_lst[[6]] <- copula$AIC
    }
  }
  name_lst[[6]] <- "cyl_rot_combine, Clayton copula, 'X'-shape"


  if (original_silent == F) {
    message("Optimizing copula 7, cyl_rot_combine(gumbelCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(gumbelCopula(4), shift = F),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(gumbelCopula(4), shift = F),
        theta = theta,
        x = x,
        parameters = c("alpha"),
        start = start_val,
        lower = 1,
        upper = start_val + 30
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
      copula_lst[[7]] <- "ERROR"
    } else{
      copula_lst[[7]] <-
        cyl_rot_combine(gumbelCopula(start_val), shift = F)
    }
    AIC_lst[[7]] <- Inf
  } else{
    copula_lst[[7]] <- copula$copula
    AIC_lst[[7]] <- copula$AIC
  }
  name_lst[[7]] <- "cyl_rot_combine, Gumbel copula, 'X'-shape"

  if (original_silent == F) {
    message("Optimizing copula 8, cyl_rot_combine(frankCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(frankCopula(2), shift = T),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(frankCopula(2), shift = T),
        theta = theta,
        x = x,
        parameters = c("alpha"),
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
      copula_lst[[8]] <- "ERROR"
    } else{
      copula_lst[[8]] <-
        cyl_rot_combine(frankCopula(start_val), shift = T)
    }
    AIC_lst[[8]] <- Inf
  } else{
    copula_lst[[8]] <- copula$copula
    AIC_lst[[8]] <- copula$AIC
  }
  name_lst[[8]] <- "cyl_rot_combine, Frank copula, 'diamond'-shape"


  if (original_silent == F) {
    message("Optimizing copula 9, cyl_rot_combine(claytonCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(claytonCopula(2), shift = T),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(claytonCopula(2), shift = T),
        theta = theta,
        x = x,
        parameters = c("alpha"),
        start = start_val,
        lower = 0.0001
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
      copula_lst[[9]] <- "ERROR"
    } else{
      copula_lst[[9]] <-
        cyl_rot_combine(claytonCopula(start_val), shift = T)
    }
    AIC_lst[[9]] <- Inf
  } else{
    if (copula$copula@parameters[1] <= 0.0005) {
      copula$copula@parameters[1] <- 0
      copula_lst[[9]] <- copula$copula
      AIC_lst[[9]] <- 2
    } else{
      copula_lst[[9]] <- copula$copula
      AIC_lst[[9]] <- copula$AIC
    }
  }
  name_lst[[9]] <-
    "cyl_rot_combine, Clayton copula, 'diamond'-shape"


  if (original_silent == F) {
    message("Optimizing copula 10, cyl_rot_combine(gumbelCopula())")
  }
  start_val <-
    tryCatch(
      cylcop::optCor(
        copula = cyl_rot_combine(gumbelCopula(4), shift = T),
        theta,
        x,
        method = "mi_cyl"
      ),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rot_combine(gumbelCopula(4), shift = T),
        theta = theta,
        x = x,
        parameters = c("alpha"),
        start = start_val,
        lower = 1,
        upper = start_val + 30,
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
      copula_lst[[10]] <- "ERROR"
    } else{
      copula_lst[[10]] <-
        cyl_rot_combine(gumbelCopula(start_val), shift = T)
    }
    AIC_lst[[10]] <- Inf
  } else{
    copula_lst[[10]] <- copula$copula
    AIC_lst[[10]] <- copula$AIC
  }
  name_lst[[10]] <-
    "cyl_rot_combine, Gumbel copula, 'diamond'-shape"



  if (original_silent == F) {
    message("Optimizing copula 11, cyl_rect_combine(frankCopula())")
  }
  start_val <-
    tryCatch(
      suppressWarnings(cylcop::optTau(
        copula = cyl_rect_combine(
          copula = frankCopula(2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
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
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rect_combine(
          copula = frankCopula(param = 2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        ),
        theta = theta,
        x = x,
        parameters = "alpha",
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
      copula_lst[[11]] <- "ERROR"
    } else{
      copula_lst[[11]] <-
        cyl_rect_combine(
          copula = frankCopula(start_val),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        )
    }
    AIC_lst[[11]] <- Inf
  } else{
    copula_lst[[11]] <- copula$copula
    AIC_lst[[11]] <- copula$AIC
  }
  name_lst[[11]] <-
    "cyl_rect_combine, Frank copula, '>'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]"



  if (original_silent == F) {
    message("Optimizing copula 12, cyl_rect_combine(claytonCopula())")
  }
  start_val <-
    tryCatch(
      suppressWarnings(cylcop::optTau(
        copula = cyl_rect_combine(
          copula = claytonCopula(2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        ),
        theta = theta,
        x = x
      )),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rect_combine(
          copula = claytonCopula(param = 2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        ),
        theta = theta,
        x = x,
        parameters = "alpha",
        start = start_val,
        lower = 0.0001
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
      copula_lst[[12]] <- "ERROR"
    } else{
      copula_lst[[12]] <-
        cyl_rect_combine(
          copula = claytonCopula(start_val),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        )
    }
    AIC_lst[[12]] <- Inf
  } else{
    if (copula$copula@parameters[1] <= 0.0005) {
      copula$copula@parameters[1] <- 0
      copula_lst[[12]] <- copula$copula
      AIC_lst[[12]] <- 2
    } else{
      copula_lst[[12]] <- copula$copula
      AIC_lst[[12]] <- copula$AIC
    }
  }
  name_lst[[12]] <-
    "cyl_rect_combine, Clayton copula, '>'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]"


  if (original_silent == F) {
    message("Optimizing copula 13, cyl_rect_combine(gumbelCopula())")
  }
  start_val <-
    tryCatch(
      suppressWarnings(cylcop::optTau(
        copula = cyl_rect_combine(
          copula = gumbelCopula(2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        ),
        theta = theta,
        x = x
      )),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <- tryCatch(
    cylcop::optML(
      copula = cyl_rect_combine(
        copula = gumbelCopula(param = 2),
        low_rect = c(0, 0.5),
        up_rect = "symmetric"
      ),
      theta = theta,
      x = x,
      parameters = "alpha",
      start = start_val,
      upper = start_val + 30
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
      copula_lst[[13]] <- "ERROR"
    } else{
      copula_lst[[13]] <-
        cyl_rect_combine(
          copula = gumbelCopula(start_val),
          low_rect = c(0, 0.5),
          up_rect = "symmetric"
        )
    }
    AIC_lst[[13]] <- Inf
  } else{
    copula_lst[[13]] <- copula$copula
    AIC_lst[[13]] <- copula$AIC
  }
  name_lst[[13]] <-
    "cyl_rect_combine, Gumbel copula, '>'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]"



  if (original_silent == F) {
    message("Optimizing copula 14, cyl_rect_combine(claytonCopula())")
  }
  start_val <-
    tryCatch(
      suppressWarnings(cylcop::optTau(
        copula = cyl_rect_combine(
          copula = claytonCopula(2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = F
        ),
        theta = theta,
        x = x
      )),
      error = function(e) {
        return(FALSE)
      }
    )
  copula <-
    tryCatch(
      cylcop::optML(
        copula = cyl_rect_combine(
          copula = claytonCopula(param = 2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = F
        ),
        theta = theta,
        x = x,
        parameters = "alpha",
        start = start_val,
        lower = 0.0001
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
      copula_lst[[14]] <- "ERROR"
    } else{
      copula_lst[[14]] <-
        cyl_rect_combine(
          copula = claytonCopula(start_val),
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = F
        )
    }
    AIC_lst[[14]] <- Inf
  } else{
    if (copula$copula@parameters[1] <= 0.0005) {
      copula$copula@parameters[1] <- 0
      copula_lst[[14]] <- copula$copula
      AIC_lst[[14]] <- 2
    } else{
      copula_lst[[14]] <- copula$copula
      AIC_lst[[14]] <- copula$AIC
    }
  }
  name_lst[[14]] <-
    "cyl_rect_combine, Clayton copula, '<'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]"


  if (original_silent == F) {
    message("Optimizing copula 15, cyl_rect_combine(gumbelCopula())")
  }
  start_val <-
    tryCatch(
      suppressWarnings(cylcop::optTau(
        copula = cyl_rect_combine(
          copula = gumbelCopula(2),
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = F
        ),
        theta = theta,
        x = x
      )),
      error = function(e) {
        return(FALSE)
      }
    )

  copula <- tryCatch(
    cylcop::optML(
      copula = cyl_rect_combine(
        copula = gumbelCopula(param = 2),
        low_rect = c(0, 0.5),
        up_rect = "symmetric",
        flip_up = F
      ),
      theta = theta,
      x = x,
      parameters = "alpha",
      start = start_val,
      upper = start_val + 30
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
      copula_lst[[15]] <- "ERROR"
    } else{
      copula_lst[[15]] <-
        cyl_rect_combine(
          copula = gumbelCopula(start_val),
          low_rect = c(0, 0.5),
          up_rect = "symmetric",
          flip_up = F
        )
    }
    AIC_lst[[15]] <- Inf
  } else{
    copula_lst[[15]] <- copula$copula
    AIC_lst[[15]] <- copula$AIC
  }
  name_lst[[15]] <-
    "cyl_rect_combine, Gumbel copula, '<'-shape, rectangles [0,0]x[0.5,1] and [0.5,0]x[1,1]"


  cylcop.env$silent <- original_silent

  if (any(is.infinite(unlist(AIC_lst)))) {
    warning(cylcop::warning_sound(),
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
