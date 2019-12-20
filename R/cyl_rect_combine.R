#' @include cyl_cop_class.R
NULL



#' An S4 class of circular-linear copulae generated from rectangular patches.
#'
#' This class contains bivariate circular-linear copulae generated from
#' linear-linear bivariate \code{Copula} objects of the \code{copula}-package or circular-linear
#' copulae of class \code{cyl_copula}. 2 non-overlapping rectangles are laid over the unit
#' square, both have width 1 in v-direction. In the area covered by the first rectangle,
#' the one covered by the second and the area not covered by any rectangle all have different copulae.
#' Rectangle 2 contains the same copula as rectangle 1, but 90 degrees rotated. The copula regions
#' are combined in a way that the overall result on the entire unit square is lso a copula.
#' See text for more rigorous explantion.
#'
#' With appropriate choices of the rectangles (again, see text) this results in copulae
#' that are periodic in u-direction (and not in v-direction) and therefore are
#' circular-linear.
#'
#' When the 2 rectangles are each other's mirror image with respect to u=0.5,
#' the resulting overall copula is symmetric with respect to u=0.5. I.e. there is
#' symmetry between positive and negative angles.
#'
#' @section Objects from the Class: Objects are created by
#'   \code{cyl_rect_combine()}.
#'
#' @slot name A character string holding the name of the copula.
#' @slot parameters A numeric vector holding the parameter values.
#' @slot param.names A character vector holding the parameter names.
#' @slot param.lowbnd A numeric vector holding the lower bounds of the
#'   parameters.
#' @slot param.upbnd A numeric vector holding the upper bounds of the
#'   parameters.
#' @slot sym.cop A \code{cyl_vonmises} or \code{Copula} object. The copula in
#'   the rectangles.
#' @slot background.cop A \code{cyl_copula} or \code{Copula} object. The copula
#'   where no rectangles overlay the unit square.
#' @slot flip_up Logical value indicating whether the copula (\code{sym.cop}) is
#'   rotated 90 degrees in the upper or lower rectangle.
#' @slot sym_rect Logical value indicating whether the upper rectangle was
#'   forced to be a mirror image of the lower one with respect to u=0.5 at the
#'   construction of the object.
#'
#' @section Extends:
#' Class \code{cyl_rect_combine} extends class \code{cyl_copula}.
#'
#' @export
#'
setClass(
  "cyl_rect_combine",
  contains = "cyl_copula",
  slots = c("sym.cop",
            "background.cop",
            "flip_up",
            "sym_rect")
)



#' Construction of \code{cyl_rect_combine} objects
#'
#' @param copula A \code{cyl_vonmises} or \code{Copula} object. The copula in
#'   the rectangles.
#' @param background A \code{cyl_copula} or \code{Copula} object. The copula
#'   where no rectangles overlay the unit square.
#' @param low_rect A numeric vector containing the lower and upper edge
#'   (u-value) of the lower rectangle.
#' @param up_rect A numeric vector containing the lower and upper edge (u-value)
#'   of the upper rectangle,  or the character string "symmetric" if it should be the mirror image
#'   (with respect to u=0.5) of the lower rectangle.
#' @param flip_up Logical value indicating whether the copula (\code{sym.cop})
#'   is rotated 90 degrees in the upper or lower rectangle
#'
#' @export
#'
#' @examples
#' cyl_rect_combine(copula = copula::frankCopula(param = 3), background =
#' cylcop::cyl_vonmises(), low_rect = c(0.1,0.4), up_rect = "symmetric", flip_up
#' = TRUE)
#'
cyl_rect_combine <-
  function(copula,
           background = indepCopula(),
           low_rect = c(0, 0.5),
           up_rect = "symmetric",
           flip_up = TRUE) {

    sym_rect = FALSE
# Checks and warnings

    if (!is.numeric(up_rect) && !identical(up_rect, "symmetric"))
      stop(
        cylcop::error_sound(),
        "up_rect must either be \"symmetric\", or provide upper and lower bound of the upper rectangle"
      )

    if (identical(up_rect, "symmetric")) {
      up_rect <- c(1 - low_rect[2], 1 - low_rect[1])
      sym_rect = TRUE
    }

    if (isTRUE(all.equal(low_rect[1], 0.0)) &&
        isTRUE(all.equal(up_rect[2], 1.0)) &&
        !isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      warning(cylcop::warning_sound(),
              "These rectangles do not lead to a periodic copula in u-direction!")
    }

    if (!isTRUE(all.equal(low_rect[1], 1 - up_rect[2])) ||
        !isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      warning(cylcop::warning_sound(), "These rectangles are not mirror images!")
    }


# Get description of Copula objects in the rectangles

    if (any(is(copula) == "cyl_vonmises")) {
      orig_copula <- copula
      name <- orig_copula@name
    }

    else if (any(is(copula) == "rotCopula")) {
      orig_copula <- copula@copula
      name <-
        class(orig_copula)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula (rotated)")
    }
    else if (any(is(copula) == "Copula") &&
             any(is(copula) != "rotCopula")) {
      orig_copula <- copula
      name <-
        class(copula)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula")
    }
    else
      stop(
        cylcop::error_sound(),
        "provide a (rotated) 'copula'-object from the 'copula'-package or a circular vonMises-object as input"
      )


# Get description of Copula objects not the rectangles (i.e. in the "background")

    if (any(is(background) == "cyl_copula")) {
      orig_copula_bg <- background
      name_bg <- orig_copula_bg@name
    }
    else if (any(is(background) == "rotCopula")) {
      orig_copula_bg <- background@copula
      name_bg <-
        class(orig_copula_bg)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula (rotated)")
      warning(cylcop::warning_sound(), "Background is not periodic in u-direction")
    }
    else if (any(is(background) == "Copula") &&
             any(is(background) != "rotCopula")) {
      orig_copula_bg <- background
      name_bg <-
        class(background)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula")
      if (!any(is(background) == "indepCopula") &&
          !any(is(copula) == "fhCopula")) {
        warning(cylcop::warning_sound(),
                "Background is not periodic in u-direction")
      }
    }
    else
      stop(
        cylcop::error_sound(),
        "provide a (rotated) 'copula'-object from the 'copula'-package or a 'cyl_copula'-object as background"
      )

    if (any(low_rect < c(0, 0)) ||
        any(up_rect > c(1, 1)))
      stop(cylcop::error_sound(), "rectangle boundaries not between 0 and 1")

# Echo when instantiating
    cat(name,if(!flip_up) "(rotated)", "put into a rectangle from \nu=",
      low_rect[1], "to", low_rect[2], "and\n",
      name,if(flip_up) "(rotated)", "put into a rectangle from \nu=",
      up_rect[1], "to", up_rect[2],
      "\noutside the rectangles, the copula is",
      name_bg, "\n"
    )

# Get parameters from copula objects
    parameters <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@parameters
      },
      if ("parameters" %in% methods::slotNames(orig_copula_bg)) {
        orig_copula_bg@parameters
      },
      low_rect, up_rect)
    param_names <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.names
      },
      if ("parameters" %in% methods::slotNames(orig_copula_bg)) {
        paste0 ("bg_", orig_copula_bg@param.names)
      },
      "low_rect1", "low_rect2", "up_rect1", "up_rect2")
    param_lowbnd <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.lowbnd
      },
      if ("parameters" %in% methods::slotNames(orig_copula_bg)) {
        orig_copula_bg@param.lowbnd
      },
      0, 0, 0.5, 0.5)
    param_upbnd <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.upbnd
      },
      if ("parameters" %in% methods::slotNames(orig_copula_bg)) {
        orig_copula_bg@param.upbnd
      },
      0.5, 0.5, 1, 1)

    new(
      "cyl_rect_combine",
      name = paste("symmetrized", name),
      parameters = parameters,     #contains all perameters of the copula in the rectangles
                                   #and of the copula outside the rectangles and the rectangles themselves
      param.names = param_names,   #dito
      param.lowbnd = param_lowbnd,
      param.upbnd = param_upbnd,
      sym.cop = copula,            # The copula in the rectangle
      background.cop = background, # The copula outside the rectangles
      flip_up = flip_up,
      sym_rect = sym_rect          # Force upper rectangle to be mirror image of lower
    )
  }



#' cyl_rect_combine show method
#' @param object A \code{cyl_rect_combine}-object
#' @describeIn cyl_copula-class What is printed
#' @export
setMethod("show", "cyl_rect_combine", function(object) {
  cat(object@name, "\n")
  for (i in seq_along(object@parameters)) {
    cat(object@param.names[i], "=", object@parameters[i], "\n")
  }
  if (object@flip_up)
    cat("upper copula is rotated 90 deg\n")
  else
    cat("lower copula is rotated 90 deg\n")
  cat("rectangles_symmetric:", object@sym_rect)
})



#' Generate random samples
#' @rdname Copula
#' @export
setMethod("rCopula", signature("numeric", "cyl_rect_combine"), function(n, copula) {
  sample <- matrix(ncol = 2,
                   nrow = n,
                   dimnames = list(c(), c("u", "v")))
  low_rect <-
    copula@parameters[match(c("low_rect1", "low_rect2"), copula@param.names)]
  up_rect <-
    copula@parameters[match(c("up_rect1", "up_rect2"), copula@param.names)]


# The case when the rectangles together cover the entire unit square is treated separatley
# becasue the equations are a lot simpler and faster to calculate (see text)

  if (isTRUE(all.equal(low_rect, c(0, 0.5))) &&
      isTRUE(all.equal(up_rect, c(0.5, 1))) &&
      any(is(copula@background.cop) == "indepCopula")) {
    for (i in 1:n) {
      #background copula
      c0 <- rCopula(1, indepCopula())
      # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
      c1 <- rCopula(1, copula@sym.cop)

      #rotate copula in the other rectangle
      if (any(is(copula@sym.cop) == "cyl_vonmises")) {
        rotated_cop <- copula@sym.cop
        rotated_cop@flip <- !copula@sym.cop@flip
        c2 <- rCopula(1, rotated_cop)
      }
      else if (any(is(copula@sym.cop) == "upfhCopula")) {
        c2 <- rCopula(1, fhCopula("lower"))
      }
      else if (any(is(copula@sym.cop) == "lowfhCopula")) {
        c2 <- rCopula(1, fhCopula("upper"))
      }
      else{
        c2 <- rCopula(1, rotCopula(copula@sym.cop, flip = c(TRUE, FALSE)))
      }

      #drawn value from background copula that is in lower rectangle
      if (c0[1] < 0.5) {
        if (copula@flip_up)
          transformed <- c(c1[1] / 2, c1[2])
        else
          transformed <- c(c2[1] / 2, c2[2])

      }
      #drawn value from background copula that is in upper rectangle
      else{
        if (copula@flip_up)
          transformed <- c((c2[1] + 1) / 2, c2[2])
        else
          transformed <- c((c1[1] + 1) / 2, c1[2])
      }
      sample[i, ] <- transformed
    }
    return(sample)
  }


# If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
# we need to do the full calculations

  else{
    for (i in 1:n) {
      lo1 <- c(low_rect[1], 0)
      lo2 <- c(low_rect[2], 1)
      up1 <- c(up_rect[1], 0)
      up2 <- c(up_rect[2], 1)

      #background copula
      c0 <- rCopula(1, copula@background.cop)

      # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
      c1 <- rCopula(1, copula@sym.cop)

       #rotate copula in the other rectangle
      if (any(is(copula) == "cyl_vonmises")) {
        rotated_cop <- copula@sym.cop
        rotated_cop@flip <- !copula@sym.cop@flip
        c2 <- rCopula(1, rotated_cop)
      }
      else if (any(is(copula@sym.cop) == "upfhCopula")) {
        c2 <- rCopula(1, fhCopula("lower"))
      }
      else if (any(is(copula@sym.cop) == "lowfhCopula")) {
        c2 <- rCopula(1, fhCopula("upper"))
      }
      else{
        c2 <- rCopula(1, rotCopula(copula@sym.cop, flip = c(TRUE, FALSE)))
      }

      #C0-volume of the lower rectangle
      Vlo <- prob(copula@background.cop, lo1, lo2)
      #C0-volume of the upper rectangle
      Vup <- prob(copula@background.cop, up1, up2)

      #See text for expalantion of psi
      psi1_rect_lo <-
        GoFKernel::inverse(
          f = function(x) {prob(copula@background.cop, lo1, c(x, 1)) / Vlo},
          lower = low_rect[1],
          upper = low_rect[2]
        )

      psi2_rect_lo <-
        GoFKernel::inverse(
          f = function(y) {prob(copula@background.cop, lo1, c(low_rect[2], y)) / Vlo},
          lower = 0,
          upper = 1
        )

      psi1_rect_up <-
        GoFKernel::inverse(
          f = function(x) {prob(copula@background.cop, up1, c(x, 1)) / Vup},
          lower = up_rect[1],
          upper = up_rect[2]
        )

      psi2_rect_up <-
        GoFKernel::inverse(
          f = function(y) {prob(copula@background.cop, up1, c(up_rect[2], y)) / Vup},
          lower = 0,
          upper = 1
        )

# Draw from background copula in lower rectangle
      if (c0[1] < low_rect[2] &&
          c0[1] > low_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_lo(c1[1]), psi2_rect_lo(c1[2]))
        else
          transformed <- c(psi1_rect_lo(c2[1]), psi2_rect_lo(c2[2]))
      }

# Draw from background copula in upper rectangle
      else if (c0[1] < up_rect[2] &&
               c0[1] > up_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_up(c2[1]), psi2_rect_up(c2[2]))
        else
          transformed <- c(psi1_rect_up(c1[1]), psi2_rect_up(c1[2]))
      }

# Draw from background copula in neither rectangle
      else{
        transformed <- c0
      }
      sample[i, ] <- transformed
    }
    return(sample)
  }
})



#' Calcualte density
#' @rdname Copula
#' @export
setMethod("dCopula", signature("matrix", "cyl_rect_combine"), function(u, copula) {
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  low_rect <-
    copula@parameters[match(c("low_rect1", "low_rect2"), copula@param.names)]
  up_rect <-
    copula@parameters[match(c("up_rect1", "up_rect2"), copula@param.names)]


# The case when the rectangles together cover the entire unit square is treated separatley
# because the equations are a lot simpler and faster to calculate (see text)

  if (isTRUE(all.equal(low_rect, c(0, 0.5))) &&
      isTRUE(all.equal(up_rect, c(0.5, 1))) &&
      any(is(copula@background.cop) == "indepCopula")) {

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1 <- copula@sym.cop

    #rotate copula in the other rectangle
    if (any(is(copula) == "cyl_vonmises")) {
      rotated_cop <- copula@sym.cop
      rotated_cop@flip <- !copula@sym.cop@flip
      c2 <- rotated_cop
    }
    else if (any(is(copula@sym.cop) == "upfhCopula")) {
      c2 <- fhCopula("lower")
    }
    else if (any(is(copula@sym.cop) == "lowfhCopula")) {
      c2 <- fhCopula("upper")
    }
    else{
      c2 <- rotCopula(copula@sym.cop, flip = c(TRUE, FALSE))
    }

    pdf <- map2_dbl(u, v, function(u, v) {
      #(u,v) is in lower rectangle
      if (u < 0.5) {
        if (copula@flip_up)
          0.5 * dCopula(c(2 * u, v), c1)
        else
          0.5 * dCopula(c(2 * u, v), c2)
      }

      #(u,v) is in upper rectangle
      else{
        if (copula@flip_up)
          0.5 * dCopula(c(2 * u - 1, v), c2)
        else
          0.5 * dCopula(c(2 * u - 1, v), c1)
      }
    })
    return(pdf)
  }


# If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
# we need to do the full calculations

  else{
    lo1 <- c(low_rect[1], 0)
    lo2 <- c(low_rect[2], 1)
    up1 <- c(up_rect[1], 0)
    up2 <- c(up_rect[2], 1)

    #background copula
    c0 <- copula@background.cop

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1 <- copula@sym.cop

    #rotate copula in the other rectangle
    if (any(is(copula) == "cyl_vonmises")) {
      rotated_cop <- copula@sym.cop
      rotated_cop@flip <- !copula@sym.cop@flip
      c2 <- rotated_cop
    }
    else if (any(is(copula@sym.cop) == "upfhCopula")) {
      c2 <- fhCopula("lower")
    }
    else if (any(is(copula@sym.cop) == "lowfhCopula")) {
      c2 <- fhCopula("upper")
    }
    else{
      c2 <- rotCopula(copula@sym.cop, flip = c(TRUE, FALSE))
    }

    #C0-volume of the lower rectangle
    Vlo <- prob(copula@background.cop, lo1, lo2)
    #C0-volume of the upper rectangle
    Vup <- prob(copula@background.cop, up1, up2)

    pdf <- map2_dbl(u, v, function(u, v) {
#(u,v) is in lower rectangle
      if (u < lo2[1] && u > lo1[1]) {
        Vlou <- prob(copula@background.cop, lo1, c(u, 1))
        Vlov <-
          prob(copula@background.cop, lo1, c(low_rect[2], v))
        integrand <- function(x) {
          dCopula(c(x[1], v), copula@background.cop)
        }
        condv <-
          stats::integrate(f = Vectorize(integrand),
                    lower = lo1[1],
                    upper = lo2[1])
        if (copula@flip_up)
          pdf <-
          (condv$value / Vlo) * dCopula(c(Vlou / Vlo, Vlov / Vlo), c1)
        else
          pdf <-
          (condv$value / Vlo) * dCopula(c(Vlou / Vlo, Vlov / Vlo), c2)
      }

#(u,v) is in upper rectangle
      else if (u < up2[1] && u > up1[1]) {
        Vupu <- prob(copula@background.cop, up1, c(u, 1))
        Vupv <-
          prob(copula@background.cop, up1, c(up_rect[2], v))
        integrand <- function(x) {
          dCopula(c(x[1], v), copula@background.cop)
        }
        condv <-
          stats::integrate(f = Vectorize(integrand),
                    lower = up1[1],
                    upper = up2[1])
        if (copula@flip_up)
          pdf <-
          (condv$value / Vup) * dCopula(c(Vupu / Vup, Vupv / Vup), c2)
        else
          pdf <-
          (condv$value / Vup) * dCopula(c(Vupu / Vup, Vupv / Vup), c1)
      }

      #(u,v) is in neither rectangle
      else {
        pdf <- dCopula(c(u, v), c0)
      }
      return(pdf)
    })
    return(pdf)
  }
})



#' Calcualte distribution
#' @rdname Copula
#' @export
setMethod("pCopula", signature("matrix", "cyl_rect_combine"), function(u, copula) {
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  low_rect <-
    copula@parameters[match(c("low_rect1", "low_rect2"), copula@param.names)]
  up_rect <-
    copula@parameters[match(c("up_rect1", "up_rect2"), copula@param.names)]


# The case when the rectangles together cover the entire unit square is treated separatley
# because the equations are a lot simpler and faster to calculate (see text)

  if (isTRUE(all.equal(low_rect, c(0, 0.5))) &&
      isTRUE(all.equal(up_rect, c(0.5, 1))) &&
      any(is(copula@background.cop) == "indepCopula")) {

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1 <- copula@sym.cop

    #rotate copula in the other rectangle
    if (any(is(copula) == "cyl_vonmises")) {
      rotated_cop <- copula@sym.cop
      rotated_cop@flip <- !copula@sym.cop@flip
      c2 <- rotated_cop
    }
    else if (any(is(copula@sym.cop) == "upfhCopula")) {
      c2 <- fhCopula("lower")
    }
    else if (any(is(copula@sym.cop) == "lowfhCopula")) {
      c2 <- fhCopula("upper")
    }
    else{
      c2 <- rotCopula(copula@sym.cop, flip = c(TRUE, FALSE))
    }

    cdf <- map2_dbl(u, v, function(u, v) {
#(u,v) is in lower rectangle
      if (u < 0.5) {
        if (copula@flip_up)
          0.5 * pCopula(c(2 * u, v), c1)
        else
          0.5 * pCopula(c(2 * u, v), c2)
      }

#(u,v) is in upper rectangle
      else{
        if (copula@flip_up)
          0.5 * pCopula(c(2 * u - 1, v), c2) + 0.5 * v
        else
          0.5 * pCopula(c(2 * u - 1, v), c1) + 0.5 * v
      }
    })
    return(cdf)
  }


# If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
# we need to do the full calculations

  else{
    lo1 <- c(low_rect[1], 0)
    lo2 <- c(low_rect[2], 1)
    up1 <- c(up_rect[1], 0)
    up2 <- c(up_rect[2], 1)

    #background copula
    c0 <- copula@background.cop

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1 <- copula@sym.cop

    #rotate copula in the other rectangle
    if (any(is(copula) == "cyl_vonmises")) {
      rotated_cop <- copula@sym.cop
      rotated_cop@flip <- !copula@sym.cop@flip
      c2 <- rotated_cop
    }
    else if (any(is(copula@sym.cop) == "upfhCopula")) {
      c2 <- fhCopula("lower")
    }
    else if (any(is(copula@sym.cop) == "lowfhCopula")) {
      c2 <- fhCopula("upper")
    }
    else{
      c2 <- rotCopula(copula@sym.cop, flip = c(TRUE, FALSE))
    }

    #C0-volume of the lower rectangle
    Vlo <- prob(copula@background.cop, lo1, lo2)
    #C0-volume of the upper rectangle
    Vup <- prob(copula@background.cop, up1, up2)

    cdf <- map2_dbl(u, v, function(u, v) {
 #(u,v) is in lower rectangle
      if (u < low_rect[2] && u > low_rect[1]) {
        Vlou <- prob(copula@background.cop, lo1, c(u, 1))
        Vlov <-
          prob(copula@background.cop, lo1, c(low_rect[2], v))
        marg1 <-
          pCopula(c(u, 0), c0) + pCopula(c(low_rect[1], v), c0) - pCopula(c(low_rect[1], 0), c0)
        if (copula@flip_up)
          cdf <- Vlo * pCopula(c(Vlou / Vlo, Vlov / Vlo), c1) + marg1
        else
          cdf <- Vlo * pCopula(c(Vlou / Vlo, Vlov / Vlo), c2) + marg1
      }

#(u,v) is in upper rectangle
      else if (u < up_rect[2] && u > up_rect[1]) {
        Vupu <- prob(copula@background.cop, up1, c(u, 1))
        Vupv <-
          prob(copula@background.cop, up1, c(up_rect[2], v))
        marg2 <-
          pCopula(c(u, 0), c0) + pCopula(c(up_rect[1], v), c0) - pCopula(c(up_rect[1], 0), c0)
        if (copula@flip_up)
          cdf <- Vup * pCopula(c(Vupu / Vup, Vupv / Vup), c2) + marg2
        else
          cdf <- Vup * pCopula(c(Vupu / Vup, Vupv / Vup), c1) + marg2
      }

#(u,v) is in neither rectangle
      else {
        cdf <- pCopula(c(u, v), c0)
      }
      return(cdf)
    })
    return(cdf)
  }
})



#-----Change attributes of existing cyl_rect_combine object.-------------------------------------------
#
#' @rdname setCopParam
#' @export
setMethod("setCopParam", "cyl_rect_combine", function(copula, param_val, param_name) {
  if (is.null(param_name))
    param_name <- copula@param.names
  param_num <- param_num_checked(copula, param_val, param_name)
  copula@parameters[param_num] <- param_val

#make sure the rectangles are still symemtrical
  if (copula@sym_rect) {
    if (any(stringr::str_starts(param_name, "up_")))
      stop(
        cylcop::error_sound(),
        "the copula has symmetric rectangles defined by the bounds of the lower one.",
        "\nChanging up_rect has no effect"
      )
    copula@parameters[match("up_rect1", copula@param.names)] <-
      1 - copula@parameters[match("low_rect2", copula@param.names)]
    copula@parameters[match("up_rect2", copula@param.names)] <-
      1 - copula@parameters[match("low_rect1", copula@param.names)]
  }


#set copula parameters

  if (any(stringr::str_starts(param_name, "bg_|up_|low_", negate = T))) {
    cop_param_val <-
      param_val[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = T))]

    if (!any(is(copula@sym.cop) == "rotCopula")) {
      cop_param_num <-
        param_name[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = T))] %>%
        match(copula@sym.cop@param.names)
      copula@sym.cop@parameters[cop_param_num] <- cop_param_val
    }
    else{
      cop_param_num <-
        param_name[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = T))] %>%
        match(copula@sym.cop@copula@param.names)
      copula@sym.cop@copula@parameters[cop_param_num] <-
        cop_param_val
    }
  }


#set background copula parameters

  if (any(stringr::str_starts(param_name, "bg_"))) {
    bg_param_val <- param_val[which(stringr::str_starts(param_name, "bg_"))]
    if (!any(is(copula@background.cop) == "rotCopula")) {
      bg_param_num <- param_name[which(stringr::str_starts(param_name, "bg_"))] %>%
        stringr::str_remove("bg_") %>%
        match(copula@background.cop@param.names)
      copula@background.cop@parameters[bg_param_num] <- bg_param_val
    }
    else{
      bg_param_num <- param_name[which(stringr::str_starts(param_name, "bg_"))] %>%
        stringr::str_remove("bg_") %>%
        match(copula@background.cop@copula@param.names)
      copula@background.cop@copula@parameters[bg_param_num] <-
        bg_param_val
    }
  }

  return(copula)
})
