#' @include cyl_cop_class.R
NULL



#' An S4 Class of Circular-Linear Copulas Generated from a Rectangular Patchwork
#'
#' This class contains bivariate circular-linear copulas generated from
#' linear-linear bivariate '\code{\linkS4class{Copula}}' objects of the package
#' '\pkg{copula}' or circular-linear copulas of class '\code{\linkS4class{cyl_copula}}'.
#' 2 non-overlapping rectangles are laid over the unit square, both have width
#' 1 in v-direction. In the area covered by the first rectangle, the copula is
#' derived from a linear-linear bivariate '\code{\linkS4class{Copula}}' object.
#' Rectangle 2 contains the same copula as rectangle 1, but 90 degrees rotated.
#' In the area not covered by the rectangles, the "background", the copula is
#' derived from a circular-linear '\code{\linkS4class{cyl_copula}}' object.
#' The copula regions are combined in a way that the overall result on the entire
#' unit square is also a copula.
#'
#' With appropriate choices of the rectangles this results in copulas
#' that are periodic in u-direction (and not in v-direction) and therefore are
#' circular-linear. When the 2 rectangles are mirror images with
#' respect to \eqn{u=0.5}, the resulting overall copula is symmetric with respect
#' to \eqn{u=0.5}, i.e. there is symmetry between positive and negative angles.
#'
#' Note that as "background copula", we can also chose a linear-linear copula,
#' the overall result will then, however, not be a symmetric circular linear copula.
#'
#' @section Objects from the Class: Objects are created by
#'   \code{\link{cyl_rect_combine}()}.
#'
#' @slot name \link[base]{character} string holding the name of the copula.
#' @slot parameters \link[base]{numeric} \link[base]{vector} holding the parameter values.
#' @slot param.names \link[base]{character} \link[base]{vector} the parameter names.
#' @slot param.lowbnd \link[base]{numeric} \link[base]{vector} holding the lower bounds of the
#'   parameters.
#' @slot param.upbnd \link[base]{numeric} \link[base]{vector} holding the upper bounds of the
#'   parameters.
#' @slot sym.cop '\code{\linkS4class{Copula}}' object of the package
#' '\pkg{copula}' or '\code{\linkS4class{cyl_vonmises}}' object. The copula in
#'   the rectangles.
#' @slot background.cop '\code{\linkS4class{cyl_vonmises}}' or
#' '\code{\linkS4class{Copula}}' object of the package '\pkg{copula}'
#' (does not lead to an overall symmetric circular-linear copula). The copula
#'   where no rectangles overlay the unit square.
#' @slot flip_up \link[base]{logical} value indicating whether the copula (\code{sym.cop}) is
#'   rotated 90 degrees in the upper or lower rectangle.
#' @slot sym_rect \link[base]{logical} value indicating whether the upper rectangle was
#'   forced to be a mirror image of the lower one with respect to u=0.5 at the
#'   construction of the object.
#'
#' @section Extends:
#' Class '\code{cyl_rect_combine}' extends class '\code{\linkS4class{Copula}}'.
#'
#' @references \insertRef{Durante2009}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
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



#' Construction of '\code{cyl_rect_combine}' Objects
#'
#' Constructs a circular-linear copula of class
#' '\code{\linkS4class{cyl_rect_combine}}' from a rectangular patchwork of copulas.
#'
#'
#' @param copula '\code{\linkS4class{Copula}}' object of the package
#' '\pkg{copula}' or '\code{\linkS4class{cyl_vonmises}}' object, the copula in
#'   the rectangles.
#' @param background '\code{\linkS4class{cyl_vonmises}}' or
#' '\code{\linkS4class{Copula}}' object of the package '\pkg{copula}'
#' (does not lead to an overall symmetric circular-linear copula), the copula
#'   where no rectangles overlay the unit square.
#' @param low_rect \link[base]{numeric} \link[base]{vector} containing the
#' lower and upper edge (u-value) of the lower rectangle.
#' @param up_rect \link[base]{numeric} \link[base]{vector} containing the lower
#' and upper edge (u-value) of the upper rectangle, or the character string
#' "symmetric" if it should be the mirror image (with respect to u=0.5) of the lower rectangle.
#' @param flip_up \link[base]{logical} value indicating whether the copula (\code{sym.cop})
#'   is rotated 90 degrees in the upper (\code{flip_up = TRUE}) or lower rectangle.
#'
#' @export
#'
#' @examples
#'
#' #symmetric rectangles spanning entire unit square
#' cop <- cyl_rect_combine(copula::frankCopula(2))
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' #symmetric rectangles, independence copula as background
#' cop <- cyl_rect_combine(copula::frankCopula(2),
#'   low_rect = c(0, 0.3),
#'   up_rect = "symmetric",
#'   flip_up = FALSE
#' )
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' #symmetric rectangles, cy_quadsec-copula as background
#' cop <- cyl_rect_combine(copula::normalCopula(0.3),
#'   low_rect = c(0.1, 0.4),
#'   up_rect = "symmetric",
#'   background = cyl_quadsec(-0.1)
#' )
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' #asymmetric rectangles, von Mises copula as background.
#' #!!Not a symmetric circular linear copula!!
#' cop <- cyl_rect_combine(copula::normalCopula(0.3), low_rect = c(0.1, 0.4),
#' up_rect = c(0.5, 0.7), background = cyl_vonmises(mu = pi, kappa = 0.3))
#' cop_plot(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#'
#' @references \insertRef{Durante2009}{cylcop}
#'
#' \insertRef{Hodelappl}{cylcop}
#'
#' \insertRef{Hodelmethod}{cylcop}
#'
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
        error_sound(),
        "up_rect must either be \"symmetric\", or provide upper and lower bound of the upper rectangle"
      )

    if (identical(up_rect, "symmetric")) {
      up_rect <- c(1 - low_rect[2], 1 - low_rect[1])
      sym_rect = TRUE
    }

    if (isTRUE(all.equal(low_rect[1], 0.0)) &&
        isTRUE(all.equal(up_rect[2], 1.0)) &&
        !isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      warning(warning_sound(),
              "These rectangles do not lead to a periodic copula in u-direction!")
    }

    if (!isTRUE(all.equal(low_rect[1], 1 - up_rect[2])) ||
        !isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      warning(warning_sound(), "These rectangles are not mirror images!")
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
        error_sound(),
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
      warning(warning_sound(), "Background is not periodic in u-direction")
    }
    else if (any(is(background) == "Copula") &&
             any(is(background) != "rotCopula")) {
      orig_copula_bg <- background
      name_bg <-
        class(background)[1] %>% stringr::str_remove("Copula") %>% stringr::str_to_sentence(locale = "en") %>% paste("copula")
      if (!any(is(background) == "indepCopula") &&
          !any(is(copula) == "fhCopula")) {
        warning(warning_sound(),
                "Background is not periodic in u-direction")
      }
    }
    else
      stop(
        error_sound(),
        "provide a (rotated) 'copula'-object from the 'copula'-package or a 'cyl_copula'-object as background"
      )

    if (any(low_rect < c(0, 0)) ||
        any(up_rect > c(1, 1)))
      stop(error_sound(), "rectangle boundaries not between 0 and 1")

# Echo when instantiating

    if(cylcop.env$silent==F){
      message(name,if(!flip_up) " (rotated)", "put into a rectangle from \nu= ",
      low_rect[1], " to ", low_rect[2], " and\n",
      name,if(flip_up) " (rotated) ", " put into a rectangle from \nu= ",
      up_rect[1], " to ", up_rect[2],
      "\noutside the rectangles, the copula is ",
      name_bg, "\n"
    )}

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
      name = paste("Rectangular patchwork of", name),
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
#' @rdname show-cyl_copula-method
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
#' @rdname Cylcop
# @describeIn cyl_rect_combine-class Generate random samples.
#' @export
setMethod("rcylcop", signature("numeric", "cyl_rect_combine"), function(n, copula) {
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

    #background copula
    c0<-rcylcop(n, indepCopula())

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1<-rcylcop(n, copula@sym.cop)

    #rotate copula in the other rectangle
    if (any(is(copula@sym.cop) == "cyl_vonmises")) {
      rotated_cop <- copula@sym.cop
      rotated_cop@flip <- !copula@sym.cop@flip
      c2 <- rcylcop(n, rotated_cop)
    }
    else if (any(is(copula@sym.cop) == "upfhCopula")) {
      c2 <- rcylcop(n, fhCopula("lower"))
    }
    else if (any(is(copula@sym.cop) == "lowfhCopula")) {
      c2 <- rcylcop(n, fhCopula("upper"))
    }
    else{
      c2 <- rcylcop(n, rotCopula(copula@sym.cop, flip = c(TRUE, FALSE)))
    }

    for (i in 1:n) {
      #drawn value from background copula that is in lower rectangle
      if (c0[i,1] < 0.5) {
        if (copula@flip_up)
          transformed <- c(c1[i,1] / 2, c1[i,2])
        else
          transformed <- c(c2[i,1] / 2, c2[i,2])

      }
      #drawn value from background copula that is in upper rectangle
      else{
        if (copula@flip_up)
          transformed <- c((c2[i,1] + 1) / 2, c2[i,2])
        else
          transformed <- c((c1[i,1] + 1) / 2, c1[i,2])
      }
      sample[i, ] <- transformed
    }
    return(sample)
  }


# If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
# we need to do the full calculations
else{
  lo1 <- c(low_rect[1], 0)
  lo2 <- c(low_rect[2], 1)
  up1 <- c(up_rect[1], 0)
  up2 <- c(up_rect[2], 1)

  #background copula
  c0 <- rcylcop(n, copula@background.cop)

  # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
  c1 <- rcylcop(n, copula@sym.cop)

  #rotate copula in the other rectangle
  if (any(is(copula) == "cyl_vonmises")) {
    rotated_cop <- copula@sym.cop
    rotated_cop@flip <- !copula@sym.cop@flip
    c2 <- rcylcop(n, rotated_cop)
  }
  else if (any(is(copula@sym.cop) == "upfhCopula")) {
    c2 <- rcylcop(n, fhCopula("lower"))
  }
  else if (any(is(copula@sym.cop) == "lowfhCopula")) {
    c2 <- rcylcop(n, fhCopula("upper"))
  }
  else{
    c2 <- rcylcop(n, rotCopula(copula@sym.cop, flip = c(TRUE, FALSE)))
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

    for (i in 1:n) {
# Draw from background copula in lower rectangle
      if (c0[i,1] < low_rect[2] &&
          c0[i,1] > low_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_lo(c1[i,1]), psi2_rect_lo(c1[i,2]))
        else
          transformed <- c(psi1_rect_lo(c2[i,1]), psi2_rect_lo(c2[i,2]))
      }

# Draw from background copula in upper rectangle
      else if (c0[i,1] < up_rect[2] &&
               c0[i,1] > up_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_up(c2[i,1]), psi2_rect_up(c2[i,2]))
        else
          transformed <- c(psi1_rect_up(c1[i,1]), psi2_rect_up(c1[i,2]))
      }

# Draw from background copula in neither rectangle
      else{
        transformed <- c0[i,]
      }
      sample[i, ] <- transformed
    }
    return(sample)
  }
})



#' Calcualte density
#' @rdname Cylcop
# @describeIn cyl_rect_combine-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_rect_combine"), function(u, copula) {
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
    if (any(is(copula@sym.cop) == "cyl_vonmises")) {
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
          dcylcop(c(2 * u, v), c1)
        else
          dcylcop(c(2 * u, v), c2)
      }

      #(u,v) is in upper rectangle
      else if (u > 0.5){
        if (copula@flip_up)
          dcylcop(c(2 * u - 1, v), c2)
        else
          dcylcop(c(2 * u - 1, v), c1)
      }
      #(u,v) is on the border, i.e. density of indepCopula
      else{
        1
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
      if (u <= lo2[1] && u >= lo1[1]) {
        Vlou <- prob(copula@background.cop, lo1, c(u, 1))
        Vlov <-
          prob(copula@background.cop, lo1, c(low_rect[2], v))

        #calculate integral as difference of conditional copulae to avoid numerical integration
        integ <- cylcop::ccylcop(c(lo2[1],v), copula@background.cop, cond_on=2, inverse=F)-
          cylcop::ccylcop(c(lo1[1],v), copula@background.cop, cond_on=2, inverse=F)
        if (copula@flip_up)
          pdf <-
          (integ / Vlo) * dcylcop(c(Vlou / Vlo, Vlov / Vlo), c1)
        else
          pdf <-
          (integ / Vlo) * dcylcop(c(Vlou / Vlo, Vlov / Vlo), c2)
      }

#(u,v) is in upper rectangle
      else if (u <= up2[1] && u >= up1[1]) {
        Vupu <- prob(copula@background.cop, up1, c(u, 1))
        Vupv <-
          prob(copula@background.cop, up1, c(up_rect[2], v))

        integ <- cylcop::ccylcop(c(up2[1],v), copula@background.cop, cond_on=2, inverse=F)-
          cylcop::ccylcop(c(up1[1],v), copula@background.cop, cond_on=2, inverse=F)


        if (copula@flip_up)
          pdf <-
          (integ / Vup) * dcylcop(c(Vupu / Vup, Vupv / Vup), c2)
        else
          pdf <-
          (integ / Vup) * dcylcop(c(Vupu / Vup, Vupv / Vup), c1)
      }

      #(u,v) is in neither rectangle
      else {
        pdf <- dcylcop(c(u, v), c0)
      }
      return(pdf)
    })
    return(pdf)
  }
})



#' Calcualte distribution
#' @rdname Cylcop
# @describeIn cyl_rect_combine-class Calculate the distribution.
#' @export
setMethod("pcylcop", signature("matrix", "cyl_rect_combine"), function(u, copula) {
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
          0.5 * pcylcop(c(2 * u, v), c1)
        else
          0.5 * pcylcop(c(2 * u, v), c2)
      }

#(u,v) is in upper rectangle
      else if (u > 0.5) {
        if (copula@flip_up)
          0.5 * pcylcop(c(2 * u - 1, v), c2) + 0.5 * v
        else
          0.5 * pcylcop(c(2 * u - 1, v), c1) + 0.5 * v
      }
#(u,v) is on border, indepCopula
      else{
        u*v
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
          pcylcop(c(u, 0), c0) + pcylcop(c(low_rect[1], v), c0) - pcylcop(c(low_rect[1], 0), c0)
        if (copula@flip_up)
          cdf <- Vlo * pcylcop(c(Vlou / Vlo, Vlov / Vlo), c1) + marg1
        else
          cdf <- Vlo * pcylcop(c(Vlou / Vlo, Vlov / Vlo), c2) + marg1
      }

#(u,v) is in upper rectangle
      else if (u < up_rect[2] && u > up_rect[1]) {
        Vupu <- prob(copula@background.cop, up1, c(u, 1))
        Vupv <-
          prob(copula@background.cop, up1, c(up_rect[2], v))
        marg2 <-
          pcylcop(c(u, 0), c0) + pcylcop(c(up_rect[1], v), c0) - pcylcop(c(up_rect[1], 0), c0)
        if (copula@flip_up)
          cdf <- Vup * pcylcop(c(Vupu / Vup, Vupv / Vup), c2) + marg2
        else
          cdf <- Vup * pcylcop(c(Vupu / Vup, Vupv / Vup), c1) + marg2
      }

#(u,v) is in neither rectangle
      else {
        cdf <- pcylcop(c(u, v), c0)
      }
      return(cdf)
    })
    return(cdf)
  }
})



#' Condtional copula
#' @rdname ccylcop
# @describeIn cyl_rect_combine-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_rect_combine"), function(u, copula, cond_on=2, inverse=F) {

  if(cond_on==2){
      if(inverse==F){
        numerical_conditional_cop(u,copula,cond_on = 2)
      }
      else{
        numerical_inv_conditional_cop(u,copula,cond_on = 2)
      }
  }
  else if(cond_on==1){
      if(inverse==F){
        numerical_conditional_cop(u,copula,cond_on = 1)
      }
      else{
        numerical_inv_conditional_cop(u,copula,cond_on = 1)
      }
  }
  else stop("cond_on must be either 1 or 2")
})



#-----Change attributes of existing cyl_rect_combine object.-------------------------------------------
#
#' @rdname setCopParam
#@describeIn cyl_rect_combine-class Change attributes of existing object.
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
        error_sound(),
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
