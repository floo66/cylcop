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
#' '\code{\linkS4class{Copula}}' object of the package '\pkg{copula}',
#' the copula where no rectangles overlay the unit square. If this copula is not
#' symmetric, the overall \code{cyl_rect_combine}-copula will also not be symmetric.
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
#' @param background '\code{\linkS4class{cyl_copula}}' or
#' '\code{\linkS4class{Copula}}' object of the package '\pkg{copula}',
#' the copula where no rectangles overlay the unit square. If this copula is not
#' symmetric, the overall \code{cyl_rect_combine}-copula will also not be symmetric.
#' @param low_rect \link[base]{numeric} \link[base]{vector} of length 2 containing the
#' lower and upper edge (u-value) of the lower rectangle.
#' @param up_rect \link[base]{numeric} \link[base]{vector} of length 2 containing the lower
#' and upper edge (u-value) of the upper rectangle, or the character string
#' "symmetric" if it should be the mirror image (with respect to u=0.5) of the lower rectangle.
#' @param flip_up \link[base]{logical} value indicating whether the copula (\code{copula})
#'   is rotated 90 degrees in the upper (\code{flip_up = TRUE}) or lower rectangle.
#'
#' @return An \R object of class '\code{\linkS4class{cyl_rect_combine}}'.
#'
#' @export
#'
#' @examples
#'
#' #symmetric rectangles spanning entire unit square
#' cop <- cyl_rect_combine(copula::frankCopula(2))
#' if(interactive()){
#'  plot_cop_surf(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#' }
#'
#' #symmetric rectangles, independence copula as background
#' cop <- cyl_rect_combine(copula::frankCopula(2),
#'   low_rect = c(0, 0.3),
#'   up_rect = "symmetric",
#'   flip_up = FALSE
#' )
#' if(interactive()){
#'  plot_cop_surf(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#' }
#'
#' #symmetric rectangles, cy_quadsec-copula as background
#' cop <- cyl_rect_combine(copula::normalCopula(0.3),
#'   low_rect = c(0.1, 0.4),
#'   up_rect = "symmetric",
#'   background = cyl_quadsec(-0.1)
#' )
#' if(interactive()){
#'  plot_cop_surf(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#' }
#'
#' #asymmetric rectangles, von Mises copula as background.
#' #!!Not a symmetric circular linear copula!!
#' cop <- cyl_rect_combine(copula::normalCopula(0.3),
#'   low_rect = c(0.1, 0.4),
#'   up_rect = c(0.5, 0.7),
#'   background = cyl_vonmises(mu = pi, kappa = 0.3)
#' )
#' if(interactive()){
#'  plot_cop_surf(copula = cop, type = "pdf", plot_type = "ggplot", resolution = 20)
#' }
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

    #validate input
    tryCatch({
      check_arg_all(list(check_argument_type(copula,
                                        type="Copula"),
                         check_argument_type(copula,
                                             type="cyl_copula")
                         )
                    ,2)
      check_arg_all(list(check_argument_type(background,
                                             type="Copula"),
                         check_argument_type(background,
                                             type="cyl_copula")
      )
      ,2)
      check_arg_all(check_argument_type(low_rect,
                                        type="numeric",
                                        length = 2,
                                        lower=0,
                                        upper=1)
                    ,1)
      check_arg_all(list(check_argument_type(up_rect,
                                             type="numeric",
                                             length = 2,
                                             lower=0,
                                             upper=1),
                         check_argument_type(up_rect,
                                             type="character",
                                             values = "symmetric")
      )
      ,2)
      check_arg_all(check_argument_type(flip_up,
                                             type="logical")
      ,1)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )

    sym_rect <-  FALSE
# Checks and warnings

    if (identical(up_rect, "symmetric")) {
      if(low_rect[2]>0.5){
        stop(
          error_sound(),
          "The rectangles are overlapping."
        )
      }
      up_rect <- c(1 - low_rect[2], 1 - low_rect[1])
      sym_rect <-  TRUE
    }

    if(!sym_rect){
      if(low_rect[2]>up_rect[1]){
        stop(
          error_sound(),
          "Rectangles must not be overlapping!"
        )
      }}

    if(low_rect[2]<low_rect[1]||up_rect[2]<up_rect[1]){
      stop(
        error_sound(),
        "Upper edge of a rectangle is below its lower edge."
      )
    }

    if (isTRUE(all.equal(low_rect[1], 0.0)) &&
        isTRUE(all.equal(up_rect[2], 1.0)) &&
        !isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      warning(warning_sound(),
              "These rectangles do not lead to a periodic copula in u-direction.")
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
        "provide a (rotated) 'copula'-object from the 'copula'-package or a circular von Mises-object as input"
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
      0, 0, 0, 0)
    param_upbnd <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.upbnd
      },
      if ("parameters" %in% methods::slotNames(orig_copula_bg)) {
        orig_copula_bg@param.upbnd
      },
      1, 1, 1, 1)

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
    c0<-runif(n,0,1)

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle

    if (!copula@flip_up){
      #rotate copula in the other rectangle
      if (any(is(copula@sym.cop) == "cyl_vonmises")) {
        sym_cop <- copula@sym.cop
        sym_cop@flip <- !copula@sym.cop@flip
      }else if (any(is(copula@sym.cop) == "upfhCopula")) {
        sym_cop <- fhCopula("lower")
      }else if (any(is(copula@sym.cop) == "lowfhCopula")) {
        sym_cop <- fhCopula("upper")
      }else{
        sym_cop <- rotCopula(copula@sym.cop, flip = c(TRUE, FALSE))
      }
    }else{
      sym_cop <- copula@sym.cop
    }

    lin_sample<-rcylcop(n, sym_cop)
    sample[,2] <- lin_sample[,2]

    ind <- which(c0 < 0.5)

    sample[ind,1] <- lin_sample[ind,1]/2

    ind <- which(c0 >= 0.5)

    sample[ind,1] <- (1-lin_sample[ind,1]/2)


    return(sample)
  }


  # If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
  # we need to do the full calculations
  else{
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
    Vlo <- low_rect[2]-low_rect[1]
    #C0-volume of the upper rectangle
    Vup <- up_rect[2]-up_rect[1]

    #See text for expalantion of psi
    psi1_rect_lo_test <- function(x){
      x*Vlo+low_rect[1]
    }
    psi2_rect_lo_test <-
      GoFKernel::inverse(
        f = function(y) {(pcylcop(c(low_rect[2], y), copula@background.cop) - pcylcop(c(low_rect[1], y), copula@background.cop)) / Vlo},
        lower = 0,
        upper = 1
      )
    psi1_rect_up_test <- function(x){
      x * Vup + up_rect[1]
    }
    psi2_rect_up_test <-
      GoFKernel::inverse(
        f = function(y) {(pcylcop(c(up_rect[2], y), copula@background.cop) - pcylcop(c(up_rect[1], y), copula@background.cop)) / Vup},
        lower = 0,
        upper = 1
      )

    for (i in seq_len(n)) {
      # Draw from background copula fell in lower rectangle
      if (c0[i,1] < low_rect[2] &&
          c0[i,1] > low_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_lo_test(c1[i,1]), psi2_rect_lo_test(c1[i,2]))
        else
          transformed <- c(psi1_rect_lo_test(c2[i,1]), psi2_rect_lo_test(c2[i,2]))
      }

      # Draw from background copula fell in upper rectangle
      else if (c0[i,1] < up_rect[2] &&
               c0[i,1] > up_rect[1]) {
        if (copula@flip_up)
          transformed <- c(psi1_rect_up_test(c2[i,1]), psi2_rect_up_test(c2[i,2]))
        else
          transformed <- c(psi1_rect_up_test(c1[i,1]), psi2_rect_up_test(c1[i,2]))
      }

      # Draw from background copula fell in neither rectangle
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
  v <- u[, 2,  drop = F]
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

    u_small_ind <- which(u<0.5)
    u_large_ind <- which(u>0.5)
    u_1 <- which(u==0.5)

    pdf <- u

    if (copula@flip_up){
      pdf[u_small_ind] <- dcylcop(cbind(2 * u[u_small_ind], v[u_small_ind]), c1)
      pdf[u_large_ind] <- dcylcop(cbind(1-(2 * u[u_large_ind]-1), v[u_large_ind]), c1)
      pdf[u_1] <- 1
    }else{
      pdf[u_small_ind] <- dcylcop(cbind(1-(2 * u[u_small_ind]), v[u_small_ind]), c1)
      pdf[u_large_ind] <- dcylcop(cbind(2 * u[u_large_ind]-1, v[u_large_ind]), c1)
      pdf[u_1] <- 1
    }
    pdf <- c(pdf)

    # pdf <- map2_dbl(u, v, function(u, v) {
    #   #(u,v) is in lower rectangle
    #   if (u < 0.5) {
    #     if (copula@flip_up)
    #       dcylcop(c(2 * u, v), c1)
    #     else
    #       dcylcop(c(2 * u, v), c2)
    #   }
    #
    #   #(u,v) is in upper rectangle
    #   else if (u > 0.5){
    #     if (copula@flip_up)
    #       dcylcop(c(2 * u - 1, v), c2)
    #     else
    #       dcylcop(c(2 * u - 1, v), c1)
    #   }
    #   #(u,v) is on the border, i.e. density of indepCopula
    #   else{
    #     1
    #   }
    # })
  }else{

    lo1 <- c(low_rect[1], 0)
    lo2 <- c(low_rect[2], 1)
    up1 <- c(up_rect[1], 0)
    up2 <- c(up_rect[2], 1)

    #background copula
    c0 <- copula@background.cop

    # non-rotated copula in the upper (flip_up==F) or lower (flip_up==T) rectangle
    c1 <- copula@sym.cop

    #C0-volume of the lower rectangle
    Vlo <- prob(copula@background.cop, lo1, lo2)
    #C0-volume of the upper rectangle
    Vup <- prob(copula@background.cop, up1, up2)

    u_low_ind <- intersect(which(u <= lo2[1]),which(u >= lo1[1]))
    u_up_ind <- intersect(which(u <= up2[1]),which(u >= up1[1]))
    u_backgr <- which(u%in%u[-c(u_low_ind,u_up_ind)])

    Vlou <- rep(NA,length(u))
    Vlov <- rep(NA,length(u))
    Vupu <- rep(NA,length(u))
    Vupv <- rep(NA,length(u))
    integ <- rep(NA,length(u))
    pdf <- rep(NA,length(u))

    Vlou[u_low_ind] <- u[u_low_ind] - low_rect[1]

    Vlov[u_low_ind] <-  pcylcop(cbind(low_rect[2], v[u_low_ind]), copula@background.cop) -
      pcylcop(cbind(low_rect[1], v[u_low_ind]), copula@background.cop)

    Vupu[u_up_ind] <- u[u_up_ind] - up_rect[1]

    Vupv[u_up_ind] <- pcylcop(cbind(up_rect[2], v[u_up_ind]), copula@background.cop) -
      pcylcop(cbind(up_rect[1], v[u_up_ind]), copula@background.cop)

    integ[u_low_ind] <- cylcop::ccylcop(cbind(lo2[1],v[u_low_ind]), copula@background.cop, cond_on=2, inverse=F)-
      cylcop::ccylcop(cbind(lo1[1],v[u_low_ind]), copula@background.cop, cond_on=2, inverse=F)
    integ[u_up_ind] <- cylcop::ccylcop(cbind(up2[1],v[u_up_ind]), copula@background.cop, cond_on=2, inverse=F)-
      cylcop::ccylcop(cbind(up1[1],v[u_up_ind]), copula@background.cop, cond_on=2, inverse=F)

    if (copula@flip_up){
      pdf[u_low_ind] <-
        (integ[u_low_ind] / Vlo) * dcylcop(cbind(Vlou[u_low_ind] / Vlo, Vlov[u_low_ind] / Vlo), c1)
      pdf[u_up_ind] <-
        (integ[u_up_ind] / Vup) * dcylcop(cbind(1-(Vupu[u_up_ind] / Vup), Vupv[u_up_ind] / Vup), c1)
      pdf[u_backgr] <- dcylcop(cbind(u[u_backgr], v[u_backgr]), c0)
    }else{
      pdf[u_low_ind] <-
        (integ[u_low_ind] / Vlo) * dcylcop(cbind(1-(Vlou[u_low_ind] / Vlo), Vlov[u_low_ind] / Vlo), c1)
      pdf[u_up_ind] <-
        (integ[u_up_ind] / Vup) * dcylcop(cbind(Vupu[u_up_ind] / Vup, Vupv[u_up_ind] / Vup), c1)
      pdf[u_backgr] <- dcylcop(cbind(u[u_backgr], v[u_backgr]), c0)
    }
  }
  return(pdf)
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


    u_small_ind <- which(u<0.5)
    u_large_ind <- which(u>0.5)
    u_1 <- which(u==0.5)

    cdf <- u

    if (copula@flip_up){
      cdf[u_small_ind] <- 0.5*pcylcop(cbind(2 * u[u_small_ind], v[u_small_ind]), c1)
      cdf[u_large_ind] <- 0.5*pcylcop(cbind(2 * u[u_large_ind]-1, v[u_large_ind]), c2)+ 0.5 * v[u_large_ind]
      cdf[u_1] <- u[u_1]*v[u_1]
    }else{
      cdf[u_small_ind] <- 0.5*pcylcop(cbind(2 * u[u_small_ind], v[u_small_ind]), c2)
      cdf[u_large_ind] <- 0.5*pcylcop(cbind(2 * u[u_large_ind]-1, v[u_large_ind]), c1) + 0.5 * v[u_large_ind]
      cdf[u_1] <- u[u_1]*v[u_1]
    }
    cdf <- c(cdf)
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

    u_low_ind <- intersect(which(u < low_rect[2]),which(u > low_rect[1]))
    u_up_ind <- intersect(which(u < up_rect[2]),which(u > up_rect[1]))
    u_backgr <- which(u%in%u[-c(u_low_ind,u_up_ind)])

    Vlou <- rep(NA,length(u))
    Vlov <- rep(NA,length(u))
    Vupu <- rep(NA,length(u))
    Vupv <- rep(NA,length(u))
    marg <- rep(NA,length(u))
    cdf <- rep(NA,length(u))



    Vlou[u_low_ind] <- u[u_low_ind] - low_rect[1]

    Vlov[u_low_ind] <-  pcylcop(cbind(low_rect[2], v[u_low_ind]), copula@background.cop) -
      pcylcop(cbind(low_rect[1], v[u_low_ind]), copula@background.cop)

    Vupu[u_up_ind] <- u[u_up_ind] - up_rect[1]

    Vupv[u_up_ind] <- pcylcop(cbind(up_rect[2], v[u_up_ind]), copula@background.cop) -
      pcylcop(cbind(up_rect[1], v[u_up_ind]), copula@background.cop)


    marg[u_low_ind] <-
      pcylcop(cbind(u[u_low_ind], 0), c0) + pcylcop(cbind(low_rect[1], v[u_low_ind]), c0) - pcylcop(c(low_rect[1], 0), c0)
    marg[u_up_ind] <-
      pcylcop(cbind(u[u_up_ind], 0), c0) + pcylcop(cbind(up_rect[1], v[u_up_ind]), c0) - pcylcop(c(up_rect[1], 0), c0)


    if (copula@flip_up){
      cdf[u_low_ind] <- Vlo * pcylcop(cbind(Vlou[u_low_ind] / Vlo, Vlov[u_low_ind] / Vlo), c1) + marg[u_low_ind]
      cdf[u_up_ind] <- Vup * pcylcop(cbind(Vupu[u_up_ind] / Vup, Vupv[u_up_ind] / Vup), c2) + marg[u_up_ind]
      cdf[u_backgr] <- pcylcop(cbind(u[u_backgr], v[u_backgr]), c0)
    }else{
      cdf[u_low_ind] <- Vlo * pcylcop(cbind(Vlou[u_low_ind] / Vlo, Vlov[u_low_ind] / Vlo), c2) + marg[u_low_ind]
      cdf[u_up_ind] <- Vup * pcylcop(cbind(Vupu[u_up_ind] / Vup, Vupv[u_up_ind] / Vup), c1) + marg[u_up_ind]
      cdf[u_backgr] <- pcylcop(cbind(u[u_backgr], v[u_backgr]), c0)
    }
  }
  return(cdf)
})



#' Condtional copula
#' @rdname ccylcop
# @describeIn cyl_rect_combine-class Calculate the conditional copula.
#' @export
setMethod("ccylcop", signature("cyl_rect_combine"), function(u,
                                                             copula,
                                                             cond_on = 2,
                                                             inverse = FALSE) {

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
#' @rdname set_cop_param
#@describeIn cyl_rect_combine-class Change attributes of existing object.
#' @export
setMethod("set_cop_param", "cyl_rect_combine", function(copula, param_val, param_name) {
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

  if (any(stringr::str_starts(param_name, "bg_|up_|low_", negate = TRUE))) {
    cop_param_val <-
      param_val[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = TRUE))]

    if (!any(is(copula@sym.cop) == "rotCopula")) {
      cop_param_num <-
        param_name[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = TRUE))] %>%
        match(copula@sym.cop@param.names)
      copula@sym.cop@parameters[cop_param_num] <- cop_param_val
    }
    else{
      cop_param_num <-
        param_name[which(stringr::str_starts(param_name, "bg_|up_|low_", negate = TRUE))] %>%
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
