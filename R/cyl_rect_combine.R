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
#' @slot cop.bg '\code{\linkS4class{cyl_vonmises}}' or
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
  slots = c("cop.lo",
            "cop.up",
            "cop.bg")
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
           copula2 = "rotated",
           background = indepCopula(),
           low_rect = c(0, 0.5),
           up_rect = "symmetric",
           flip = TRUE) {

    #validate input
    tryCatch({
      check_arg_all(list(check_argument_type(copula,
                                             type="Copula",
                                             dimension=2),
                         check_argument_type(copula,
                                             type="cyl_copula")
      )
      ,2)
      check_arg_all(list(check_argument_type(copula2,
                                             type="Copula",
                                             dimension=2),
                         check_argument_type(copula2,
                                             type="cyl_copula"),
                         check_argument_type(copula2,
                                             type="character",
                                             values = "rotated")
      )
      ,3)
      check_arg_all(list(check_argument_type(background,
                                             type="Copula",
                                             dimension=2),
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
      check_arg_all(list(check_argument_type(flip,
                                             type="logical"),
                         check_argument_type(flip,
                                             type="NULL")
      )
      ,2)
    },
    error = function(e) {
      error_sound()
      rlang::abort(conditionMessage(e))
    }
    )

    sym_rect <-  FALSE
    sym_bg <- TRUE
    sym_cop <- FALSE
    period_bg <- TRUE
    period_cop <- FALSE


    # Check rectangles

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

      if(low_rect[2]>up_rect[1]){
        stop(
          error_sound(),
          "Rectangles must not be overlapping!"
        )
      }

    if(low_rect[2]<low_rect[1]||up_rect[2]<up_rect[1]){
      stop(
        error_sound(),
        "Upper edge of a rectangle is below its lower edge."
      )
    }

    if (isTRUE(all.equal(low_rect[1], 1 - up_rect[2])) &&
        isTRUE(all.equal(low_rect[2], 1 - up_rect[1]))) {
      sym_rect <-  TRUE
    }


    # Get description of Copula objects in the rectangles

    #helper function to flip copulas
    flip_cop <- function(copula){

      #rotate copula in the other rectangle
      if (any(is(copula) == "cyl_vonmises")) {
        copula@flip <- !copula@flip
      }else if(any(is(copula) == "Copula")){
        if (any(is(copula) == "upfhCopula")) {
          copula <- fhCopula("lower")
        }else if (any(is(copula) == "lowfhCopula")) {
          copula <- fhCopula("upper")
        }else if(any(is(copula) == "rotCopula")){
        copula@flip[1] <- !copula@flip[1]
        if(all(copula@flip==c(F,F))){
          copula <- copula@copula
        }
          }else{
          copula <- rotCopula(copula, flip = c(TRUE, FALSE))
        }
      }else{
        stop()
      }
      return(copula)
    }


    name <- get_cop_name(copula)

    if (identical(copula2, "rotated")) {
      sym_cop <- TRUE
      name2 <- NULL
      copula2 <-     tryCatch(flip_cop(copula),
      error = function(e) {
        stop(
          error_sound(),
          paste("Rotating a", name, "object does not make sense. Please enter copula2 manually")
        )
      }
      )

      if(flip){
        temp <- copula2
        copula2 <- copula
        copula <- temp
      }

    }else{
      name2 <- get_cop_name(copula2)

      flipped_copula <- tryCatch(flip_cop(copula),
                                 error = function(e) {
                                   NULL
                                 }
      )
      if(identical(copula2,flipped_copula)){
        sym_cop <- TRUE
      }
      if(any(is(copula) == "cyl_quadsec")||
         any(is(copula) == "cyl_cubsec")||
         any(is(copula) == "cyl_rot_combine")){
        if(identical(copula2,copula)){
          sym_cop <- TRUE
        }
      }

    if(flip){
      if(!is.null(flipped_copula)){
        copula <- flipped_copula
      }else{
        stop(
          error_sound(),
          paste("Flipping a", name, "object does not make sense. Please flip the copulas manually")
        )
      }
      copula2 <-tryCatch(flip_cop(copula2),
                                      error = function(e) {
                                        stop(
                                          error_sound(),
                                          paste("Flipping a", name2, "object (copula2) does not make sense. Please flip the copulas manually")
                                        )
                                      }
      )
    }
    }

    # Get description of Copula objects not in the rectangles (i.e. in the "background")

    name_bg <- get_cop_name(background)



    if(any(is(background) == "Copula")){
      if (!any(is(background) == "indepCopula")) {
        period_bg <- FALSE
      }
    }else if(any(is(background) == "cyl_vonmises")){
      sym_bg <- FALSE
    }

    period_overall <- FALSE
    if((low_rect[1] >0) &&
       (up_rect[2] <1)){
      if(period_bg){
        period_overall <- T
      }
    }

    if(isTRUE(all.equal(low_rect[1], 0)) &&
       isTRUE(all.equal(up_rect[2], 1))){
      if(sym_cop){
        period_overall <- T
      }
    }

    if(!period_overall){
      if(period_bg){
      warning(warning_sound(),
              "Copulas in the rectangles are not symmetric to each other. Therefore the entire copula is not periodic.")
      }else{
        warning(warning_sound(),
                "Background is not periodic in u-direction. Therefore the entire copula is not periodic.")
      }
    }

    sym_overall <- T
    warning_text <- NULL
    if(!sym_bg){
      warning_text <- paste0(warning_text,"Background copula is not symmetric.\n" )
      sym_overall <- F
    }
    if(!sym_rect){
      warning_text <- paste0(warning_text,"Rectangles are not symmetric.\n" )
      sym_overall <- F
    }
    if(!sym_cop){
      warning_text <- paste0(warning_text,"Copulas are not mirror images of each other.\n" )
      sym_overall <- F
    }

    if(!sym_overall && period_overall){
      warning(warning_sound(),
              paste0(warning_text,"Therefore the entire copula is not symmetric"))
    }


    # Get parameters from copula objects

    #helper_function
    get_orig_cop <- function(copula){
      if(any(is(copula) == "rotCopula")){
        copula <- copula@copula
      }
      return(copula)
    }

    orig_copula <- get_orig_cop(copula)
    orig_copula2 <- get_orig_cop(copula2)
    orig_background <- get_orig_cop(background)
    parameters <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@parameters
      },
      if ("parameters" %in% methods::slotNames(orig_copula2)) {
        orig_copula2@parameters
      },
      if ("parameters" %in% methods::slotNames(orig_background)) {
        orig_background@parameters
      },
      low_rect, up_rect)
    param_names <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        paste0(orig_copula@param.names,"_lo")
      },
      if ("parameters" %in% methods::slotNames(orig_copula2)) {
        paste0(orig_copula2@param.names,"_up")
      },
      if ("parameters" %in% methods::slotNames(orig_background)) {
        paste0(orig_background@param.names,"_bg")
      },
      "low_rect1", "low_rect2", "up_rect1", "up_rect2")
    param_lowbnd <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.lowbnd
      },
      if ("parameters" %in% methods::slotNames(orig_copula2)) {
        orig_copula2@param.lowbnd
      },
      if ("parameters" %in% methods::slotNames(orig_background)) {
        orig_background@param.lowbnd
      },
      0, 0, 0, 0)
    param_upbnd <-
      c(if ("parameters" %in% methods::slotNames(orig_copula)) {
        orig_copula@param.upbnd
      },
      if ("parameters" %in% methods::slotNames(orig_copula2)) {
        orig_copula2@param.upbnd
      },
      if ("parameters" %in% methods::slotNames(orig_background)) {
        orig_background@param.upbnd
      },
      1, 1, 1, 1)


if(is.null(name2)){
  if(flip){
    full_name <- paste0("Rectangular patchwork of\nlower rectangle: upper copula rotated 90 degrees \n",
                        "upper rectangle: ",name, "\nbackground: ",name_bg)
  }else{
  full_name <- paste0("Rectangular patchwork of\nlower rectangle: ",
                      name,
                      "\nupper rectangle: lower copula rotated 90 degrees \nbackground: ",name_bg)
  }
}else{
  full_name <- paste0("Rectangular patchwork of\nlower rectangle: ", name, "\nupper rectangle: ",name2)
  if(flip){
    full_name <- paste0(full_name,"\nboth copulas are flipped")
  }
  full_name <- paste0(full_name,"\nbackground: ",name_bg)
}




    new(
      "cyl_rect_combine",
      name = full_name,
      parameters = parameters,     #contains all perameters of the copula in the rectangles
      #and of the copula outside the rectangles and the rectangles themselves
      param.names = param_names,   #dito
      param.lowbnd = param_lowbnd,
      param.upbnd = param_upbnd,
      cop.lo = copula,            # The copula in the lower rectangle
      cop.up = copula2,           # The copula in the upper rectangle
      cop.bg = background # The copula outside the rectangles
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
      any(is(copula@cop.bg) == "indepCopula")) {

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
    c0 <- rcylcop(n, copula@cop.bg)

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
    V_lo_rec <- low_rect[2]-low_rect[1]
    #C0-volume of the upper rectangle
    V_up_rec <- up_rect[2]-up_rect[1]

    #See text for expalantion of psi
    psi1_rect_lo_test <- function(x){
      x*V_lo_rec+low_rect[1]
    }
    psi1_rect_up_test <- function(x){
      x * V_up_rec + up_rect[1]
    }
    if(any(is(copula@cop.bg) == "indepCopula")){
      psi2_rect_lo_test <- function(x){
        x
      }
      psi2_rect_up_test <- psi2_rect_lo_test

    }else{
      psi2_rect_lo_test <-
        GoFKernel::inverse(
          f = function(y) {(pcylcop(c(low_rect[2], y), copula@cop.bg) - pcylcop(c(low_rect[1], y), copula@cop.bg)) / V_lo_rec},
          lower = 0,
          upper = 1
        )

      psi2_rect_up_test <-
        GoFKernel::inverse(
          f = function(y) {(pcylcop(c(up_rect[2], y), copula@cop.bg) - pcylcop(c(up_rect[1], y), copula@cop.bg)) / V_up_rec},
          lower = 0,
          upper = 1
        )
    }

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

precalc_cyl_rect <- function(u, v, copula){

  low_rect <-
    copula@parameters[match(c("low_rect1", "low_rect2"), copula@param.names)]
  up_rect <-
    copula@parameters[match(c("up_rect1", "up_rect2"), copula@param.names)]
  #background copula
  c0 <- copula@cop.bg

  # copula in lower rectangle
  c_lo <- copula@cop.lo

  #copula in upper rectangle
  c_up <- copula@cop.up

  if (isTRUE(all.equal(low_rect, c(0, 0.5))) &&
      isTRUE(all.equal(up_rect, c(0.5, 1))) &&
      (any(is(c0) == "indepCopula")||
       any(is(c0) == "cyl_quadsec")||
       any(is(c0) == "cyl_cubsec")||
       any(is(c0) == "cyl_rot_combine"))) {

    #Assign u on the border of the rectangle to rectangle 2.
    #doesn't mater, same value anyway
    u_low_ind <- which(u<0.5)
    u_up_ind <- which(u >= 0.5)
    u_backgr <- NULL

    lst <- list(
      simplified = T,
      u_low_ind = u_low_ind,
      u_up_ind = u_up_ind,
      u_backgr = u_backgr,
      low_rect = NULL,
      up_rect = NULL,
      c0 = c0,
      c_lo = c_lo,
      c_up = c_up,
      V_u = NULL,
      V_v = NULL,
      V_lo_rec = NULL,
      V_up_rec = NULL)

  }else{

    lo1 <- c(low_rect[1], 0)
    lo2 <- c(low_rect[2], 1)
    up1 <- c(up_rect[1], 0)
    up2 <- c(up_rect[2], 1)


    #C0-volume of the lower rectangle
    V_lo_rec <- low_rect[2]-low_rect[1]
    #C0-volume of the upper rectangle
    V_up_rec <- up_rect[2]-up_rect[1]

    u_low_ind <- intersect(which(u <= low_rect[2]),which(u >= low_rect[1]))
    u_up_ind <- intersect(which(u <= up_rect[2]),which(u >= up_rect[1]))
    non_backgr <- c(u_low_ind,u_up_ind)
    if(length(non_backgr)>0){
      u_backgr <- which(u%in%u[-non_backgr])
    }else{
      u_backgr <- seq(1,length(u))
    }

    V_u <- rep(NA,length(u))
    V_v <- rep(NA,length(u))
    marg <- rep(NA,length(u))


    if(length(u_low_ind)>0){
      V_u[u_low_ind] <- u[u_low_ind] - low_rect[1]
      V_v[u_low_ind] <-  pcylcop(cbind(low_rect[2], v[u_low_ind]), c0) -
        pcylcop(cbind(low_rect[1], v[u_low_ind]), c0)
    }

    if(length(u_up_ind)>0){
      V_u[u_up_ind] <- u[u_up_ind] - up_rect[1]
      V_v[u_up_ind] <- pcylcop(cbind(up_rect[2], v[u_up_ind]), c0) -
        pcylcop(cbind(up_rect[1], v[u_up_ind]), c0)
    }
    lst <- list(
      simplified = F,
      u_low_ind = u_low_ind,
      u_up_ind = u_up_ind,
      u_backgr = u_backgr,
      low_rect = low_rect,
      up_rect = up_rect,
      c0 = c0,
      c_lo = c_lo,
      c_up = c_up,
      V_u = V_u,
      V_v = V_v,
      V_lo_rec = V_lo_rec,
      V_up_rec = V_up_rec)
  }
  return(lst)
}


#' Calcualte density
#' @rdname Cylcop
# @describeIn cyl_rect_combine-class Calculate the density.
#' @export
setMethod("dcylcop", signature("matrix", "cyl_rect_combine"), function(u, copula) {
  v <- u[, 2,  drop = F]
  u <- u[, 1, drop = F]

  precalc <- precalc_cyl_rect(u ,v, copula)

  u_low_ind <- precalc$u_low_ind
  u_up_ind <- precalc$u_up_ind
  u_backgr <- precalc$u_backgr
  c_lo <- precalc$c_lo
  c_up <- precalc$c_up
  # The case when the rectangles together cover the entire unit square is treated separatley
  # because the equations are a lot simpler and faster to calculate (see text)

  if (precalc$simplified) {
    pdf <- u

    pdf[u_low_ind] <- dcylcop(cbind(2 * u[u_low_ind], v[u_low_ind]), c_lo)
    # pdf[u_up_ind] <- dcylcop(cbind(2-(2 * u[u_up_ind]), v[u_up_ind]), c_lo)
    pdf[u_up_ind] <- dcylcop(cbind(2 * u[u_up_ind]-1, v[u_up_ind]), c_up)
    pdf[u_backgr] <- 1

    pdf <- c(pdf)

  }else{
    c0 <- precalc$c0
    c_up <- precalc$c_up
    V_u <- precalc$V_u
    V_v <- precalc$V_v
    V_up_rec <- precalc$V_up_rec
    V_lo_rec <- precalc$V_lo_rec
    low_rect <- precalc$low_rect
    up_rect <- precalc$up_rect

    integ <- rep(NA,length(u))
    pdf <- rep(NA,length(u))

    if(length(u_low_ind)>0){
      integ <- cylcop::ccylcop(cbind(low_rect[2],v[u_low_ind]), c0, cond_on=2, inverse=F)-
        cylcop::ccylcop(cbind(low_rect[1],v[u_low_ind]), c0, cond_on=2, inverse=F)
      pdf[u_low_ind] <-
        (integ / V_lo_rec) * dcylcop(cbind(V_u[u_low_ind] / V_lo_rec, V_v[u_low_ind] / V_lo_rec), c_lo)
    }
    if(length(u_up_ind)>0){
      integ <- cylcop::ccylcop(cbind(up_rect[2],v[u_up_ind]), c0, cond_on=2, inverse=F)-
        cylcop::ccylcop(cbind(up_rect[1],v[u_up_ind]), c0, cond_on=2, inverse=F)
      pdf[u_up_ind] <-
        (integ / V_up_rec) * dcylcop(cbind((V_u[u_up_ind] / V_up_rec), V_v[u_up_ind] / V_up_rec), c_up)
    }
    if(length(u_backgr)>0){
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

  precalc <- precalc_cyl_rect(u ,v, copula)

  u_low_ind <- precalc$u_low_ind
  u_up_ind <- precalc$u_up_ind
  u_backgr <- precalc$u_backgr
  c_lo <- precalc$c_lo
  c_up <- precalc$c_up

  # The case when the rectangles together cover the entire unit square is treated separatley
  # because the equations are a lot simpler and faster to calculate (see text)

  if (precalc$simplified) {
    cdf <- u

    cdf[u_low_ind] <- 0.5*pcylcop(cbind(2 * u[u_low_ind], v[u_low_ind]), c_lo)
    cdf[u_up_ind] <- 0.5*pcylcop(cbind(2 * u[u_up_ind]-1, v[u_up_ind]), c_up)+0.5*v[u_up_ind]
    cdf[u_backgr] <- u[u_backgr]*v[u_backgr]

    cdf <- c(cdf)
  }


  # If background is not independent copula or rectangles are not from 0 to 0.5 and from 0.5 to 1
  # we need to do the full calculations

  else{
    c0 <- precalc$c0
    c_up <- precalc$c_up
    V_u <- precalc$V_u
    V_v <- precalc$V_v
    V_up_rec <- precalc$V_up_rec
    V_lo_rec <- precalc$V_lo_rec
    low_rect <- precalc$low_rect
    up_rect <- precalc$up_rect

    marg <- rep(NA,length(u))
    cdf <- rep(NA,length(u))

    if(length(u_low_ind)>0){
      marg <-
        pcylcop(cbind(low_rect[1], v[u_low_ind]), c0)
      cdf[u_low_ind] <- V_lo_rec * pcylcop(cbind(V_u[u_low_ind] / V_lo_rec, V_v[u_low_ind] / V_lo_rec), c_lo) + marg
    }
    if(length(u_up_ind)>0){
      marg <-
        pcylcop(cbind(up_rect[1], v[u_up_ind]), c0)
      cdf[u_up_ind] <- V_up_rec * pcylcop(cbind(V_u[u_up_ind] / V_up_rec, V_v[u_up_ind] / V_up_rec), c_up) + marg
    }
    if(length(u_backgr)>0){
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
  v <- u[, 2, drop = F]
  u <- u[, 1, drop = F]

  precalc <- precalc_cyl_rect(u ,v, copula)

  u_low_ind <- precalc$u_low_ind
  u_up_ind <- precalc$u_up_ind
  u_backgr <- precalc$u_backgr
  c_lo <- precalc$c_lo
  c0 <- precalc$c0

  # The case when the rectangles together cover the entire unit square is treated separatley
  # because the equations are a lot simpler and faster to calculate (see text)

  if (precalc$simplified) {
    cond <- u

    if(length(u_low_ind)>0){
      if(cond_on==1){
        cond[u_low_ind] <- ccylcop(cbind(2*u[u_low_ind],v[u_low_ind]), copula = c_lo, cond_on = 1,inverse = inverse)
      }else{
        cond[u_low_ind] <- 0.5*ccylcop(cbind(2*u[u_low_ind],v[u_low_ind]), copula = c_lo, cond_on = 2,inverse = inverse)
      }
    }
    if(length(u_up_ind)>0){
      if(cond_on==1){
        cond[u_up_ind] <- ccylcop(cbind(2-2*u[u_up_ind],v[u_up_ind]), copula = c_lo, cond_on = 1,inverse = inverse)
      }else{
        cond[u_up_ind] <- 1-0.5*ccylcop(cbind(2-2*u[u_up_ind],v[u_up_ind]), copula = c_lo, cond_on = 2,inverse = inverse)
      }
    }
    if(length(u_backgr)>0){
      cond[u_backgr] <- ccylcop(cbind(u[u_backgr],v[u_backgr]), copula = c0, cond_on = cond_on, inverse = inverse)
    }
    cond <- c(cond)
  }else{


    # If background is independent copula but rectangles are not from 0 to 0.5 and from 0.5 to 1
    c0 <- precalc$c0
    c_up <- precalc$c_up
    V_u <- precalc$V_u
    V_v <- precalc$V_v
    V_up_rec <- precalc$V_up_rec
    V_lo_rec <- precalc$V_lo_rec
    low_rect <- precalc$low_rect
    up_rect <- precalc$up_rect

    cond <- rep(NA,length(u))


    if(any(is(c0) %in% "indepCopula")){
      if(length(u_low_ind)>0){
        if(cond_on ==1){
          cond[u_low_ind] <- ccylcop(u=cbind((V_u[u_low_ind] / V_lo_rec), v[u_low_ind]), copula = c_lo,cond_on = 1, inverse = inverse)
        }else{
          cond[u_low_ind] <-  V_lo_rec *
            ccylcop(u=cbind((V_u[u_low_ind] / V_lo_rec), v[u_low_ind]), copula = c_lo, cond_on = 2, inverse = inverse)+
            low_rect[1]
        }
      }
      if(length(u_up_ind)>0){
        if(cond_on ==1){
          cond[u_up_ind] <- ccylcop(u=cbind(V_u[u_up_ind] / V_up_rec, v[u_up_ind]), copula=c_up, cond_on=1, inverse = inverse)
        }else{
          cond[u_up_ind] <-
            V_up_rec *
            ccylcop(u=cbind(V_u[u_up_ind] / V_up_rec, v[u_up_ind]), copula = c_up, cond_on = 2, inverse = inverse) +
            up_rect[1]
        }
      }
      if(length(u_backgr)>0){
        cond[u_backgr] <- cbind(u[u_backgr], v[u_backgr])[,3-cond_on,drop=T]
      }
    }else{
      numerical_flag <- F
      if(length(u_low_ind)>0){
        if(cond_on ==1){
          if(!inverse){
            cond[u_low_ind] <- ccylcop(u=cbind(V_u[u_low_ind] / V_lo_rec, V_v[u_low_ind] / V_lo_rec), copula = c_lo,cond_on = 1, inverse = F)
          }else{
            psi2_rect_lo_test <-
              GoFKernel::inverse(
                f = function(y) {(pcylcop(c(low_rect[2], y), c0) - pcylcop(c(low_rect[1], y), c0)) / V_lo_rec},
                lower = 0,
                upper = 1
              )
            for (i in u_low_ind) {
              cond[i] <-  psi2_rect_lo_test(ccylcop(u=cbind(V_u[i] / V_lo_rec, v[i]), copula = c_lo,cond_on = 1, inverse = T))
            }
          }
        }else{
          if(!inverse){
            cond[u_low_ind] <- (ccylcop(u=cbind(low_rect[2], v[u_low_ind]), copula = c0,cond_on = 2, inverse = inverse) -
                                  ccylcop(u=cbind(low_rect[1], v[u_low_ind]), copula = c0,cond_on = 2, inverse = inverse ))*
              ccylcop(u=cbind((V_u[u_low_ind] / V_lo_rec), (V_v[u_low_ind] / V_lo_rec)), copula = c_lo, cond_on = 2, inverse = inverse) +
              ccylcop(u=cbind(low_rect[1], v[u_low_ind]), copula = c0,cond_on = 2, inverse = inverse)
          }else{
            numerical_flag <- T
          }
        }
      }
      if(length(u_up_ind)>0){
        if(cond_on ==1){
          if(!inverse){
            cond[u_up_ind] <- ccylcop(u=cbind(V_u[u_up_ind] / V_up_rec, V_v[u_up_ind] / V_up_rec), copula=c_up, cond_on=1, inverse = F)
          }else{

            psi2_rect_lo_test <-
              GoFKernel::inverse(
                f = function(y) {(pcylcop(c(up_rect[2], y), c0) - pcylcop(c(up_rect[1], y), c0)) / V_up_rec},
                lower = 0,
                upper = 1
              )
            for (i in u_up_ind) {
              cond[i] <-  psi2_rect_lo_test(ccylcop(u=cbind(V_u[i] / V_up_rec, v[i]), copula = c_up,cond_on = 1, inverse = T))
            }

          }
        }else{
          if(!inverse){
            cond[u_up_ind] <- (ccylcop(u=cbind(up_rect[2], v[u_up_ind]), copula = c0,cond_on = 2, inverse = inverse ) -
                                 ccylcop(u=cbind(up_rect[1], v[u_up_ind]), copula = c0,cond_on = 2, inverse = inverse ))*
              ccylcop(u=cbind(V_u[u_up_ind] / V_up_rec, V_v[u_up_ind] / V_up_rec), copula = c_up,cond_on = 2, inverse = inverse) +
              ccylcop(u=cbind(up_rect[1], v[u_up_ind]), copula = c0,cond_on = 2, inverse = inverse)
          }else{
            numerical_flag <- T
          }
        }
      }
      if(length(u_backgr)>0){
        if(inverse && cond_on==2){
          numerical_flag <- T
        }else{
          cond[u_backgr] <- ccylcop(u=cbind(u[u_backgr], v[u_backgr]), copula=c0, cond_on=cond_on, inverse = inverse)
        }
      }
      if(numerical_flag){
        cond <- numerical_inv_conditional_cop(cbind(u,v),copula = copula, cond_on=2)
      }
    }
  }
  return(cond)

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
    if (!any(is(copula@cop.bg) == "rotCopula")) {
      bg_param_num <- param_name[which(stringr::str_starts(param_name, "bg_"))] %>%
        stringr::str_remove("bg_") %>%
        match(copula@cop.bg@param.names)
      copula@cop.bg@parameters[bg_param_num] <- bg_param_val
    }
    else{
      bg_param_num <- param_name[which(stringr::str_starts(param_name, "bg_"))] %>%
        stringr::str_remove("bg_") %>%
        match(copula@cop.bg@copula@param.names)
      copula@cop.bg@copula@parameters[bg_param_num] <-
        bg_param_val
    }
  }

  return(copula)
})
