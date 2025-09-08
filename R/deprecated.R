#' Deprecated functions
#'
#' These functions are provided for compatibility with older version of
#' the cylcop package. They may eventually be completely
#' @rdname cylcop-deprecated
#' @name cylcop-deprecated
#' @param ... Parameters to be passed to the new versions of the functions
#' @docType package
#' @export  fit_steplength fit_angle scat_plot traj_plot circ_plot cop_scat_plot cop_plot make_traj qmixedvonmises mle.mixedvonmises
#' @aliases fit_steplength fit_angle scat_plot traj_plot circ_plot cop_scat_plot cop_plot make_traj qmixedvonmises mle.mixedvonmises
#' @section Details:
#' \tabular{rl}{
#'   \code{fit_angle()} \tab is replaced by \code{\link{fit_circ_param}()} and \code{\link{fit_circ_np}()}\cr
#'   \code{fit_steplength()} \tab is replaced by \code{\link{fit_lin_param}()} and \code{\link{fit_lin_np}()}\cr
#'   \code{scat_plot()} \tab is replaced by \code{\link{plot_joint_scat}()}\cr
#'   \code{traj_plot()} \tab is replaced by \code{\link{plot_track}()}\cr
#'   \code{circ_plot()} \tab is replaced by \code{\link{plot_joint_circ}()}\cr
#'   \code{cop_scat_plot()} \tab is replaced by \code{\link{plot_cop_scat}()}\cr
#'   \code{cop_plot()} \tab is replaced by \code{\link{plot_cop_surf}()}\cr
#'   \code{make_traj()} \tab is replaced by \code{\link{traj_sim}()}\cr
#'   \code{qmixedvonmises()} \tab is replaced by \code{\link{qvonmisesmix}()}\cr
#'   \code{mle.mixedvonmises()} \tab is replaced by \code{\link{mle.vonmisesmix}()}\cr
#' }


#helper function to get a list containing the default
#and non-default input parameters
get_inp_lst <- function(mc,defaults, non_default_names){

  # Define the ordered list of required and optional parameters
  param_order <- c(non_default_names, names(defaults))

  if(is.null(names(mc))){
    unnamed_args <- seq(1,length(mc))
  }else{
    unnamed_args <- which(names(mc) == "")
  }

  # Determine which parameters are available for unnamed arguments
  unnamed_avail <- param_order[which(!param_order %in%names(mc))]

  # Loop through unnamed arguments and assign them the appropriate names
  name_vec <- names(mc)
  j <- 1
  for (i in unnamed_args) {
    name_vec[i] <- unnamed_avail[j]
    j <- j+1
  }

  #Update the names of the captured arguments
  names(mc) <- name_vec

 #Ensure the required parameters are provided
  for(i in seq_along(non_default_names)){
    if (!non_default_names[i] %in% names(mc)) {
      stop(paste0(non_default_names[i], " is a required parameter."))
    }
  }

  # Combine the defaults with the provided arguments,
  # giving precedence to user-specified values
  mc <- modifyList(defaults, mc)

  return(mc)
}


fit_steplength <- function(...){
  # Capture the arguments passed via ...
  mc <- list(...)
  mc <- get_inp_lst(mc,
                    defaults = list(start=NULL, bandwidth=NULL, ncomp=2),
                    non_default_names=c("x", "parametric"))

  if(isFALSE(mc$parametric)){
    .Deprecated("fit_lin_np", package = "cylcop")

    fit_lin_np(x = mc$x, bandwidth = mc$bandwidth, limits = c(0,Inf))
  }else{
    .Deprecated("fit_lin_param", package = "cylcop")
    fit_lin_param(x = mc$x, densfun = mc$parametric, start = mc$start, ncomp = mc$ncomp)
  }
}

fit_angle <- function(...){
  mc <- list(...)
  mc <- get_inp_lst(mc,
                    defaults = list(bandwidth=NULL, mu=NULL, ncomp=2),
                    non_default_names=c("theta", "parametric"))
  if(isFALSE(mc$parametric)){
    .Deprecated("fit_circ_np", package = "cylcop")
    fit_circ_np(theta = mc$theta, bandwidth = mc$bandwidth)
  }else{
    .Deprecated("fit_circ_param", package = "cylcop")
    fit_circ_param(theta = mc$theta, densfun = mc$parametric, mu = mc$mu, ncomp = mc$ncomp)
  }
}

scat_plot <- function(...) {
  .Deprecated("plot_joint_scat", package = "cylcop")
  mc <- list(...)
  traj <- mc$traj
  if(is.null(traj)){
    traj <- mc[[1]]
    }
  plot_joint_scat(traj = traj, periodic = TRUE)
}


traj_plot <- function(...) {
  .Deprecated("plot_track", package = "cylcop")
  mc <- list(...)
  traj <- mc$traj
  if(is.null(traj)){
    traj <- mc[[1]]
  }
  plot_track(traj=traj)
}


circ_plot <- function(...) {
  .Deprecated("plot_joint_circ", package = "cylcop")
  mc <- list(...)
  traj <- mc$traj
  if(is.null(traj)){
    traj <- mc[[1]]
  }
  plot_joint_circ(traj=traj)
}

cop_scat_plot <- function(...) {

    .Deprecated("plot_cop_scat", package = "cylcop")
  mc <- list(...)
  input <- mc[[1]]

  if ((any(is(input) == "cyl_copula"))||(any(is(input) == "Copula"))){
    sample <- rcylcop(10000,input)
    traj <- data.frame(cop_u=sample[,1], cop_v=sample[,2])
  } else if (is.data.frame(input)){
    if(!all(c("cop_u","cop_v") %in% colnames(input))){
      stop(error_sound(), "trajectory must contain the columns 'cop_u' and 'cop_v'")
    }
    traj <- input
  }
  else{
    stop(error_sound(), "Provide either a (cyl_)copula object or a trajectory data.frame")
  }

  plot_theme <- list(
    geom_point(size=0.01,alpha=0.5),
    theme(legend.position = "none"),
    theme_bw(),
    xlab("v"),
    ylab("u"),
    theme(
      axis.title = element_text(size = 12, colour = "black"),
      axis.text = element_text(size = 10, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    )
  )

  p <-
    ggplot(traj, aes(x = .data$cop_v, y = .data$cop_u)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    plot_theme +
    coord_fixed()
  return(suppressWarnings(cowplot::ggdraw(p)))
}


  cop_plot <- function(...) {
    .Deprecated("plot_cop_surf", package = "cylcop")
    plot_cop_surf(...)
  }


  make_traj <-
    function(...) {

      .Deprecated("traj_sim", package = "cylcop")
    }


  qmixedvonmises <- function(...) {
    .Deprecated("qvonmisesmix", package = "cylcop")
  }


  mle.mixedvonmises <-  function(...)  {
    .Deprecated("mle.vonmisesmix", package = "cylcop")
  }



