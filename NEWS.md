
# cylcop 0.2.1

## Bug fixes

* Several plotting routines used the function `expression()` as a text label in
ggplot. This is not allowed and was historically possible in older versions of `ggplot2`
only by accident. It has been replaced by using a regular character vector and
`parse=T`.

* Missing package anchors were added to Rd \link{} targets in the documentations.



# cylcop 0.2.0

## New features

* New plotting functions: circular boxplots with `plot_joint_box()` and circular
histograms with `plot_circ_hist()`.

* New functions to assess the goodness of fit of copulas.
  + `cramer_vonmises()` calculates the Cramér-von-Mises criterion and the 
  associated p-value to compare a parametric copula to the empirical copula 
  of the data.
  + `wasserstein()` calculates the Wasserstein distance between two  copula PDFs
  or between a copula PDF and pseudo-observations.

* New functions for density, distribution, quantiles and random number generation of 
linear mixture distributions with an arbitrary number of components.
  + Mixed Weibull `dweibullmix()`, `pweibullmix()`, `qweibullmix()`, and `rweibullmix()`.
  + Mixed gamma `dgammamix()`, `pgammamix()`, `qgammamix()`, and `rgammamix()`.
  + Mixed normal `dnormmix()`, `pnormmix()`, `qnormmix()`, and `rnormmix()`.
  + Mixed log-normal `dlnormmix()`, `plnormmix()`, `qlnormmix()`, and `rlnormmix()`.

* New functions for density, distribution, quantiles and random number generation of 
von Mises mixture distributions with an arbitrary number of components:
`dvonmisesmix()`, `pvonmisesmix()`, `qvonmisesmix()`, and `rvonmisesmix()`.
`qvonmisesmix()` replaces the function `qmixedvonmises()`, which is deprecated.

* New functions to calculate density and distribution and to generate samples 
of bivariate joint distributions that are defined in terms of two marginal 
distributions and a copula: `djoint()`, `pjoint()`, and `rjoint()`.

* New function to obtain trajectories from x- and y-coordinates, `traj_get()`

## Enhancement

* New checking routines that check every argument of every function for validity.

* The plotting functions have new names and added and changed functionalities.
The old functions are deprecated.
  + `traj_plot()` is now `plot_track()`.
  It can now plot multiple trajectories and also separate x- and 
  y-coordinates.
  + `scat_plot()` is now `plot_joint_scat()`.
  It can take custom kernel density estimates for the margins as 
  input and the plots have an improved appearance. It can now take as input 
  trajectories as well as separate angles and step lengths
  + `cop_scat_plot()` is now `plot_cop_scat()`.
  It has new input arguments but with the same functionalities.
  + `circ_plot` is now `plot_joint_circ()`.
  It can now take as input 
  trajectories as well as separate angles and step lengths.
  + `cop_plot` is now `plot_cop_surf()`.
  
  * `make_traj()` (deprecated) is now called `traj_sim()`. 
  It can now take as argument
  directly the output of `fit_steplength() `or `fit_angle()`.

* The performance of many density, distribution, quantile and random number 
  generating functions is now improved. Especially for copulas of classes 
 `cyl_vonmises` and `cyl_rect_combine` and `cyl_cubsec`.

* The names of the arguments of `opt_lin_bw()` and `opt_circ_bw()`, have changed
 and are now more intuitive.

## Bug fixes

* `mi_cyl()` gives now exactly 1 for copulas with perfect correlation if
`normalize=TRUE` and `symmetrize=TRUE` are selected.

* `rdens()` (random number generation from kernel density estimates) 
now works correctly with circular density estimates.

* `fit_angle()` and `fit_steplength()` now work correctly with all possible 
parametric distributions and return the correct distribution names.

* When you plot a copula using `plot()` and provide additional arguments that go 
to the generic plotting function `base::plot()` there was an error that is now 
fixed.

# cylcop 0.1.0

First public release.
