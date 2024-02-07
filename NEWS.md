# DiscreteFDR 1.3.7

* Introduction of `lifecycle` mechanisms.
* Marked `fast.Discrete()`, `fisher.pvalues.support()`, `match.pvals()`,
  `kernel_*()` and `amnesia` dataset as deprecated.
* Various documentation updates.
* Removal of links to `discreteMTP` packages, since it was removed from CRAN.


# DiscreteFDR 1.3.6

* Fixed a problem with `fisher.pvalues.support` that could cause p-values to 
be wrong or NA (Tanks to Iqraa Meah).
* Added GitHub.


# DiscreteFDR 1.3.5

* Fixed a problem with `fisher.pvalues.support` that could cause an infinite
loop when using `alternative = two.sided` (Thanks to Lukas Jansen).
* Changed version scheme from `x.y-z` to `x.y.z`


# DiscreteFDR 1.3-4

* Added a `NEWS.md` file to track changes to the package.
* Corrected a bug in `plot.DiscreteFDR` function that produced a false legend.
* Added argument plausibility checks to `discrete.BH` and `DBR` functions.