#' @title
#' FDR-based Multiple Testing Procedures with Adaptation for Discrete Tests
#'
#' @description
#' This package implements the \[HSU\], \[HSD\],
#' \[AHSU\], \[AHSD\] and \[HBR-\eqn{\lambda}\] procedures for
#' discrete tests (see References). 
#'
#' @docType package
#' @import Rcpp
#' @useDynLib DiscreteFDRTest
#' @name DiscreteFDRTest
#'
#' @details
#' The functions are reorganized from the reference paper in the following way.
#' [`discrete.BH()`] (for Discrete Benjamini-Hochberg) implements
#' \[HSU\], \[HSD\], \[AHSU\] and \[AHSD\], while [`DBR()`] (for Discrete 
#' Blanchard-Roquain) implements \[HBR-\eqn{\lambda}\]. [`DBH()`] and [`ADBH()`]
#' are wrapper functions for [`discrete.BH()`] to access \[HSU\] and \[HSD\], as
#' well as \[AHSU\] and \[AHSD\] directly.
#' 
#' This package is part of a package family to which the
#' [`DiscreteDatasets`][DiscreteDatasets::DiscreteDatasets-package] and
#' [`DiscreteTests`][DiscreteTests::DiscreteTests-package] packages also
#' belong. The latter allows to compute p-values and their respective supports
#' for various tests. The objects that contain these results can be used
#' directly by the [`discrete.BH()`], [`DBH()`], [`ADBH()`] and [`DBR()`]
#' functions. Alternatively, these functions also accept a vector of raw
#' observed p-values and a list of the same length, whose elements are the
#' discrete supports of the CDFs of the p-values.
#' 
#' **Note**: The function [`fisher.pvalues.support()`], which allows to compute
#' such p-values and supports in the framework of a Fisher's exact test, is now
#' deprecated and should not be used anymore. It has been replaced by
#' [`generate.pvalues()`].
#' 
#' The same applies for the function [`fast.Discrete()`], which is a wrapper for
#' [`fisher.pvalues.support()`] and [`discrete.BH()`] and allows to apply
#' discrete procedures directly to a data set of contingency tables and do some
#' transformations before p-values are computed. It has been replaced by
#' [`direct.discrete.BH()`], but for more flexibility, users should employ
#' pipelines, e.g.\cr
#' `data |>`\cr
#' `  DiscreteDatasets::reconstruct_*(<args>) |>`\cr
#' `  DiscreteTests::*.test.pv(<args>) |>`\cr
#' `  discrete.BH(<args>)`.
#' 
#' @references
#' Döhler, S., Durand, G., & Roquain, E. (2018). New FDR bounds for discrete
#'   and heterogeneous tests. *Electronic Journal of Statistics*, *12*(1),
#'   pp. 1867-1900. \doi{10.1214/18-EJS1441}
#'   
#' G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control
#'   under independence and dependence. *Journal of Machine Learning Research*,
#'   *10*, pp. 2837-2871.
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
