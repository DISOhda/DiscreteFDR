#' @title
#' DiscreteFDRTest
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
#' [discrete.BH] (for Discrete Benjamini-Hochberg) implements
#' \[HSU\], \[HSD\], \[AHSU\] and \[AHSD\], while [DBR] (for Discrete 
#' Blanchard-Roquain) implements \[HBR-\eqn{\lambda}\]. [DBH] and [ADBH] are
#' wrapper functions for [discrete.BH] to access \[HSU\] and \[HSD\], as well as
#' \[AHSU\] and \[AHSD\] directly. Their main arguments are a vector of raw
#' observed p-values, and a list of the same length, whose elements are the
#' discrete supports of the CDFs of the p-values.
#' 
#' The function [fisher.pvalues.support] allows to compute such p-values and
#' support in the framework of a Fisher's exact test of association. It has been
#' inspired by a help page of the package `discreteMTP`, which is no longer
#' available on CRAN.
#' 
#' The function [fast.Discrete] is a wrapper for [fisher.pvalues.support] and
#' [discrete.BH] which allows to apply discrete procedures directly to a data
#' set of contingency tables.
#' 
#' We also provide the `amnesia` data set, used in our examples and in our
#' paper. It is basically the `amnesia` data set of package `discreteMTP` (no
#' longer on CRAN), but slightly reformatted.
#' 
#' No other function of the package should be used directly, as they are only
#' internal functions called by the main ones.
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
