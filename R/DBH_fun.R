#' @name DBH
#' 
#' @title
#' Wrapper functions for the Discrete Benjamini-Hochberg procedure
#' 
#' @description
#' 
#' `DBH()` is a wrapper function of [discrete.BH()] for computing \[HSU\] and
#' \[HSD\]. It simply passes its arguments to [discrete.BH()] with fixed
#' `adaptive = FALSE`.
#' 
#' @template details_crit
#' 
#' @seealso
#' [`discrete.BH()`], [`ADBH()`], [`DBR()`]
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar test.results TRUE
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar threshold TRUE
#' @templateVar pCDFlist.indices TRUE
#' @templateVar triple.dots TRUE
#' @template param
#' 
#' @templateVar DBR FALSE
#' @template return
#' 
#' @references
#' Döhler, S., Durand, G., & Roquain, E. (2018). New FDR bounds for discrete
#'   and heterogeneous tests. *Electronic Journal of Statistics*, *12*(1),
#'   pp. 1867-1900. \doi{10.1214/18-EJS1441}
#'   
#' @template exampleGPV
#' @examples
#' # DBH (SU) without critical values; using extracted p-values and supports
#' DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#' summary(DBH.su.fast)
#' 
#' # DBH (SD) without critical values; using extracted p-values and supports
#' DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#' summary(DBH.sd.fast)
#' 
#' # DBH (SU) with critical values; using test results
#' DBH.su.crit <- DBH(test.result, ret.crit.consts = TRUE)
#' summary(DBH.su.crit)
#' 
#' # DBH (SD) with critical values; using test results
#' DBH.sd.crit <- DBH(test.result, direction = "sd", ret.crit.consts = TRUE)
#' summary(DBH.sd.crit)
#' 
#' @export
DBH <- function(test.results, ...) UseMethod("DBH")

#' @rdname DBH
#' @export
DBH.default <- function(
  test.results,
  pCDFlist,
  alpha = 0.05,
  direction = "su",
  ret.crit.consts = FALSE,
  threshold = 1,
  pCDFlist.indices = NULL, 
  ...
) {
  out <- discrete.BH.default(
    test.results, 
    pCDFlist, 
    alpha, 
    direction, 
    adaptive = FALSE, 
    ret.crit.consts, 
    threshold, 
    pCDFlist.indices
  )
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(pCDFlist))
  )
  
  return(out)
}

#' @rdname DBH
#' @export
DBH.DiscreteTestResults <- function(
  test.results,
  alpha = 0.05,
  direction = "su",
  ret.crit.consts = FALSE,
  threshold = 1, 
  ...
) {
  out <- discrete.BH.DiscreteTestResults(test.results, alpha, direction, adaptive = FALSE, ret.crit.consts, threshold)
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}