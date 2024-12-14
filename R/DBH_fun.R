#' @name DBH
#' 
#' @title
#' Wrapper Functions for the Discrete Benjamini-Hochberg Procedure
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
#' [`discrete.BH()`], [`ADBH()`], [`DBR()`], [`DBY()`]
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar test.results TRUE
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar select.threshold TRUE
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
#' # DBH (step-up) without critical values; using test results object
#' DBH.su.fast <- DBH(test.result)
#' summary(DBH.su.fast)
#' 
#' # DBH (step-down) without critical values; using extracted p-values 
#' # and supports
#' DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#' summary(DBH.sd.fast)
#' 
#' # DBH (step-up) with critical values; using extracted p-values and supports
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBH.su.crit)
#' 
#' # DBH (step-down) with critical values; using test results object
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
  select.threshold = 1,
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
    select.threshold, 
    pCDFlist.indices,
    ...
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
  select.threshold = 1, 
  ...
) {
  out <- discrete.BH.DiscreteTestResults(
    test.results,
    alpha,
    direction,
    adaptive = FALSE,
    ret.crit.consts,
    select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}