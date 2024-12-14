#' @name ADBH
#' 
#' @title
#' Wrapper Functions for the Adaptive Discrete Benjamini-Hochberg Procedure
#' 
#' @description
#' `ADBH()` is a wrapper function of [discrete.BH()] for computing \[AHSU\] and
#' \[AHSD\], which are more powerful than \[HSU\] and \[HSD\], respectively. It
#' simply passes its arguments to [discrete.BH()] with fixed `adaptive = TRUE`
#' and is computationally more demanding than [DBH()].
#' 
#' @template details_crit
#'  
#' @seealso
#' [`discrete.BH()`], [`DBH()`], [`DBR()`], [`DBY()`]
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
#' # ADBH (step-up) without critical values; using test results object
#' ADBH.su.fast <- ADBH(test.result)
#' summary(ADBH.su.fast)
#' 
#' # ADBH (step-down) without critical values; using extracted p-values 
#' # and supports
#' ADBH.sd.fast <- ADBH(raw.pvalues, pCDFlist, direction = "sd")
#' summary(ADBH.sd.fast)
#' 
#' # ADBH (step-up) with critical values; using extracted p-values and supports
#' ADBH.su.crit <- ADBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(ADBH.su.crit)
#' 
#' # ADBH (step-down) with critical values; using test results object
#' ADBH.sd.crit <- ADBH(test.result, direction = "sd", ret.crit.consts = TRUE)
#' summary(ADBH.sd.crit)
#' 
#' @export
ADBH <- function(test.results, ...) UseMethod("ADBH")

#' @rdname ADBH
#' @export
ADBH.default <- function(
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
    adaptive = TRUE,
    ret.crit.consts,
    select.threshold,
    pCDFlist.indices
  )
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(pCDFlist))
  )
  
  return(out)
}

#' @rdname ADBH
#' @export
ADBH.DiscreteTestResults <- function(
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
    adaptive = TRUE,
    ret.crit.consts,
    select.threshold
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}