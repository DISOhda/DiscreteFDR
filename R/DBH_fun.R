#' @name DBH
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
#' [discrete.BH()], [ADBH()], [DBR()]
#' 
#' @templateVar x TRUE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar threshold TRUE
#' @templateVar pCDFlist.indices TRUE
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
#' @template example
#' @examples
#' 
#' DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#' summary(DBH.su.fast)
#' DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#' summary(DBH.sd.fast)
#' 
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBH.su.crit)
#' DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd",
#'                    ret.crit.consts = TRUE)
#' summary(DBH.sd.crit)
#' 
#' @export
DBH <- function(x, ...) UseMethod("DBH")

#' @rdname DBH
#' @export
DBH.default <- function(
  raw.pvalues,
  pCDFlist,
  alpha = 0.05,
  direction = "su",
  ret.crit.consts = FALSE,
  threshold = 1,
  pCDFlist.indices = NULL
) {
  out <- discrete.BH.default(
    raw.pvalues, 
    pCDFlist, 
    alpha, 
    direction, 
    adaptive = FALSE, 
    ret.crit.consts, 
    threshold, 
    pCDFlist.indices
  )
  
  out$Data$Data.name <- paste(
    deparse(substitute(raw.pvalues)),
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
  threshold = 1
) {
  out <- discrete.BH.DiscreteTestResults(test.results, alpha, direction, adaptive = FALSE, ret.crit.consts, threshold)
  out$Data$Data.name <- deparse(substitute(test.results))
  return(out)
}