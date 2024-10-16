#' @name match.pvals
#' 
#' @keywords internal
#' 
#' @title
#' Matching Raw P-Values with Supports 
#' 
#' @description
#' Constructs the observed p-values from the raw observed p-values, by rounding
#' them to their nearest neighbor matching with the supports of their
#' respective CDFs (as in function `p.discrete.adjust()` of package
#' `discreteMTP`, which is no longer available on CRAN).
#' 
#' **Note**: This is an internal function and has to be called directly via
#' `:::`, i.e. `DiscreteFDR:::match.pvals()`.
#' 
#' @details
#' Well computed raw p-values should already belong to their respective CDF
#' support. So this function is called at the beginning of [`discrete.BH()`],
#' [`DBH()`], [`ADBH()`] and [`DBR()`], just in case raw p-values are biased.
#'
#' For each raw p-value that needs to be rounded, a warning is issued.
#'
#' @seealso
#' [`discrete.BH()`], [`DBR()`]
#'
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar pCDFlist.indices TRUE
#' @template param
#' 
#' @examples \dontrun{
#' toyList <- list(c(0.3,0.7,1),c(0.1,0.65,1))
#' toyRaw1 <- c(0.3,0.65)
#' match.pvals(toyRaw1, toyList)
#' toyRaw2 <- c(0.31,0.6)
#' match.pvals(toyRaw2, toyList)
#' }
#'
#' @return
#' A vector where each raw p-value has been replaced by its nearest neighbor, if
#' necessary.
#'
match.pvals <- function(test.results, pCDFlist, pCDFlist.indices = NULL) {
  m <- length(test.results)
  if(!is.null(pCDFlist.indices)) {
    idx <- unlist(pCDFlist.indices)
    counts <- sapply(pCDFlist.indices, length)
    pCDFlist <- rep(pCDFlist, counts)[order(idx)]
  }
  n <- length(pCDFlist)
  if(m > 0 && m == n) {
    pvec <- test.results
    in.CDF <- numeric(m)
    for(k in seq_len(m)) {
      in.CDF[k] <- match(pvec[k], pCDFlist[[k]])
      if(is.na(in.CDF[k])) {
        in.CDF[k] <- which.min(abs(pCDFlist[[k]] - pvec[k]))
        pvec[k] <- pCDFlist[[k]][in.CDF[k]]
        ordinal <- "th"
        if(k %% 10 == 1) ordinal <- "st"
        if(k %% 10 == 2) ordinal <- "nd"
        if(k %% 10 == 3) ordinal <- "rd"
        if(k %% 100 - k %% 10 == 10) ordinal <- "th"
        warning("Since ", test.results[k], 
                " is not a value of the CDF of the ", k, ordinal, " p-value,\n",
                "  the p-value is rounded to be ", pCDFlist[[k]][in.CDF[k]],
                call. = F)
      }
    }
    return(pvec)
  } else {
    stop("'pCDFlist' and 'test.results' do not match")
  }
}
