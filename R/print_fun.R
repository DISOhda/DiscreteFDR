#' @title
#' Printing DiscreteFDR results
#' 
#' @description
#' Prints the results of discrete FDR analysis, stored in a `DiscreteFDR` class
#' object.
#' 
#' @param x          an object of class "`DiscreteFDR`".
#' @param ...        further arguments to be passed to or from other methods.
#'                   They are ignored in this function.
#' 
#' @return
#' The input object `x` is invisibly returned via `invisible(x)`.
#' 
#' @template exampleGPV
#' @examples
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, direction = "su", ret.crit.consts = TRUE)
#' print(DBH.su.crit)
#' 
#' @importFrom stats p.adjust
#' @method print DiscreteFDR
#' @export
## S3 method for class 'DiscreteFDR'
print.DiscreteFDR <- function(x, ...){
  if(!any(c("DiscreteFDR", "summary.DiscreteFDR") %in% class(x)))
    return(print(x))
  
  n <- length(x$Data$raw.pvalues)
  k <- x$Num.rejected
  BH <- p.adjust(x$Data$raw.pvalues, "BH")
  
  # print title (i.e. algorithm)
  cat("\n")
  cat("\t", x$Data$Method, "\n")
  
  # print dataset name(s)
  cat("\n")
  cat("Data: ", x$Data$Data.name, "\n")
  
  # print short results overview
  if(!exists('Select', x))
    cat("Number of tests =", n, "\n") else{
      cat("Number of selected tests =", x$Select$Number, "out of", n, "\n")
      cat("Selection threshold =", x$Select$Threshold, "\n")
    }
  
  cat("Number of rejections =", k, "at global FDR level", x$Data$FDR.level, "\n")
  cat("(Original BH rejections = ", sum(BH <= x$Data$FDR.level), ")\n", sep = "")
  if(k && !exists('Select', x)) cat("Largest rejected p value: ", max(x$Rejected), "\n")
  
  cat("\n")
  invisible(x)
}
