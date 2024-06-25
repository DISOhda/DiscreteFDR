#' @title
#' Summarizing Discrete FDR Results
#' 
#' @description
#' `summary` method for class `DiscreteFDR`.
#'
#' @param object       an object of class `DiscreteFDR`.
#' @param x            an object of class `summary.DiscreteFDR`.
#' @param max          numeric or `NULL`, specifying the maximal number of
#'                     *rows* of the p-value table to be printed. By default,
#'                     when `NULL`, `getOption("max.print")` is used.
#' @param ...          further arguments passed to or from other methods.
#'
#' @details
#' `summary.DiscreteFDR` objects contain all data of an `DiscreteFDR` object,
#' but also include an additional table which includes the raw p-values,
#' their indices, the respective critical values (if present), the adjusted
#' p-values (if present) and a logical column to indicate rejection. The table
#' is sorted in ascending order by the raw p-values.
#' 
#' `print.summary.DiscreteFDR` simply prints the same output as
#' `print.DiscreteFDR`, but also prints the p-value table.
#' 
#' @return
#' `summary.DiscreteFDR` computes and returns a list that includes all the
#' data of an input `DiscreteFDR` object, plus
#' \item{Table}{`data.frame`, sorted by the raw p-values, that contains the
#'              indices, the raw p-values themselves, their respective critical
#'              values (if present), their adjusted p-values (if present) and a
#'              logical column to indicate rejection.}
#' 
#' @template exampleGPV
#' @examples 
#' DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#' summary(DBH.sd.crit)
#' 
#' @rdname summary.DiscreteFDR
#' @export
## S3 method for class 'DiscreteFDR'
summary.DiscreteFDR <- function(object, ...){
  if(!("DiscreteFDR" %in% class(object)))
    return(summary(object))
    
  # determine if selection as performed
  select <- exists('Select', object)
  if(select) m <- object$Select$Number
  
  # number of tests
  n <- length(object$Data$raw.pvalues)
  # determine order of raw p-values
  o <- order(object$Data$raw.pvalues)
  # ordered indices
  i <- 1:n
  # determine for each p-value if its corresponding null hypothesis is rejected
  r <- i %in% object$Indices #if(!select) o %in% object$Indices else o %in% object$Select.Indices[object$Indices]
  
  # create output object with summary table; include all data of object (DiscreteFDR)
  out <- c(object, list(Table = data.frame('Index' = i, 'P.value' = object$Data$raw.pvalues)))
  if(select){
    out$Table$Selected <- i %in% object$Select$Indices #rep(c(TRUE, FALSE), c(m, n - m))
    out$Table$Scaled <- NA
    out$Table$Scaled[out$Table$Selected] <- object$Select$Scaled
  }
  out$Table <- out$Table[o, ]
  if(exists('Critical.values', object)) {
    if(select) {
      out$Table$Critical.value <- NA
      out$Table$Critical.value[1:m] <- object$Critical.values[1:m][order(order(out$Table$Scaled[1:m]))]
    } else out$Table$Critical.value <- object$Critical.values
  }
  if(exists('Adjusted', object)){
    out$Table$Adjusted <- object$Adjusted[o]
  }
  out$Table <- data.frame(out$Table, 'Rejected' = r[o])
  rownames(out$Table) <- i
  
  # return output object
  class(out) <- "summary.DiscreteFDR" # basically a 'DiscreteFDR' object, but with a summary table (just like 'lm' and 'summary.lm' classes)
  return(out)
}

#'@rdname summary.DiscreteFDR
#'@export
## S3 method for class 'summary.DiscreteFDR'
print.summary.DiscreteFDR <- function(x, max = NULL, ...){
  if(!("summary.DiscreteFDR" %in% class(x)))
    return(print(x))
  
  # print 'DiscreteFDR' part of the object
  print.DiscreteFDR(x)
  
  # rows to print: number of rejections + 5 (if not requested otherwise)
  max <- if(!is.null(max)) ncol(x$Table) * max else getOption("max.print")
  
  # print additional summary table
  print(x$Table, max = max, ...)
  
  cat("\n")
  invisible(x)
}
