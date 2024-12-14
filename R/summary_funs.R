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
#' DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd",
#'                    ret.crit.consts = TRUE)
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
  n <- length(object$Data$Raw.pvalues)
  # determine order of raw p-values
  o <- order(object$Data$Raw.pvalues)
  # ordered indices
  i <- seq_len(n)
  # determine for each p-value if its corresponding null hypothesis is rejected
  r <- i %in% object$Indices
  
  # create summary table
  tab <- data.frame('Index' = i, 'P.value' = object$Data$Raw.pvalues)
  if(select) {
    tab$Selected <- i %in% object$Select$Indices
    tab$Scaled <- NA
    tab$Scaled[tab$Selected] <- object$Select$Scaled
  }
  tab <- tab[o, ]
  if(exists('Critical.values', object)) {
    if(select) {
      ro <- order(order(tab$Scaled[seq_len(m)]))
      tab$Critical.value <- NA
      tab$Critical.value[seq_len(m)] <- object$Critical.values[seq_len(m)][ro]
    } else tab$Critical.value <- object$Critical.values
  }
  if(exists('Adjusted', object)){
    tab$Adjusted <- object$Adjusted[o]
  }
  tab <- data.frame(tab, 'Rejected' = r[o])
  # if row names are numbers, rearrange them to represent sort order
  if(all(rownames(tab) == tab$Index)) rownames(tab) <- i
  
  # return output object
  out <- c(object, list(Table = tab))
  class(out) <- "summary.DiscreteFDR"
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
