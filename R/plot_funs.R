#' @name hist.DiscreteFDR
#' 
#' @title
#' Histogram of Raw p-Values
#' 
#' @description
#' Computes a histogram of the raw p-values of a `DiscreteFDR` object.
#' 
#' @param x          an object of class "`DiscreteFDR`".
#' @param breaks     as in [hist]; here, the Friedman-Diaconis algorithm
#'                   (`"FD"`) is used as default.
#' @param plot       a boolean; if `TRUE` (the default), a histogram is plotted,
#'                   otherwise a list of breaks and counts is returned.
#' @param ...        further arguments to [hist] or [plot.histogram],
#'                   respectively.
#' 
#' @details
#' This method simply calls [hist] and passes the raw p-values of the object.
#' 
#' @return
#' An object of class `histogram`.
#'
#' @template example
#' @examples
#' 
#' DBH <- DBH(raw.pvalues, pCDFlist)
#' hist(DBH)
#'
#' @importFrom graphics hist
#' @export
hist.DiscreteFDR <- function(x, breaks = "FD", plot = TRUE, ...){
  # necessary to appease automated R CMD check on CRAN
  main <- xlab <- NULL
  lst <- list(...)
  
  # call 'hist' function with raw p.values (default breaks: "FD"); "..." passes
  # all additional 'hist' arguments to this call
  if(plot && !exists('main', where = lst) && !exists('xlab', where = lst))
    r <- hist(x$Data$raw.pvalues, breaks = breaks, main = "Histogram of raw p-values", xlab = "Raw p-values", plot = plot,  ...)
  else if(plot && !exists('main', where = lst))
    r <- hist(x$Data$raw.pvalues, breaks = breaks, main = "Histogram of raw p-values", plot = plot,  ...)
  else if(plot && !exists('xlab', where = lst))
    r <- hist(x$Data$raw.pvalues, breaks = breaks, xlab = "Raw p-values", plot = plot,  ...)
  else r <- hist(x$Data$raw.pvalues, breaks = breaks, plot = plot,  ...)
  
  r$xname <- deparse(substitute(x))
  
  if(plot) return(invisible(r)) else return(r)
}


#' @name plot.DiscreteFDR
#' @title Plot Method for `DiscreteFDR` objects
#'
#' @description
#' Plots raw p-values of a `DiscreteFDR` object and highlights rejected and
#' accepted p-values. If present, the critical values are plotted, too.
#'
#' @param x          an object of class "`DiscreteFDR`".
#' @param col        a numeric or character vector of length 3 indicating the
#'                   colors of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present).
#'                   }
#' @param pch        a numeric or character vector of length 3 indicating the
#'                   point characters of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present and `type.crit`
#'                           is a plot type like `'p'`, `'b'` etc.).
#'                   }
#' @param lwd        a numeric vector of length 3 indicating the thickness of
#'                   the points and lines.
#' @param type.crit  1-character string giving the type of plot desired for the
#'                   critical values (e.g.: `'p'`, `'l'` etc; see [plot]).
#' @param legend     if `NULL`, no legend is plotted; otherwise expecting a
#'                   character string like `"topleft"` etc. or a numeric vector
#'                   of two elements indicating (x, y) coordinates.
#' @param ...        further arguments to [plot.default].
#' 
#' @template example
#' @examples
#' 
#' DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#' DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#' 
#' plot(DBH.sd.fast)
#' plot(DBH.sd.crit, xlim = c(1, 5), ylim = c(0, 0.4))
#' plot(DBH.su.fast, col = c(2, 4), pch = c(2, 3), lwd = c(2, 2), 
#'      legend = "topleft", xlim = c(1, 5), ylim = c(0, 0.4))
#' plot(DBH.su.crit, col = c(2, 4, 1), pch = c(1, 1, 4), lwd = c(1, 1, 2), 
#'      type.crit = 'o', legend = c(1, 0.4), lty = 1, xlim = c(1, 5), 
#'      ylim = c(0, 0.4))
#' 
#' @importFrom graphics plot lines points title
#' @export
plot.DiscreteFDR <- function(x, col = c(2, 4, 1, 8), pch = rep(1, 4), lwd = rep(1, 4), type.crit = 'b', legend = NULL, ...){
  # determine number of tests, selections and rejections
  n <- length(x$Data$raw.pvalues)
  select <- exists('Select$Number', x)
  m <- if(select) x$Select$Number else n
  
  # replicate shorter plot parameter vectors to avoid errors
  col <- rep_len(col, 4)
  pch <- rep_len(pch, 4)
  lwd <- rep_len(lwd, 4)
  
  # get values of ...-arguments
  lst <- list(...)
  
  # start plotting with empty area
  if(exists('ylim', where = lst))
    plot(x$Data$raw.pvalues, col = NA, ...) else
      plot(x$Data$raw.pvalues, ylim = c(0, 1), col = NA, ...)
  if(exists('main', where = lst)) title(lst$main) else title(x$Method)
  if(exists('ylab', where = lst)) 
    title(ylab = lst$ylab) else title(ylab = "Sorted raw p-values")
  
  # plot critical values (if present and plotting is requested by the user)
  if(exists('Critical.values', where = x) && type.crit != 'n'){
    lines(x$Critical.values * if(select) x$Select$Effective.Thresholds else 1, 
          col = col[3], lwd = lwd[3], pch = pch[3], type = type.crit, ...)
  }
  
  idx <- match(x$Rejected, sort(x$Data$raw.pvalues))
  
  # plot accepted p-values
  points(setdiff(1:n, idx), sort(x$Data$raw.pvalues[setdiff(1:n, x$Indices)]), 
         col = col[2], pch = pch[2], lwd = lwd[2], ...)
  
  # plot rejected p-values
  if(x$Num.rejected)
    points(match(x$Rejected, sort(x$Data$raw.pvalues)), x$Rejected, 
           col = col[1], pch = pch[1], lwd = lwd[1], ...)
  
  # plot selection area (if threshold < 1)
  if(select){
    lines(rep(m + 0.05, 2), c(x$Select$Threshold, -10), col = "darkgrey", lwd = lwd[4], lty = 3)
    lines(c(-1000, m + 0.05), rep(x$Select$Threshold, 2), col = "darkgrey", lwd = lwd[4], lty = 3)
  }
  
  # plot legend
  if(!is.null(legend)){
    idx <- 1:2
    lt <- rep(0, 4)
    if(exists('Critical.values', x) && type.crit != "n"){
      idx <- c(idx, 3)
      if(!(type.crit %in% c("p", "o", "b"))) pch[3] <- NA
      lt[3] <- if(type.crit %in% c('b', 'l', 'o')){
        if(exists('lty', where = lst)) lst$lty else 1
      }else 0
    }
    if(select){
      idx <- c(idx, 4)
      pch[4] <- NA
      lt[4] <- 3
    }
    len <- length(legend)
    if(len <= 2 & len >= 1){
      if(len == 1){
        x <- legend
        y <- NULL
      }else{
        x <- legend[1]
        y <- legend[2]
      }
      legend(x, y, c("Rejected", "Accepted", "Critical values", "Selection area")[idx], col = c(col[idx], "darkgrey"), pch = pch[idx], lty = lt[idx], lwd = lwd[idx])
    }else warning("Expecting character string or numeric vector of one or two elements for creating a legend")
  }
}
