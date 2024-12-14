#' @name hist.DiscreteFDR
#' 
#' @title
#' Histogram of Raw P-Values
#' 
#' @description
#' Computes a histogram of the raw p-values of a `DiscreteFDR` object.
#' 
#' @param x          an object of class `DiscreteFDR`.
#' @param breaks     as in [`graphics::hist()`]; here, the Friedman-Diaconis
#'                   algorithm (`"FD"`) is used as default.
#' @param mode       single character string specifying for which $p$-values the
#'                   histogram is to be generated; must either be `"raw"` or
#'                   `"selected"`.
#' @param ...        further arguments to [`graphics::hist()`] or
#'                   [`graphics::plot.histogram()`], respectively.
#' 
#' @details
#' If `x` does not contain results of a selection approach, a warning is issued
#' and a histogram of the raw p-values is drawn. 
#' 
#' @return
#' An object of class `histogram`.
#'
#' @template exampleGPV
#' @examples
#' # DBH (SU)
#' DBH <- DBH(raw.pvalues, pCDFlist)
#' hist(DBH)
#'
#' @importFrom graphics hist
#' @export
hist.DiscreteFDR <- function(
    x,
    breaks = "FD",
    mode = c("raw", "selected"),
    ...
) {
  if(!("DiscreteFDR" %in% class(x)))
    stop("'x' must be an object of class DiscreteFDR")
  
  mode <- match.arg(tolower(mode), c("raw", "selected"))
  
  # determine if selection was performed
  select <- exists('Select', x)
  select.mode <- mode == "selected"
  # warn if histogram is for selected p-values that were never computed
  if(select.mode && !select)
    warning("No selected p-values present in object. Using raw p-values.")
  
  # get values of ...-arguments
  lst <- list(...)
  
  # p-value type
  pv.type <- ifelse(
    test = select.mode, 
    yes = ifelse(select, "Selected and Scaled", "Raw"),
    no = "Raw"
  )
  
  # labels
  if(!exists('main', where = lst))
    lst$main <- paste("Histogram of", pv.type, "P-Values")
  if(!exists('xlab', where = lst))
    lst$xlab <- paste(pv.type, "P-Values")
  
  # add 'breaks' and 'plot' to 'lst' for 'do.call'
  lst$breaks <- breaks
  lst$x <- if(mode == "selected" && select)
    x$Select$Scaled else
      x$Data$Raw.pvalues
  
  # call 'hist'
  fig <- do.call(hist, lst)
  
  # add input data variable name to result
  fig$xname <- deparse(substitute(x))
  
  # output
  if(ifelse(exists('plot', lst), lst$plot, TRUE))
    return(invisible(fig)) else
      return(fig)
}


#' @name plot.DiscreteFDR
#' @title Plot Method for `DiscreteFDR` objects
#'
#' @description
#' Plots raw p-values of a `DiscreteFDR` object and highlights rejected and
#' accepted p-values. If present, the critical values are plotted, too.
#'
#' @param x          object of class `DiscreteFDR`.
#' @param col        numeric or character vector of length 3 indicating the
#'                   colors of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present).
#'                   }
#' @param pch        numeric or character vector of length 3 indicating the
#'                   point characters of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present and `type.crit`
#'                           is a plot type like `'p'`, `'b'` etc.).
#'                   }
#' @param lwd        numeric vector of length 3 indicating the thickness of the
#'                   points and lines; defaults to current `par()$lwd` setting.
#' @param type.crit  1-character string giving the type of plot desired for the
#'                   critical values (e.g.: `'p'`, `'l'` etc; see [`plot()`]).
#' @param legend     if `NULL`, no legend is plotted; otherwise expecting a
#'                   character string like `"topleft"` etc. or a numeric vector
#'                   of two elements indicating (x, y) coordinates.
#' @param cex        numeric vector of length 3 indicating the size of point
#'                   characters or lines of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present).
#'                   }
#'                   defaults to current `par()$cex` setting.
#' @param ...        further arguments to [`plot.default()`].
#' 
#' @template exampleGPV
#' @examples 
#' DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' DBH.sd.fast <- DBH(test.result, direction = "sd")
#' DBH.sd.crit <- DBH(test.result, direction = "sd", ret.crit.consts = TRUE)
#' 
#' plot(DBH.sd.fast)
#' plot(DBH.sd.crit, xlim = c(1, 5), ylim = c(0, 0.4))
#' plot(DBH.su.fast, col = c(2, 4), pch = c(2, 3), lwd = c(2, 2), 
#'      legend = "topleft", xlim = c(1, 5), ylim = c(0, 0.4))
#' plot(DBH.su.crit, col = c(2, 4, 1), pch = c(1, 1, 4), lwd = c(1, 1, 2), 
#'      type.crit = 'o', legend = c(1, 0.4), lty = 1, xlim = c(1, 5), 
#'      ylim = c(0, 0.4))
#' 
#' @importFrom graphics plot lines points par
#' @importFrom checkmate assert check_character check_choice check_numeric
#' @importFrom checkmate assert_string
#' @export
plot.DiscreteFDR <- function(
  x,
  col = c(2, 4, 1),
  pch = c(20, 20, 17),
  lwd = rep(par()$lwd, 3),
  cex = rep(par()$cex, 3),
  type.crit = 'b',
  legend = NULL,
  ...
) {
  if(!("DiscreteFDR" %in% class(x)))
    stop("'x' must be an object of class DiscreteFDR")
  
  # make sure 'col' includes integers or color strings
  assert(
    check_character(col, min.len = 1, max.len = 3, null.ok = TRUE),
    check_numeric(col, min.len = 1, max.len = 3, null.ok = TRUE)
  )
  
  # make sure 'pch' includes integers or color strings
  assert(
    check_character(col, min.len = 1, max.len = 3, null.ok = TRUE),
    check_numeric(col, upper = 25, min.len = 1, max.len = 3, null.ok = TRUE)
  )
  
  # make sure 'type.crit' is a single character string
  assert_string(type.crit, null.ok = TRUE)
  
  # make sure 'legend' is a single string or a numerical vector of two values
  assert(
    check_choice(
      x = legend, 
      choices = c("bottomright", "bottom", "bottomleft", "left",
                  "topleft", "top", "topright", "right", "center"),
      null.ok = TRUE
    ),
    check_numeric(legend, len = 2, any.missing = FALSE)
  )
  
  # determine number of tests, selections and rejections
  n <- length(x$Data$Raw.pvalues)
  select <- exists('Select', x)
  m <- if(select) x$Select$Number else n
  
  # replace NAs in plot parameters with current par() settings
  col[is.na(col)] <- par()$col
  pch[is.na(pch)] <- par()$pch
  lwd[is.na(lwd)] <- par()$lwd
  cex[is.na(cex)] <- par()$cex
  
  # replicate shorter plot parameter vectors to avoid errors
  len <- 2L + as.integer(exists('Critical.values', where = x))
  col <- if(!is.null(col)) rep_len(col, len) else c( 2,  4,  1)[seq_len(len)]
  pch <- if(!is.null(pch)) rep_len(pch, len) else c(20, 20, 17)[seq_len(len)]
  lwd <- if(!is.null(lwd)) rep_len(lwd, len) else rep(par()$lwd, len)
  cex <- if(!is.null(cex)) rep_len(cex, len) else rep(par()$cex, len)
  
  # get values of ...-arguments
  lst <- list(...)
  
  # labels
  if(!exists('main', where = lst)) lst$main <- x$Data$Method
  if(!exists('xlab', where = lst)) {
    lst$xlab <- "Index"
    if(select) lst$xlab <- paste(lst$xlab, "(selected)")
  }
  if(!exists('ylab', where = lst)) {
    if(exists('Critical.values', where = x))
      lst$ylab <- "Sorted raw p-values / Critical values" else
        lst$ylab <- "Sorted raw p-values"
    if(select) lst$ylab <- paste(lst$ylab, "(selected, scaled)")
  }
  
  # start plotting with empty area
  lst$x <- if(select) x$Select$Scaled else x$Data$Raw.pvalues
  lst$col <- "NA"
  do.call(plot, lst)
  
  # plot accepted p-values
  if(select) {
    idx <- which(!(x$Select$Pvalues %in% x$Rejected))
    if(length(idx)) {
      y_acc <- sort(x$Select$Scaled[idx])
      x_acc <- (m - length(y_acc) + 1):m
    }
  } else {
    idx <- setdiff(1:n, x$Indices)
    if(length(idx)) {
      y_acc <- sort(x$Data$Raw.pvalues[idx])
      x_acc <- setdiff(seq_len(n), seq_len(x$Num.rejected))
    }
  }
  if(length(idx))
    points(x_acc, y_acc, col = col[2], pch = pch[2],
           lwd = lwd[2], cex = cex[2], ...)
  
  # plot rejected p-values
  if(x$Num.rejected) {
    y_rej <- if(select) sort(x$Select$Scaled[-idx]) else sort(x$Rejected)
    points(1:x$Num.rejected, y_rej, col = col[1], pch = pch[1],
           lwd = lwd[1], cex = cex[1], ...)
  }
  
  # plot critical values (if present and plotting is requested by the user)
  if(exists('Critical.values', where = x) && type.crit != 'n'){
    lines(x$Critical.values, col = col[3], lwd = lwd[3], pch = pch[3],
          type = type.crit, cex = cex[3], ...)
  }
  
  # plot legend
  if(!is.null(legend)) {
    idx <- 1:2
    lt <- rep(0, 3)
    if(exists('Critical.values', x) && type.crit != "n") {
      idx <- c(idx, 3)
      if(!(type.crit %in% c("p", "o", "b"))) pch[3] <- NA
      lt[3] <- if(type.crit %in% c('b', 'l', 'o')) {
        if(exists('lty', where = lst)) lst$lty else 1
      } else 0
    }
    len <- length(legend)
    if(len <= 2 & len >= 1) {
      if(len == 1){
        x <- legend
        y <- NULL
      } else {
        x <- legend[1]
        y <- legend[2]
      }
      legend(x, y, c("Rejected", "Accepted", "Critical values")[idx],
             col = c(col[idx], "darkgrey"), pch = pch[idx],
             lty = lt[idx], lwd = lwd[idx])
    } else 
      warning("Expecting character string or numeric vector of one or two ",
              "elements for creating a legend")
  }
}
