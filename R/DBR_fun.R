#' @title
#' \[HBR-\eqn{\lambda}\] procedure
#' 
#' @description
#' Apply the \[HBR-\eqn{\lambda}\] procedure, with or without computing the
#' critical constants, to a set of p-values and their discrete support.
#' 
#' @details
#' \[DBR-\eqn{\lambda}\] is the discrete version of the 
#' \[Blanchard-Roquain-\eqn{\lambda}\] procedure (see References). The authors
#' of the latter suggest to take `lambda = alpha` (see their Proposition 17),
#' which explains the choice of the default value here. 
#' 
#' @section References:
#' G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control
#'   under independence and dependence. *Journal of Machine Learning Research*,
#'   *10*, pp. 2837-2871.
#'
#' @seealso
#' [discrete.BH]
#' 
#' @templateVar pvalues FALSE
#' @templateVar stepUp FALSE
#' @templateVar alpha TRUE
#' @templateVar sorted_pv FALSE
#' @templateVar support FALSE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar direction FALSE
#' @templateVar ret.crit.consts TRUE
#' @templateVar lambda TRUE
#' @templateVar adaptive FALSE
#' @templateVar threshold TRUE
#' @template param 
#' 
#' @template example
#' @examples
#' 
#' DBR.fast <- DBR(raw.pvalues, pCDFlist)
#' summary(DBR.fast)
#' DBR.crit <- DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBR.crit)
#' 
#' @templateVar DBR TRUE
#' @template return
#' 
#' @name DBR
#' @export
DBR <- function(raw.pvalues, pCDFlist, alpha = 0.05, lambda = NULL, ret.crit.consts = FALSE, threshold = 1){
  # check arguments
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  if (is.null(lambda)){
    # if lambda is not provided, set lambda = alpha
    lambda <- alpha
  }else{
    if(is.na(lambda) || !is.numeric(lambda) || lambda < 0 || lambda > 1)
      stop("'lambda' must be a probability between 0 and 1!")
  }
  if(is.null(threshold) || is.na(threshold) || !is.numeric(threshold) || threshold <= 0 || threshold > 1)
    stop("'threshold' must be a probability greater than 0 and limited to 1!")
  
  n <- length(raw.pvalues)
  if(n != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!")
  
  for(i in 1:n){
    if(!is.numeric(pCDFlist[[i]])){
      stop("All elements of 'pCDFlist' must be numeric vectors!")
    }#else pCDFlist[[i]] <- sort(unique(pmin(1, c(0, pCDFlist[[i]], 1))))
  }
  
  # prepare output data
  output.Data <- list()
  output.Data$raw.pvalues <- raw.pvalues
  output.Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output.Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(pCDFlist, raw.pvalues)
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1){
    select <- which(pvec <= threshold)
    m <- length(select)
    pCDFlist <- pCDFlist[select]
    F_thresh <- sapply(pCDFlist, function(X) X[max(which(X <= threshold))])
    pCDFlist <- sapply(1:m, function(k) pCDFlist[[k]] / F_thresh[k])
    pvec <- pvec[select] / F_thresh
  }else{
    select <- 1:n
    F_thresh <- rep(1, n)
    m <- n
  }
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)                                                                                       
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(pmin(as.numeric(unlist(pCDFlist)), threshold/max(F_thresh))))
  #--------------------------------------------
  #        Compute [DBR-lambda] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(ret.crit.consts){
    # compute transformed support
    y <- kernel_DBR_crit(pCDFlist, pv.list.all, sorted.pvals, lambda, alpha)
    # find critical constants
    crit.constants <- y$crit.consts
    y <- y$pval.transf
    idx <- which(sorted.pvals <= crit.constants)   
  }
  else{
    # compute transformed p-values
    y <- kernel_DBR_fast(pCDFlist, sorted.pvals, lambda)
    idx <- which(y <= alpha)
  }
  m.rej <- length(idx)
  if(m.rej){
    idx <- which(pvec <= sorted.pvals[m.rej]) 
    pvec.rej <- raw.pvalues[select][idx]
  }else{
    idx <- integer(0)
    pvec.rej <- numeric(0)
  }
  #--------------------------------------------
  #       Create output list
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = select[idx], Alpha = m.rej * alpha / m, Num.rejected = m.rej, Lambda = lambda)
  # compute adjusted p-values
  pv.adj <- rev(cummin(rev(pmin(y, 1))))
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted <- numeric(n)
  output$Adjusted[select]  <- pv.adj[ro]
  output$Adjusted[-select] <- NA
  # add critical values to output list
  if(ret.crit.consts) output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include details of the used algorithm as strings
  output$Method <- paste("Discrete Blanchard-Roquain procedure (lambda = ", lambda, ")", sep = "")
  output$Signif.level <- alpha
  output$Tuning <- lambda
  if(threshold < 1){
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F_thresh
    output$Select$Pvalues <- output.Data$raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec[order(output$Select$Pvalues)]
    output$Select$Number <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- output.Data
  
  class(output) <- "DiscreteFDR"
  return(output)
  
  return(output)
}