#' @title
#' \[HSU\], \[HSD\], \[AHSU\] and \[AHSD\] procedures
#' 
#' @description
#' Apply the \[HSU\], \[HSD\], \[AHSU\] and \[AHSD\] procedures, with or without
#' computing the critical constants, to a set of p-values and their discrete
#' support.
#' 
#' @details
#' `DBH` and `ADBH` are wrapper functions for `discrete.BH`. `DBH` simply passes
#' all its parameters to `discrete.BH` with `adaptive = FALSE`. `ADBH` does the
#' same with `adaptive = TRUE`.
#' 
#' @seealso
#' [fast.Discrete], [DBR]
#' 
#' @templateVar pvalues FALSE
#' @templateVar stepUp FALSE
#' @templateVar alpha TRUE
#' @templateVar sorted_pv FALSE
#' @templateVar support FALSE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar direction TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar lambda FALSE
#' @templateVar adaptive TRUE
#' @templateVar threshold TRUE
#' @template param 
#' 
#' @references
#' Döhler, S., Durand, G., & Roquain, E. (2018). New FDR bounds for discrete
#'   and heterogeneous tests. *Electronic Journal of Statistics*, *12*(1),
#'   pp. 1867-1900. \doi{10.1214/18-EJS1441}
#'   
#' @template example
#' @examples 
#' DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#' summary(DBH.su.fast)
#' DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#' DBH.sd.fast$Adjusted
#' summary(DBH.sd.fast)
#' 
#' DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBH.su.crit)
#' DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd",
#'                    ret.crit.consts = TRUE)
#' DBH.sd.crit$Adjusted
#' summary(DBH.sd.crit)
#' 
#' ADBH.su.fast <- ADBH(raw.pvalues, pCDFlist)
#' summary(ADBH.su.fast)
#' ADBH.sd.fast <- ADBH(raw.pvalues, pCDFlist, direction = "sd")
#' ADBH.sd.fast$Adjusted
#' summary(ADBH.sd.fast)
#'
#' ADBH.su.crit <- ADBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(ADBH.su.crit)
#' ADBH.sd.crit <- ADBH(raw.pvalues, pCDFlist, direction = "sd",
#'                      ret.crit.consts = TRUE)
#' ADBH.sd.crit$Adjusted
#' summary(ADBH.sd.crit)
#' 
#' @templateVar DBR FALSE
#' @template return
#' 
#' @name discrete.BH
NULL

#' @rdname discrete.BH
#' @export
discrete.BH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", adaptive = FALSE, ret.crit.consts = FALSE, threshold = 1){
  # check arguments
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(threshold) || is.na(threshold) || !is.numeric(threshold) || threshold <= 0 || threshold > 1)
    stop("'threshold' must be a probability greater than 0 and limited to 1!")
  
  n <- length(raw.pvalues)
  if(n != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!")
  
  for(i in 1:n){
    if(!is.numeric(pCDFlist[[i]])){
      stop("All elements of 'pCDFlist' must be numeric vectors!")
    }#else pCDFlist[[i]] <- sort(unique(pmin(1, c(0, pCDFlist[[i]], 1))))
  }
  
  # check and match 'direction'
  direction <- match.arg(direction, c("su", "sd"))
  
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
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  direction <- match.arg(tolower(direction), c("su", "sd"))
  if(direction == "su"){
    # SU case
    if(ret.crit.consts){
      if(adaptive){
        # compute transformed support
        y <- kernel_ADBH_crit(pCDFlist, pv.list.all, sorted.pvals, TRUE, alpha)
      }
      else{
        # compute transformed support
        y <- kernel_DBH_crit(pCDFlist, pv.list.all, sorted.pvals, TRUE, alpha)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      if(adaptive){
        # compute transformed observed p-values
        y <- kernel_ADBH_fast(pCDFlist, sorted.pvals, TRUE, alpha, pv.list.all)
      }
      else{
        # compute transformed observed p-values
        y <- kernel_DBH_fast(pCDFlist, sorted.pvals, TRUE, alpha, pv.list.all)
      }
      # determine significant (transformed) p-values
      len.y <- length(y)
      if(len.y){
        idx <- which(y <= 1:len.y * alpha)
      }else{
        idx <- integer(0)
      }
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[select][idx]
    }
    else{
      m.rej <- 0
      idx <- integer(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    # SD case
    if(ret.crit.consts){
      if(adaptive){
        # compute transformed support
        y <- kernel_ADBH_crit(pCDFlist, pv.list.all, sorted.pvals, FALSE, alpha)
      }
      else{
        # compute transformed support
        y <- kernel_DBH_crit(pCDFlist, pv.list.all, sorted.pvals, FALSE, alpha)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      if(adaptive){
        # compute transformed sorted p-values
        y <- kernel_ADBH_fast(pCDFlist, sorted.pvals, FALSE)
      }
      else{
        # compute transformed sorted p-values
        y <- kernel_DBH_fast(pCDFlist, sorted.pvals, FALSE)
      }
      # determine significant (transformed) p-values
      idx <- which(y > 1:m * alpha) 
    }
    if(length(idx)){
      m.rej <- min(idx) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- raw.pvalues[select][idx]
      }
      else{
        idx <- numeric(0)
        pvec.rej <- numeric(0)
      }
    }
    else{
      m.rej <- m
      idx <- 1:m
      pvec.rej <- raw.pvalues[select]
    }
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = select[idx], Num.rejected = m.rej)
  if(direction == "sd"){
    if(ret.crit.consts){
      y <- y$pval.transf
    }
    # compute adjusted p-values
    pv.adj <- cummax(pmin(y / 1:m, 1))
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted <- numeric(n)
    output$Adjusted[select]  <- pv.adj[ro]
    output$Adjusted[-select] <- NA
  }
  # add critical values to output list
  if(ret.crit.consts) output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include details of the used algorithm as strings
  alg <- "Discrete Benjamini-Hochberg procedure"
  alg <- if(adaptive) paste("Adaptive", alg) else alg
  output$Method <- paste(alg, switch(direction, su = "(step-up)", sd = "(step-down)"))
  output$Signif.level <- alpha
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
}

#'@rdname discrete.BH
#'@export
DBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE, threshold = 1){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = FALSE, ret.crit.consts, threshold))
}

#'@rdname discrete.BH
#'@export
ADBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE, threshold = 1){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = TRUE, ret.crit.consts, threshold))
}
