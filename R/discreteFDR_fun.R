discrete.fdr.int <- function(
  pvec,
  pCDFlist,
  pCDFlist.indices,
  method = c("ADBH", "DBH", "DBR"),
  alpha = 0.05,
  method.parameter = NULL,
  crit.consts = FALSE,
  threshold = 1,
  data.name = NULL
) {
  # original number of hypotheses
  n <- length(pvec)
  
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input.data <- list()
  if(method != "DBR") {
    alg <- "Discrete Benjamini-Hochberg procedure"
    alg <- ifelse(method == "ADBH", paste("Adaptive", alg), alg)
    input.data$Method <- paste(alg, ifelse(method.parameter, "(step-up)", "(step-down)"))
  } else {
    input.data$Method <- paste0(
      "Discrete Blanchard-Roquain procedure (lambda = ", method.parameter, ")"
    )
  }
  input.data$raw.pvalues <- pvec
  if(length(pCDFlist) == n) {
    input.data$pCDFlist <- pCDFlist
  } else {
    idx <- unlist(pCDFlist.indices)
    pCDFlist.counts <- sapply(pCDFlist.indices, length)
    input.data$pCDFlist <- rep(pCDFlist, pCDFlist.counts)[order(idx)]
  }
  #input.data$pCDFlist    <- pCDFlist
  #if(!is.null(pCDFlist.indices)) {
  #  input.data$pCDFlist.indices <- pCDFlist.indices
  #}
  input.data$FDR.level   <- alpha
  input.data$Data.name   <- ifelse(
    !is.null(data.name),
    data.name,
    paste(deparse(substitute(pvec)), "and", deparse(substitute(pCDFlist)))
  )
  if(method == "DBR") input.data$DBR.Tuning <- method.parameter
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1) {
    # which p-values are not above threshold?
    select <- which(pvec <= threshold)
    # number of selected p-values
    m <- length(select)
    # filter pCDFs, indices and counts of selected p-values
    if(is.null(pCDFlist.indices)) {
      pCDFlist.indices <- as.list(select)
      pCDFlist.counts  <- rep(1, m)
      pCDFlist         <- pCDFlist[select]
    } else{
      pCDFlist.indices <- lapply(pCDFlist.indices, setdiff, seq_len(n)[-select])
      pCDFlist.counts  <- sapply(pCDFlist.indices, length)
      idx.nonempty     <- which(pCDFlist.counts > 0)
      pCDFlist         <- pCDFlist[idx.nonempty]
      pCDFlist.indices <- pCDFlist.indices[idx.nonempty]
      pCDFlist.counts  <- pCDFlist.counts[idx.nonempty]
    }
    pCDFlist.idx <- order(unlist(pCDFlist.indices))
    # rescale pCDFs
    F.thresh <- sapply(pCDFlist, function(X){t <- which(X <= threshold); if(length(t)) X[max(t)] else 0})
    pCDFlist <- sapply(seq_along(pCDFlist), function(k) pCDFlist[[k]] / F.thresh[k])
    # rescale selected p-values
    pvec <- pvec[select] / rep(F.thresh, pCDFlist.counts)[pCDFlist.idx]
  } else {
    # all p-values were selected
    select <- seq_len(n)
    m <- n
    # use original counts (or 1 for all, if all pCDFs are unique)
    pCDFlist.counts <- if(is.null(pCDFlist.indices)){
      rep(1, m)
    } else sapply(pCDFlist.indices, length)
    # F_i(1) = 1 for all i = 1, ..., n
    F.thresh <- rep(1.0, n)
  }
  
  #--------------------------------------------
  #       determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  support <- sort(unique(pmin(as.numeric(unlist(pCDFlist)), 1.0)))
  
  #--------------------------------------------
  #        compute significant p-values, their
  #        indices and the number of rejections
  #--------------------------------------------
  if(crit.consts) {
    switch(
      EXPR = method,
      ADBH = {
        res <- kernel_ADBH_crit(pCDFlist, support, sorted.pvals, method.parameter, alpha, pCDFlist.counts)
        crit.constants <- res$crit.consts
        idx.rej <- if(method.parameter) 
          which(sorted.pvals <= crit.constants) else 
            which(sorted.pvals > crit.constants)
      },
      DBH = {
        res <- kernel_DBH_crit(pCDFlist, support, sorted.pvals, method.parameter, alpha, pCDFlist.counts)
        crit.constants <- res$crit.consts
        idx.rej <- if(method.parameter) 
          which(sorted.pvals <= crit.constants) else 
            which(sorted.pvals > crit.constants)
      },
      DBR = {
        res <- kernel_DBR_crit(pCDFlist, support, sorted.pvals, method.parameter, alpha, pCDFlist.counts)
        crit.constants <- res$crit.consts
        idx.rej <- which(sorted.pvals <= crit.constants)
      }
    )
  } else {
    switch(
      EXPR = method,
      ADBH = {
        res <- kernel_ADBH_fast(pCDFlist, sorted.pvals, method.parameter, alpha, support, pCDFlist.counts)
        idx.rej <- if(method.parameter) 
          which(res <= seq_along(res) * alpha) else 
            which(res > seq_along(res) * alpha)
      },
      DBH = {
        res <- kernel_DBH_fast(pCDFlist, sorted.pvals, method.parameter, alpha, support, pCDFlist.counts)
        idx.rej <- if(method.parameter) 
          which(res <= seq_along(res) * alpha) else 
            which(res > seq_along(res) * alpha)
      },
      DBR = {
        res <- kernel_DBR_fast(pCDFlist, sorted.pvals, method.parameter, pCDFlist.counts)
        idx.rej <- which(res <= alpha)
      }
    )
  }
  
  k <- length(idx.rej)
  
  if(method == "DBR" || (method != "DBR" && method.parameter)) {
    if(k > 0) {
      m.rej <- if(method == "DBR") k else max(idx.rej)
      # determine significant (observed) p-values in sorted.pvals
      idx.rej <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- input.data$raw.pvalues[select][idx.rej]
    } else {
      m.rej <- 0
      idx.rej <- integer(0)
      pvec.rej <- numeric(0)
    }
  } else {
    if(k > 0) {
      m.rej <- min(idx.rej) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx.rej <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- input.data$raw.pvalues[select][idx.rej]
      } else {
        idx.rej <- numeric(0)
        pvec.rej <- numeric(0)
      }
    } else {
      m.rej <- m
      idx.rej <- 1:m
      pvec.rej <- input.data$raw.pvalues[select]
    }
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected = pvec.rej,
    Indices = select[idx.rej],
    Num.rejected = m.rej
  )
  
  # adjusted p-values
  if(method == "DBR" || (method != "DBR" && !method.parameter)){
    if(crit.consts){
      res <- res$pval.transf
    }
    # compute adjusted p-values
    pv.adj <- switch(
      EXPR = method,
      DBH  = cummax(pmin(res / 1:m, 1)),
      ADBH = cummax(pmin(res / 1:m, 1)),
      DBR  = rev(cummin(rev(pmin(res, 1))))
    ) 
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted          <- numeric(n)
    output$Adjusted[select]  <- pv.adj[ro]
    output$Adjusted[-select] <- NA
  }
  # add critical values to output list
  if(crit.consts) output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include selection data, if selection was applied
  if(threshold < 1){
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F.thresh
    output$Select$Pvalues <- input.data$raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec
    output$Select$Number <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input.data
  
  class(output) <- c("DiscreteFDR", class(output))
  return(output)
}
