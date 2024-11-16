discrete.fdr.int <- function(
  pvec,
  pCDFlist,
  pCDFlist.indices,
  method = c("ADBH", "DBH", "DBR", "DBY"),
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
  switch(
    EXPR = substr(method, nchar(method) - 2, nchar(method)),
    DBH = {
      alg <- "Discrete Benjamini-Hochberg procedure"
      alg <- ifelse(method == "ADBH", paste("Adaptive", alg), alg)
      input.data$Method <- paste(alg, ifelse(method.parameter, "(step-up)", "(step-down)"))
    },
    DBR = {
      input.data$Method <- paste0(
        "Discrete Blanchard-Roquain procedure (lambda = ", method.parameter, ")"
      )
    },
    DBY = {
      input.data$Method <- "Discrete Benjamini-Yekutieli procedure"
    }
  )
  input.data$Raw.pvalues <- pvec
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
      # keep CDFs of selected p-values
      pCDFlist         <- pCDFlist[select]
      # create indices according to selction
      pCDFlist_indices <- as.list(seq_len(m))
      # set counts (each CDF exists exactly once)
      pCDFlist_counts  <- rep(1, m)
    } else{
      # remove indices that were not selected
      pCDFlist_indices <- lapply(pCDFlist_indices, setdiff, seq_len(n)[-select])
      # determine counts
      pCDFlist_counts  <- sapply(pCDFlist_indices, length)
      # determine CDFs with non-zero counts
      idx_nonempty     <- which(pCDFlist_counts > 0)
      # keep CDFs with non-zero counts
      pCDFlist         <- pCDFlist[idx_nonempty]
      pCDFlist_indices <- pCDFlist_indices[idx_nonempty]
      pCDFlist_counts  <- pCDFlist_counts[idx_nonempty]
      # determine new indices (selection changes numbering!)
      new_idx          <- rep(NA, n)
      new_idx[select]  <- seq_len(m)
      # change indices according to selection
      pCDFlist_indices <- lapply(pCDFlist_indices, function(l) new_idx[l])
    }
    pCDFlist.idx <- order(unlist(pCDFlist.indices))
    # rescale pCDFs
    F.thresh <- sapply(pCDFlist, function(X) {t <- which(X <= threshold); if(length(t)) X[max(t)] else 0})
    pCDFlist <- sapply(seq_along(pCDFlist), function(k) pCDFlist[[k]] / F.thresh[k])
    # rescale selected p-values
    pvec <- pvec[select] / rep(F.thresh, pCDFlist.counts)[pCDFlist.idx]
  } else {
    # all p-values were selected
    select <- seq_len(n)
    m <- n
    # use original counts (or 1 for all, if all pCDFs are unique)
    pCDFlist.counts <- if(is.null(pCDFlist.indices)) {
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
      },
      DBY = {
        res <- kernel_DBY_crit(pCDFlist, support, sorted.pvals, alpha, pCDFlist.counts)
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
        res <- kernel_DBH_fast(pCDFlist, sorted.pvals, method.parameter, NULL, alpha, support, pCDFlist.counts)
        idx.rej <- if(method.parameter) 
          which(res <= seq_along(res) * alpha) else 
            which(res > seq_along(res) * alpha)
      },
      DBR = {
        res <- kernel_DBR_fast(pCDFlist, sorted.pvals, method.parameter, pCDFlist.counts)
        idx.rej <- which(res <= alpha)
      },
      DBY = {
        res <- kernel_DBY_fast(pCDFlist, sorted.pvals, pCDFlist.counts)
        idx.rej <- which(res <= seq_along(res) * alpha)
      }
    )
  }
  
  k <- length(idx.rej)
  
  if(method %in% c("DBR", "DBY") || (!(method %in% c("DBR", "DBY")) && method.parameter)) {
    if(k > 0) {
      m.rej <- if(method == "DBR") k else max(idx.rej)
      # determine significant (observed) p-values in sorted.pvals
      idx.rej <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- input.data$Raw.pvalues[select][idx.rej]
    } else {
      m.rej <- 0
      idx.rej <- integer(0)
      pvec.rej <- numeric(0)
    }
  } else {
    if(k > 0) {
      m.rej <- min(idx.rej) - 1
      if(m.rej) {
        # determine significant (observed) p-values in sorted.pvals
        idx.rej <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- input.data$Raw.pvalues[select][idx.rej]
      } else {
        idx.rej <- numeric(0)
        pvec.rej <- numeric(0)
      }
    } else {
      m.rej <- m
      idx.rej <- seq_len(m)
      pvec.rej <- input.data$Raw.pvalues[select]
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
  if(method %in% c("DBR", "DBY") || (!(method %in% c("DBR", "DBY")) && !method.parameter)) {
    if(crit.consts) {
      res <- res$pval.transf
    }
    # compute adjusted p-values
    pv.adj <- switch(
      EXPR = method,
      DBH  = cummax(pmin(res / seq_len(m), 1)),
      ADBH = cummax(pmin(res / seq_len(m), 1)),
      DBR  = rev(cummin(rev(pmin(res, 1)))),
      DBY  = rev(cummin(rev(pmin(res / seq_len(m), 1))))
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
  if(threshold < 1) {
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F.thresh
    output$Select$Pvalues <- input.data$Raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec
    output$Select$Number <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input.data
  
  class(output) <- c("DiscreteFDR", class(output))
  return(output)
}
