discrete_fdr_int <- function(
  pvec,
  pCDFlist,
  pCDFlist_indices,
  method = c("ADBH", "DBH", "DBR"),
  alpha = 0.05,
  method_parameter = NULL,
  crit_consts = FALSE,
  threshold = 1,
  data_name = NULL,
  output_R6 = FALSE
) {
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input_data             <- list()
  if(method != "DBR") {
    alg <- "Discrete Benjamini-Hochberg procedure"
    alg <- ifelse(method == "ADBH", paste("Adaptive", alg), alg)
    input_data$Method <- paste(alg, ifelse(method_parameter, "(step-up)", "(step-down)"))
  } else {
    input_data$Method <- paste0(
      "Discrete Blanchard-Roquain procedure (lambda = ", method_parameter, ")"
    )
  }
  input_data$raw.pvalues <- pvec
  input_data$pCDFlist    <- pCDFlist
  if(!is.null(pCDFlist_indices)) {
    input_data$pCDFlist.indices <- pCDFlist_indices
  }
  input_data$FDR.level   <- alpha
  input_data$Data.name   <- ifelse(
    !is.null(data_name),
    data_name,
    paste(deparse(substitute(pvec)), "and", deparse(substitute(pCDFlist)))
  )
  if(method == "DBR") input_data$DBR.Tuning <- method_parameter
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  # original number of hypotheses
  n <- length(pvec)
  if(threshold < 1){
    # which p-values are not above threshold?
    select <- which(pvec <= threshold)
    # number of selected p-values
    m <- length(select)
    # filter pCDFs, indices and counts of selected p-values
    if(is.null(pCDFlist_indices)){
      pCDFlist_indices <- as.list(select)
      pCDFlist_counts  <- rep(1, m)
      pCDFlist         <- pCDFlist[select]
    } else{
      pCDFlist_indices <- lapply(pCDFlist_indices, setdiff, seq_len(n)[-select])
      pCDFlist_counts  <- sapply(pCDFlist_indices, length)
      idx_nonempty     <- which(pCDFlist_counts > 0)
      pCDFlist         <- pCDFlist[idx_nonempty]
      pCDFlist_indices <- pCDFlist_indices[idx_nonempty]
      pCDFlist_counts  <- pCDFlist_counts[idx_nonempty]
    }
    pCDFlist_idx <- order(unlist(pCDFlist_indices))
    # rescale pCDFs
    F_thresh <- sapply(pCDFlist, function(X){t <- which(X <= threshold); if(length(t)) X[max(t)] else 0})
    pCDFlist <- sapply(seq_along(pCDFlist), function(k) pCDFlist[[k]] / F_thresh[k])
    # rescale selected p-values
    pvec <- pvec[select] / rep(F_thresh, pCDFlist_counts)[pCDFlist_idx]
  }else{
    # all p-values were selected
    select <- seq_len(n)
    m <- n
    # use original counts (or 1 for all, if all pCDFs are unique)
    pCDFlist_counts <- if(is.null(pCDFlist_indices)){
      rep(1, m)
    } else sapply(pCDFlist_indices, length)
    # F_i(1) = 1 for all i = 1, ..., n
    F_thresh <- rep(1.0, n)
  }
  
  #--------------------------------------------
  #       determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted_pvals <- pvec[o]
  
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  support <- sort(unique(pmin(as.numeric(unlist(pCDFlist)), 1.0)))
  
  #--------------------------------------------
  #        compute significant p-values, their
  #        indices and the number of rejections
  #--------------------------------------------
  if(crit_consts) {
    switch(
      EXPR = method,
      ADBH = {
        res <- kernel_ADBH_crit(pCDFlist, support, sorted_pvals, method_parameter, alpha, pCDFlist_counts)
        crit_constants <- res$crit.consts
        idx_rej <- if(method_parameter) 
          which(sorted_pvals <= crit_constants) else 
            which(sorted_pvals > crit_constants)
      },
      DBH = {
        res <- kernel_DBH_crit(pCDFlist, support, sorted_pvals, method_parameter, alpha, pCDFlist_counts)
        crit_constants <- res$crit.consts
        idx_rej <- if(method_parameter) 
          which(sorted_pvals <= crit_constants) else 
            which(sorted_pvals > crit_constants)
      },
      DBR = {
        res <- kernel_DBR_crit(pCDFlist, support, sorted_pvals, method_parameter, alpha, pCDFlist_counts)
        crit_constants <- res$crit.consts
        idx_rej <- which(sorted_pvals <= crit_constants)
      }
    )
  } else {
    switch(
      EXPR = method,
      ADBH = {
        res <- kernel_ADBH_fast(pCDFlist, sorted_pvals, method_parameter, alpha, support, pCDFlist_counts)
        idx_rej <- if(method_parameter) 
          which(res <= seq_along(res) * alpha) else 
            which(res > seq_along(res) * alpha)
      },
      DBH = {
        res <- kernel_DBH_fast(pCDFlist, sorted_pvals, method_parameter, alpha, support, pCDFlist_counts)
        idx_rej <- if(method_parameter) 
          which(res <= seq_along(res) / m * alpha) else 
            which(res > seq_along(res) / m * alpha)
      },
      DBR = {
        res <- kernel_DBR_fast(pCDFlist, sorted_pvals, method_parameter, pCDFlist_counts)
        idx_rej <- which(res <= alpha)
      }
    )
  }
  
  k <- length(idx_rej)
  
  if(method == "DBR" || (method != "DBR" && method_parameter)) {
    if(k > 0) {
      m_rej <- if(method == "DBR") k else max(idx_rej)
      # determine significant (observed) p-values in sorted_pvals
      idx_rej <- which(pvec <= sorted_pvals[m_rej]) 
      pvec_rej <- input_data$raw.pvalues[select][idx_rej]
    } else {
      m_rej <- 0
      idx_rej <- integer(0)
      pvec_rej <- numeric(0)
    }
  } else {
    if(k > 0) {
      m_rej <- min(idx_rej) - 1
      if(m_rej){
        # determine significant (observed) p-values in sorted_pvals
        idx_rej <- which(pvec <= sorted_pvals[m_rej])
        pvec_rej <- input_data$raw.pvalues[select][idx_rej]
      } else {
        idx_rej <- numeric(0)
        pvec_rej <- numeric(0)
      }
    } else {
      m_rej <- m
      idx_rej <- 1:m
      pvec_rej <- input_data$raw.pvalues[select]
    }
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected = pvec_rej,
    Indices = select[idx_rej],
    Num.rejected = m_rej
  )
  
  # adjusted p-values
  if(method == "DBR" || (method != "DBR" && !method_parameter)){
    if(crit_consts){
      res <- res$pval.transf
    }
    # compute adjusted p-values
    pv_adj <- switch(
      EXPR = method,
      DBH  = cummax(pmin(res * m / 1:m, 1)),
      ADBH = cummax(pmin(res / 1:m, 1)),
      DBR  = rev(cummin(rev(pmin(res, 1))))
    ) 
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted          <- numeric(n)
    output$Adjusted[select]  <- pv_adj[ro]
    output$Adjusted[-select] <- NA
  }
  # add critical values to output list
  if(crit_consts) output$Critical.values <- c(crit_constants, rep(NA, n - m))
  
  # include selection data, if selection was applied
  if(threshold < 1){
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F_thresh
    output$Select$Pvalues <- input_data$raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec
    output$Select$Number <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input_data
  
  class(output) <- "DiscreteFDR"
  return(output)
}
