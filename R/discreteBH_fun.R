#' @title
#' Discrete Benjamini-Hochberg procedure
#' 
#' @description
#' Apply the \[HSU\], \[HSD\], \[AHSU\] and \[AHSD\] procedures at a given FDR
#' level, with or without computing the critical constants, to a set of p-values
#' and their respective discrete supports.
#' 
#' @details
#' The adaptive variants \[AHSU\] and \[AHSD\], which are executed via
#' `adaptive = TRUE`, are often slightly more powerful than \[HSU\] and \[HSD\],
#' respectively. But they are also computationally more demanding.
#' @template details_crit
#' 
#' @seealso
#' [DBH()], [ADBH()], [DBR()]
#' 
#' @templateVar x TRUE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar threshold TRUE
#' @templateVar pCDFlist.indices TRUE
#' @template param
#' 
#' @templateVar DBR FALSE
#' @template return
#' 
#' @references
#' Döhler, S., Durand, G., & Roquain, E. (2018). New FDR bounds for discrete
#'   and heterogeneous tests. *Electronic Journal of Statistics*, *12*(1),
#'   pp. 1867-1900. \doi{10.1214/18-EJS1441}
#'   
#' @template example
#' @examples
#' 
#' DBH.su.fast <- discrete.BH(raw.pvalues, pCDFlist)
#' summary(DBH.su.fast)
#' DBH.sd.fast <- discrete.BH(raw.pvalues, pCDFlist, direction = "sd")
#' summary(DBH.sd.fast)
#' 
#' DBH.su.crit <- discrete.BH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBH.su.crit)
#' DBH.sd.crit <- discrete.BH(raw.pvalues, pCDFlist, direction = "sd",
#'                            ret.crit.consts = TRUE)
#' summary(DBH.sd.crit)
#' 
#' ADBH.su.fast <- discrete.BH(raw.pvalues, pCDFlist, adaptive = TRUE)
#' summary(ADBH.su.fast)
#' ADBH.sd.fast <- discrete.BH(raw.pvalues, pCDFlist, direction = "sd",
#'                             adaptive = TRUE)
#' summary(ADBH.sd.fast)
#'
#' ADBH.su.crit <- discrete.BH(raw.pvalues, pCDFlist, adaptive = TRUE, 
#'                             ret.crit.consts = TRUE)
#' summary(ADBH.su.crit)
#' ADBH.sd.crit <- discrete.BH(raw.pvalues, pCDFlist, direction = "sd", 
#'                             adaptive = TRUE, ret.crit.consts = TRUE)
#' summary(ADBH.sd.crit)
#' 
#' @name discrete.BH
#' @export
discrete.BH <- function(x, ...) UseMethod("discrete.BH")

#' @rdname discrete.BH
#' @importFrom checkmate assert_character assert_integerish assert_list assert_numeric qassert
#' @export
discrete.BH.default <- function(
  raw.pvalues,
  pCDFlist,
  alpha = 0.05,
  direction = "su",
  adaptive = FALSE,
  ret.crit.consts = FALSE,
  threshold = 1,
  pCDFlist.indices = NULL
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # raw p-values
  qassert(x = raw.pvalues, rules = "N+[0, 1]")
  n <- length(raw.pvalues)
  
  # list structure of p-value distributions
  assert_list(
    x = pCDFlist,
    types = "numeric",
    any.missing = FALSE,
    min.len = 1,
    max.len = n
  )
  # individual p-value distributions
  for(i in seq_along(pCDFlist)){
    assert_numeric(
      x = pCDFlist[[i]],
      lower = 0,
      upper = 1,
      any.missing = FALSE,
      min.len = 1,
      sorted = TRUE
    )
    if(max(pCDFlist[[i]]) != 1)
      stop("Last value of each vector in 'pCDFlist' must be 1!")
  }
  m <- length(pCDFlist)
  
  # significance level
  qassert(x = alpha, rules = "N1(0, 1]")
  
  # step-up/step-down direction
  assert_character(
    x = direction,
    n.chars = 2,
    len = 1,
    any.missing = FALSE
  )
  direction <- match.arg(tolower(direction), c("su", "sd"))
  
  # adaptiveness
  qassert(adaptive, "B1")
  
  # compute and return critical values?
  qassert(ret.crit.consts, "B1")
  
  # selection threshold
  qassert(x = threshold, rules = "N1(0, 1]")
  
  # list structure of indices
  assert_list(
    x = pCDFlist.indices,
    types = "numeric",
    any.missing = FALSE,
    len = m,
    unique = TRUE,
    null.ok = TRUE
  )
  # individual index vectors (if not NULL)
  if(is.null(pCDFlist.indices)){
    if(n != m){
      stop(
        paste(
          "If no counts for the p-value CDFs are provided, the lengths of",
          "'raw.pvalues' and 'pCDFlist' must be equal!"
        )
      )
    }
    pCDFlist.indices <- as.list(1:n)
    pCDFlist.counts <- rep(1, n)
  } else {
    set <- 1L:n
    for(i in seq_along(pCDFlist.indices)){
      pCDFlist.indices[[i]] <- assert_integerish(
        x = pCDFlist.indices[[i]],
        lower = 1,
        upper = n,
        any.missing = FALSE,
        min.len = 1,
        max.len = n,
        unique = TRUE,
        sorted = TRUE,
        coerce = TRUE
      )
      set <- setdiff(set, pCDFlist.indices[[i]])
    }
    if(length(set))
      stop("'pCDFlist.indices' must contain each p-value index exactly once!")
    pCDFlist.counts <- sapply(pCDFlist.indices, length)
  }
  
  #----------------------------------------------------
  #       check and prepare p-values for processing
  #----------------------------------------------------
  pvec <- match.pvals(pCDFlist, raw.pvalues, pCDFlist.indices)
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete_fdr_int(
    pvec = pvec,
    pCDFlist = pCDFlist,
    pCDFlist_indices = pCDFlist.indices,
    method = ifelse(adaptive, "ADBH", "DBH"),
    alpha = alpha,
    method_parameter = (direction == "su"),
    crit_consts = ret.crit.consts,
    threshold = threshold,
    data_name = paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  )
  
  return(output)
}

#' @rdname discrete.BH
#' @importFrom checkmate assert_character assert_r6 qassert
#' @export
discrete.BH.DiscreteTestResults <- function(
  test.results,
  alpha = 0.05,
  direction = "su",
  adaptive = FALSE,
  ret.crit.consts = FALSE,
  threshold = 1
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # discrete test results object
  assert_r6(
    x = test.results,
    classes = "DiscreteTestResults",
    public = c("get_pvalues", "get_pvalue_supports", "get_support_indices")
  )
  
  # significance level
  qassert(x = alpha, rules = "N1(0, 1]")
  
  # step-up/step-down direction
  assert_character(
    x = direction,
    n.chars = 2,
    len = 1,
    any.missing = FALSE
  )
  direction <- match.arg(tolower(direction), c("su", "sd"))
  
  # adaptiveness
  qassert(adaptive, "B1")
  
  # compute and return critical values?
  qassert(ret.crit.consts, "B1")
  
  # selection threshold
  qassert(x = threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete_fdr_int(
    pvec             = test.results$get_pvalues(),
    pCDFlist         = test.results$get_pvalue_supports(unique = TRUE),
    pCDFlist_indices = test.results$get_support_indices(),
    method           = ifelse(adaptive, "ADBH", "DBH"),
    alpha            = alpha,
    method_parameter = (direction == "su"),
    crit_consts      = ret.crit.consts,
    threshold        = threshold,
    data_name        = deparse(substitute(test.results))
  )
  
  return(output)
}