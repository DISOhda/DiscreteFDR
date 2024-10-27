#' @name DBR
#' 
#' @title
#' The Discrete Blanchard-Roquain Procedure
#' 
#' @description
#' Applies the \[HBR-\eqn{\lambda}\] procedure, with or without computing the
#' critical constants, to a set of p-values and their respective discrete
#' supports.
#' 
#' @details
#' \[DBR-\eqn{\lambda}\] is the discrete version of the 
#' \[Blanchard-Roquain-\eqn{\lambda}\] procedure (see References). The authors
#' of the latter suggest to take `lambda = alpha` (see their Proposition 17),
#' which explains the choice of the default value here.
#' @template details_crit 
#' 
#' @references:
#' G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control
#'   under independence and dependence. *Journal of Machine Learning Research*,
#'   *10*, pp. 2837-2871. \doi{10.48550/arXiv.0707.0536}
#'
#' @seealso
#' [`discrete.BH()`], [`DBH()`], [`ADBH()`], [`DBY()`]
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar test.results TRUE
#' @templateVar alpha TRUE
#' @templateVar lambda TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar select.threshold TRUE
#' @templateVar pCDFlist.indices TRUE
#' @templateVar triple.dots TRUE
#' @template param
#' 
#' @templateVar BR TRUE
#' @template return
#' 
#' @template exampleGPV
#' @examples
#' # DBR without critical values; using test results object
#' DBR.fast <- DBR(test.result)
#' summary(DBR.fast)
#' 
#' # DBR with critical values; using extracted p-values and supports
#' DBR.crit <- DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#' summary(DBR.crit)
#' 
#' @export
DBR <- function(test.results, ...) UseMethod("DBR")

#' @rdname DBR
#' @importFrom checkmate assert_character assert_integerish assert_list assert_numeric qassert
#' @export
DBR.default <- function(
  test.results,
  pCDFlist,
  alpha = 0.05,
  lambda = NULL,
  ret.crit.consts = FALSE,
  select.threshold = 1,
  pCDFlist.indices = NULL, 
  ...
) {
  # check arguments
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  # raw p-values
  qassert(x = test.results, rules = "N+[0, 1]")
  n <- length(test.results)
  
  # list structure of p-value distributions
  assert_list(
    x = pCDFlist,
    types = "numeric",
    any.missing = FALSE,
    min.len = 1,
    max.len = n
  )
  # individual p-value distributions
  for(i in seq_along(pCDFlist)) {
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
  
  # lambda
  if(is.null(lambda) || is.na(lambda)) {
    # if lambda is not provided, set lambda = alpha
    lambda <- alpha
  } else qassert(x = lambda, rules = "N1[0, 1]")
  
  # compute and return critical values?
  qassert(ret.crit.consts, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
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
  if(is.null(pCDFlist.indices)) {
    if(n != m) {
      stop(
        paste(
          "If no counts for the p-value CDFs are provided, the lengths of",
          "'test.results' and 'pCDFlist' must be equal!"
        )
      )
    }
    pCDFlist.indices <- as.list(seq_len(n))
    pCDFlist.counts <- rep(1, n)
  } else {
    set <- seq_len(n)
    for(i in seq_along(pCDFlist.indices)) {
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
  pvec <- match.pvals(test.results, pCDFlist, pCDFlist.indices)
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdr.int(
    pvec             = pvec,
    pCDFlist         = pCDFlist,
    pCDFlist.indices = pCDFlist.indices,
    method           = "DBR",
    alpha            = alpha,
    method.parameter = lambda,
    crit.consts      = ret.crit.consts,
    threshold        = select.threshold,
    data.name        = paste(
                         deparse(substitute(test.results)),
                         "and",
                         deparse(substitute(pCDFlist))
                       )
  )
  
  return(output)
}

#' @rdname DBR
#' @importFrom checkmate assert_character assert_r6 qassert
#' @export
DBR.DiscreteTestResults <- function(
  test.results,
  alpha = 0.05,
  lambda = NULL,
  ret.crit.consts = FALSE,
  select.threshold = 1, 
  ...
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
  
  # lambda
  if(is.null(lambda) || is.na(lambda)){
    # if lambda is not provided, set lambda = alpha
    lambda <- alpha
  }else qassert(x = lambda, rules = "N1[0, 1]")
  
  # compute and return critical values?
  qassert(ret.crit.consts, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdr.int(
    pvec             = test.results$get_pvalues(),
    pCDFlist         = test.results$get_pvalue_supports(unique = TRUE),
    pCDFlist.indices = test.results$get_support_indices(),
    method           = "DBR",
    alpha            = alpha,
    method.parameter = lambda,
    crit.consts      = ret.crit.consts,
    threshold        = select.threshold,
    data.name        = deparse(substitute(test.results))
  )
  
  return(output)
}