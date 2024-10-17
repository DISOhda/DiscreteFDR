#' @title
#' Computing Discrete P-Values and Their Supports for Fisher's Exact Test
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Computes discrete raw p-values and their support for Fisher's exact test
#' applied to 2x2 contingency tables summarizing counts coming from two
#' categorical measurements.
#'
#' **Note**: This function is deprecated and will be removed in a future
#' version. Please use [`generate.pvalues()`] with
#' `test.fun = DiscreteTests::fisher.test.pv` and (optional) 
#' `preprocess.fun = DiscreteDatasets::reconstruct_two` or 
#' `preprocess.fun = DiscreteDatasets::reconstruct_four` instead. Alternatively,
#' use a pipeline like\cr
#' `data |>`\cr
#' `  DiscreteDatasets::reconstruct_*(<args>) |>`\cr
#' `  DiscreteTests::fisher.test.pv(<args>)`
#' 
#' @details
#' Assume that each contingency tables compares two variables and resumes the
#' counts of association or not with a condition. This can be resumed in the
#' following table:
#' \tabular{lccc}{
#' \tab Association \tab No association  \tab      Total      \cr
#'      Variable 1  \tab    \eqn{X_1}    \tab    \eqn{Y_1}    \tab \eqn{N_1} \cr
#'      Variable 2  \tab    \eqn{X_2}    \tab    \eqn{Y_2}    \tab \eqn{N_2} \cr
#'      Total       \tab \eqn{X_1 + X_2} \tab \eqn{Y_1 + Y_2} \tab \eqn{N_1 + N_2}
#' }
#' If `input="noassoc"`, `counts` has four columns which respectively contain,
#' \eqn{X_1}, \eqn{Y_1}, \eqn{X_2} and \eqn{Y_2}. If `input="marginal"`,
#' `counts` has four columns which respectively contain \eqn{X_1}, \eqn{N_1},
#' \eqn{X_2} and \eqn{N_2}.
#' 
#' If `input="HG2011"`, we are in the situation of the `amnesia` data set as
#' in Heller & Gur (2011, see References). Each contingency table is obtained
#' from one variable which is compared to all other variables of the study. That
#' is, counts for "second variable" are replaced by the sum of the counts of the
#' other variables:
#' \tabular{lccc}{
#' \tab Association            \tab No association            \tab Total                     \cr
#'      Variable \eqn{j}       \tab \eqn{X_j}                 \tab \eqn{Y_j}                 \tab \eqn{N_j} \cr
#'      Variables \eqn{\neq j} \tab \eqn{\sum_{i \neq j} X_i} \tab \eqn{\sum_{i \neq j} Y_i} \tab \eqn{\sum_{i \neq j} N_i} \cr
#'      Total                  \tab \eqn{\sum X_i}            \tab \eqn{\sum Y_i}            \tab \eqn{\sum N_i}
#' }
#' Hence `counts` needs to have only two columns which respectively contain \eqn{X_j} and \eqn{Y_j}.
#'
#' The code for the computation of the p-values of Fisher's exact test is
#' inspired by the example in the help page of `p.discrete.adjust` of package
#' `discreteMTP`, which is no longer available on CRAN.
#'
#' See the Wikipedia article about Fisher's exact test, paragraph Example, for
#' a good depiction of what the code does for each possible value of
#' `alternative`.
#'
#' @seealso
#' [`fisher.test()`]
#' 
#' @param counts        a data frame of two or four columns and any number of
#'                      lines; each line represents a 2x2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument, see
#'                      Details.
#' @param alternative   same argument as in [`stats::fisher.test()`]. The three
#'                      possible values are `"greater"` (default),
#'                      `"two.sided"` or `"less"` and you can specify
#'                      just the initial letter.
#' @param input         the format of the input data frame, see Details. The
#'                      three possible values are `"noassoc"` (default),
#'                      `"marginal"` or `"HG2011"` and you can specify
#'                      just the initial letter.
#' 
#' @template exampleFPV
#' 
#' @return
#' A list of two elements:
#' \item{raw}{raw discrete p-values.}
#' \item{support}{a list of the supports of the CDFs of the p-values.
#' Each support is represented by a vector in increasing order.}
#' 
#' @references
#' R. Heller and H. Gur (2011). False discovery rate controlling procedures for
#'   discrete tests. arXiv preprint.
#'   [arXiv:1112.4627v2](https://arxiv.org/abs/1112.4627v2).
#'
#' "Fisher's exact test", Wikipedia, The Free Encyclopedia, accessed 2018-03-20,
#' [link](https://en.wikipedia.org/w/index.php?title=Fisher's_exact_test&oldid=823327889).
#' 
#' @importFrom stats dhyper phyper pbinom
#' @importFrom lifecycle deprecate_soft
#' @importFrom checkmate assert_data_frame qassert
#' @importFrom DiscreteDatasets reconstruct_four reconstruct_two
#' @import DiscreteTests
#' @export
fisher.pvalues.support <- function(counts, alternative = "greater", input = "noassoc"){
  deprecate_soft("2.0.0", "fisher.pvalues.support()", "generate.pvalues()")
  
  # 'counts' must be a non-empty data frame or a non-empty matrix 
  qassert(x = counts, rules = c("M+", "D+"))
  # convert to data frame, if it is a matrix
  if(is.matrix(counts))
    counts <- as.data.frame(counts)
  # check if it contains only integer-like values and between two and four columns
  assert_data_frame(
    x = counts,
    types = "integerish",
    any.missing = FALSE,
    min.rows = 1,
    min.cols = 2,
    max.cols = 4
  )
  # round to integer
  counts <- round(counts)
  # ensure that the data frame has exactly two or exactly four columns
  qassert(x = counts, rules = c("D2", "D4"))
  
  # check and match 'input' parameter
  qassert(x = input, rules = "S1")
  input <- match.arg(input, c("noassoc", "marginal", "HG2011"))
  
  # transform 'counts' according to 'input'
  counts <- switch(
    EXPR = input,
    HG2011 = reconstruct_two(counts),
    marginal = reconstruct_four(counts, c(2, 4)),
    counts
  )
  
  # compute p-values and supports
  res <- fisher.test.pv(counts, alternative)
  
  # return list of results
  return(
    list(
      raw = res$get_pvalues(),
      support = res$get_pvalue_supports(unique = FALSE)
    )
  )
}

#' @title 
#' Generation of P-Values and Their Supports After Data Transformations
#' 
#' @description
#' Simple wrapper for generating p-values of discrete tests and their supports
#' after pre-processing the input data. The user only has to provide 1.) a
#' function that generates p-values and supports and 2.) an optional function
#' that pre-processes (i.e. transforms) the input data (if necessary) before it
#' can be used for p-value calculations. The respective arguments are provided
#' 
#' @templateVar dat TRUE
#' @templateVar test.fun TRUE
#' @templateVar test.args TRUE
#' @templateVar preprocess.fun TRUE
#' @templateVar preprocess.args TRUE
#' @template param
#'
#' @return
#' A [DiscreteTestResults][DiscreteTests::DiscreteTestResults] R6 class object.
#' 
#' @template exampleGPV
#' @examples
#' # Compute p-values and their supports of Fisher's exact test with pre-processing
#' df2 <- data.frame(X1, N1, X2, N2)
#' generate.pvalues(
#'   dat = df2,
#'   test.fun = "fisher_test_pv",
#'   preprocess.fun = function(tab) {
#'     for(col in c(2, 4)) tab[, col] <- tab[, col] - tab[, col - 1]
#'     return(tab)
#'   }
#' )
#' 
#' # Compute p-values and their supports of a binomial test with pre-processing
#' generate.pvalues(
#'   dat = rbind(c(5, 2, 7), c(3, 4, 0)), 
#'   test.fun = "binom_test_pv",
#'   test.args = list(n = c(9, 8, 11), p = 0.6, alternative = "two.sided"),
#'   preprocess.fun = colSums
#' )
#' 
#' @importFrom checkmate assert assert_choice assert_list check_function check_string
#' @import DiscreteTests
#' @export
generate.pvalues <- function(
  dat,
  test.fun,
  test.args = NULL,
  preprocess.fun = NULL,
  preprocess.args = NULL
) {
  #  make sure test function originates from package 'DiscreteTests'
  ## make sure it is a function or a string
  assert(
    check_function(test.fun),
    check_string(test.fun)
  )
  ## get all functions from package 'DiscreteTests'
  funs <- ls(asNamespace("DiscreteTests"))
  ## extract available test functions
  funs <- funs[endsWith(funs, "_test_pv")]
  ## make sure input 'test.fun' matches an available test function
  if(is.character(test.fun)) {
    ### match input 'test.fun' string to available test functions
    test.fun <- match.arg(tolower(test.fun), funs)
    ### convert string to actual function
    test.fun <- eval(parse(text = paste0("DiscreteTests::", test.fun)))
  } else {
    ### make sure input 'test.fun' matches an available test function
    OK <- FALSE
    for(fun in funs) {
      pkg_fun <- eval(parse(text = paste0("DiscreteTests::", fun)))
      if(all(all.equal(test.fun, pkg_fun) == TRUE)) {
        OK <- TRUE
        break;
      }
    }
    if(!OK) stop(paste("'test.fun' must be one of the '*_test_pv' functions of",
                       "package 'DiscreteTests'."))
  }
  
  # make sure parameters for 'test.fun' are in a named list or NULL
  assert_list(test.args, names = "unique", null.ok = TRUE)
  
  # make sure date pre-processing function is a function, string or NULL
  assert(
    check_function(preprocess.fun, null.ok = TRUE),
    check_string(preprocess.fun)
  )
  ## convert string to actual function (if necessary)
  if(is.character(preprocess.fun))
    preprocess.fun <- eval(parse(text = preprocess.fun))
  
  # make sure parameters for 'preprocess.fun' are in a named list or NULL
  assert_list(preprocess.args, names = "unique", null.ok = TRUE)
  
  if(!is.null(preprocess.fun)){
    # prepend data to pre-processing function's arguments list
    preprocess.args <- c(list(dat), preprocess.args)
    # set data parameter name according to first parameter of 'preprocess.fun'
    names(preprocess.args)[1] <- names(as.list(args(preprocess.fun)))[1]
    # perform pre-processing
    dat <- do.call(preprocess.fun, preprocess.args)
  }
  
  # get original data name
  data.name <- sapply(match.call(), deparse1)["dat"]
  # assign data to new variable with original name
  assign(data.name, dat)
  
  # prepend pre-processed data to test function's arguments list
  test.args <- c(list(as.name(data.name)), test.args)
  # set data parameter name according to first parameter of 'test.fun'
  names(test.args)[1] <- names(as.list(args(test.fun)))[1]
  # perform test(s)
  res <- do.call(test.fun, test.args)
  
  # output results
  return(res)
}
