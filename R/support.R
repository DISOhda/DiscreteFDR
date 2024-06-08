#' @title
#' Computing discrete p-values and their support for binomial and Fisher's
#' exact tests
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Computes discrete raw p-values and their support for Fisher's exact test
#' applied to 2x2 contingency tables summarizing counts coming from two
#' categorical measurements.
#'
#' **Note**: This function is deprecated. Please use
#' [`DiscreteTests::fisher.test.pv()`] after 
#' [`DiscreteDatasets::reconstruct_two()`] or
#' [`DiscreteDatasets::reconstruct_four()`].
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
#' If `input="HG2011"`, we are in the situation of the [amnesia] data set as
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
#' [fisher.test]
#' 
#' @param counts        a data frame of two or four columns and any number of
#'                      lines; each line represents a 2x2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument, see
#'                      Details.
#' @param alternative   same argument as in [stats::fisher.test]. The three
#'                      possible values are `"greater"` (default),
#'                      `"two.sided"` or `"less"` and you can specify
#'                      just the initial letter.
#' @param input         the format of the input data frame, see Details. The
#'                      three possible values are `"noassoc"` (default),
#'                      `"marginal"` or `"HG2011"` and you can specify
#'                      just the initial letter.
#' 
#' @template example
#' @template exampleHG
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
#' @importFrom DiscreteTests fisher.test.pv
#' @export
fisher.pvalues.support <- function(counts, alternative = "greater", input = "noassoc"){
  deprecate_soft("1.3.7", "fisher.pvalues.support()",
                 details = paste("Please use",
                                 "`DiscreteDatasets::reconstruct_*() %>%",
                                 "DiscreteTests::fisher.test.pv()` or just",
                                 "`DiscreteTests::fisher.test.pv()` instead."))
  
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
    marginal = reconstruct_four(counts),
    counts
  )
  
  # compute p-values and supports
  res <- fisher.test.pv(counts, alternative)
  
  # return list of results
  return(list(raw = res$get_pvalues(),
              support = res$get_pvalue_supports(unique = FALSE)))
}

#' @importFrom checkmate assert_function assert_list test_r6
#' @export
generate.pvalues <- function(
  data,
  test_fun,
  test_args = NULL,
  preprocess_fun = NULL,
  preprocess_args = NULL
) {
  assert_function(test_fun)
  assert_list(test_args, names = "unique", null.ok = TRUE)
  assert_function(preprocess_fun, null.ok = TRUE)
  assert_list(preprocess_args, names = "unique", null.ok = TRUE)
  
  if(!is.null(preprocess_fun)){
    xname <- names(as.list(args(preprocess_fun)))[1]
    preprocess_args <- c(list(x = data), preprocess_args)
    names(preprocess_args)[1] <- xname
    
    data <- do.call(preprocess_fun, preprocess_args)
  }
  
  xname <- names(as.list(args(test_fun)))[1]
  test_args <- c(list(x = data), test_args)
  names(test_args)[1] <- xname
  
  res <- do.call(test_fun, test_args)
  if(!test_r6(res, "DiscreteTestResults"))
    warning(
      paste("Result must be an R6 object of class 'DiscreteTestResults'.",
            "Use a test function from package 'DiscreteTests' as 'test_fun'."#,
            #sep = "\n  "
      )
    )
  
  return(res)
}
