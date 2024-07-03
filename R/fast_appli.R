#' @title Fast Application of Discrete Multiple Testing Procedures
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Apply the \[HSU\], \[HSD\], \[AHSU\] or \[AHSD\] procedure,
#' without computing the critical constants, to a data set of 2x2 contingency
#' tables which may have to be pre-processed in order to have the correct
#' structure for computing p-values using Fisher's exact test.
#' 
#' **Note**: This function is deprecated and will be removed in a future
#' version. Please use [`direct.discrete.BH()`] with
#' `test.fun = DiscreteTests::fisher.test.pv` and (optional) 
#' `preprocess.fun = DiscreteDatasets::reconstruct_two` or 
#' `preprocess.fun = DiscreteDatasets::reconstruct_four` instead. Alternatively,
#' use a pipeline, e.g.\cr
#' `data |>`\cr
#' `  DiscreteDatasets::reconstruct_*(<args>) |>`\cr
#' `  DiscreteTests::*.test.pv(<args>) |>`\cr
#' `  discrete.BH(<args>)`.
#'
#' @param counts        a data frame of two or four columns and any number of
#'                      lines; each line representing a 2x2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument (see Details
#'                      section of [`fisher.pvalues.support()`]).
#' @param alternative   same argument as in [`stats::fisher.test()`]. The three
#'                      possible values are `"greater"` (default), `"two.sided"`
#'                      or `"less"` (may be abbreviated).
#' @param input         the format of the input data frame (see Details section
#'                      of [`fisher.pvalues.support()`]. The three possible
#'                      values are `"noassoc"` (default), `"marginal"` or
#'                      `"HG2011"` (may be abbreviated).
#'
#' @template param
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
#' @templateVar select.threshold TRUE
#' 
#' @template return
#' 
#' @seealso
#' [`fisher.pvalues.support()`], [`discrete.BH()`]
#' 
#' @template example
#' @examples
#' DBH.su <- fast.Discrete(df, input = "noassoc", direction = "su")
#' summary(DBH.su)
#' 
#' DBH.sd <- fast.Discrete(df, input = "noassoc", direction = "sd")
#' DBH.sd$Adjusted
#' summary(DBH.sd)
#' 
#' ADBH.su <- fast.Discrete(df, input = "noassoc", direction = "su", adaptive = TRUE)
#' summary(ADBH.su)
#' 
#' ADBH.sd <- fast.Discrete(df, input = "noassoc", direction = "sd", adaptive = TRUE)
#' ADBH.sd$Adjusted
#' summary(ADBH.sd)
#' 
#' @importFrom lifecycle deprecate_soft
#' @export
fast.Discrete <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, direction = "su", adaptive = FALSE, select.threshold = 1){
  deprecate_soft("1.3.7", "fast.Discrete()", "direct.discrete.BH()")
  
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive, FALSE, select.threshold)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}


#' @title 
#' Direct Application of Multiple Testing Procedures to Dataset
#' 
#' @description
#' Apply the \[HSU\], \[HSD\], \[AHSU\] or \[AHSD\] procedure, with or without
#' computing the critical constants,
#' to a data set of 2x2 contingency tables using Fisher's exact tests which
#' may have to be transformed before computing p-values.
#' 
#' @templateVar dat TRUE
#' @templateVar test.fun TRUE
#' @templateVar test.args TRUE
#' @templateVar alpha TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
#' @templateVar ret.crit.consts TRUE
#' @templateVar select.threshold TRUE
#' @templateVar preprocess.fun TRUE
#' @templateVar preprocess.args TRUE
#' @template param
#' 
#' @template example
#' @examples
#' DBH.su <- direct.discrete.BH(df, "fisher", direction = "su")
#' summary(DBH.su)
#' 
#' DBH.sd <- direct.discrete.BH(df, "fisher", direction = "sd")
#' DBH.sd$Adjusted
#' summary(DBH.sd)
#' 
#' ADBH.su <- direct.discrete.BH(df, "fisher", direction = "su", adaptive = TRUE)
#' summary(ADBH.su)
#' 
#' ADBH.sd <- direct.discrete.BH(df, "fisher", direction = "sd", adaptive = TRUE)
#' ADBH.sd$Adjusted
#' summary(ADBH.sd)
#' 
#' @export
direct.discrete.BH <- function(
  dat,
  test.fun, 
  test.args = NULL,
  alpha = 0.05, 
  direction = "su",
  adaptive = FALSE,
  ret.crit.consts = FALSE,
  select.threshold = 1,
  preprocess.fun = NULL, 
  preprocess.args = NULL
) {
  out <- discrete.BH.DiscreteTestResults(
    test.results = generate.pvalues(
      dat             = dat,
      test.fun        = test.fun,
      test.args       = test.args,
      preprocess.fun  = preprocess.fun,
      preprocess.args = preprocess.args
    ),
    alpha            = alpha,
    direction        = direction,
    adaptive         = adaptive,
    ret.crit.consts  = ret.crit.consts,
    select.threshold = select.threshold
  )
  
  out$Data$Data.name <- deparse(substitute(dat))
  
  return(out)
}
