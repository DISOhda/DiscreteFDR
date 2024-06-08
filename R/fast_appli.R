#' @title Fast application of discrete procedures
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Apply the \[HSU\], \[HSD\], \[AHSU\] or \[AHSD\] procedure,
#' without computing the critical constants,
#' to a data set of 2x2 contingency tables using Fisher's exact tests.
#' 
#' **Note**: In future versions, this function will be removed and replaced by a
#'  more flexible one, which will not be limited to Fisher's exact test.
#'
#' @param counts        a data frame of two or four columns and any number of
#'                      lines; each line representing a 2x2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument (see Details
#'                      section of [fisher.pvalues.support]).
#' @param alternative   same argument as in [fisher.test]. The three
#'                      possible values are `"greater"` (default), `"two.sided"`
#'                      or `"less"` (may be abbreviated).
#' @param input         the format of the input data frame (see Details section
#'                      of [fisher.pvalues.support]. The three possible values
#'                      are `"noassoc"` (default), `"marginal"` or `"HG2011"`
#'                      (may be abbreviated).
#'
#' @templateVar stepf FALSE
#' @templateVar pv.numer FALSE
#' @templateVar pv.denom FALSE
#' @templateVar alpha TRUE
#' @templateVar sorted.pv FALSE
#' @templateVar pCDFlist FALSE
#' @templateVar raw.pvalues FALSE
#' @templateVar direction TRUE
#' @templateVar ret.crit.consts FALSE
#' @templateVar sorted.num FALSE
#' @templateVar t FALSE
#' @templateVar lambda FALSE
#' @templateVar adaptive TRUE
#' @templateVar threshold TRUE
#' @template param
#' 
#' @return 
#' A `DiscreteFDR` S3 class object whose elements are:
#' \item{Rejected}{Rejected raw p-values}
#' \item{Indices}{Indices of rejected hypotheses}
#' \item{Num.rejected}{Number of rejections}
#' \item{Adjusted}{Adjusted p-values (only for step-down direction).}
#' \item{Method}{Character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
#' \item{Signif.level}{Significance level `alpha`}
#' \item{Data$raw.pvalues}{The values of `raw.pvalues`}
#' \item{Data$pCDFlist}{The values of `pCDFlist`}
#' \item{Data$data.name}{The variable name of the `counts` dataset}
#' 
#' @seealso
#' [fisher.pvalues.support], [discrete.BH]
#' 
#' @examples 
#' X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#' N1 <- rep(148, 9)
#' N2 <- rep(132, 9)
#' Y1 <- N1 - X1
#' Y2 <- N2 - X2
#' df <- data.frame(X1, Y1, X2, Y2)
#' df
#' 
#' DBH.su <- fast.Discrete(counts = df, input = "noassoc", direction = "su")
#' summary(DBH.su)
#' 
#' DBH.sd <- fast.Discrete(counts = df, input = "noassoc", direction = "sd")
#' DBH.sd$Adjusted
#' summary(DBH.sd)
#' 
#' ADBH.su <- fast.Discrete(counts = df, input = "noassoc", direction = "su", adaptive = TRUE)
#' summary(ADBH.su)
#' 
#' ADBH.sd <- fast.Discrete(counts = df, input = "noassoc", direction = "sd", adaptive = TRUE)
#' ADBH.sd$Adjusted
#' summary(ADBH.sd)
#' 
#' @name fast.Discrete
#' @importFrom lifecycle deprecate_soft
#' @export
fast.Discrete <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, direction = "su", adaptive = FALSE, threshold = 1){
  deprecate_soft("1.3.7", "fast.Discrete()",
                 details = paste("In future versions of this package, this",
                                 "function will be removed and replaced by a",
                                 "more flexible one, which will not be",
                                 "limited to Fisher's exact test."))
  
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive, FALSE, threshold)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}
