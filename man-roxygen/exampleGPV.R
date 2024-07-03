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
#' # Compute p-values and their supports of Fisher's exact test
#' test.result <- generate.pvalues(df, "fisher")
#' raw.pvalues <- test.result$get_pvalues()
#' pCDFlist <- test.result$get_pvalue_supports()
#' 
