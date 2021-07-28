#'@return
#'A \code{DiscreteFDR} S3 class object whose elements are:
#'\item{Rejected}{Rejected raw p-values}
#'\item{Indices}{Indices of rejected hypotheses}
#'\item{Num.rejected}{Number of rejections}
#'\item{Adjusted}{Adjusted p-values <%=ifelse(exists("DBR") && DBR, "","(only for step-down direction).") %>}
#'\item{Critical.constants}{Critical constants (if requested)}
#'\item{Method}{Character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
#'\item{Signif.level}{Significance level \code{alpha}}
#'<%=ifelse(exists("DBR") && DBR, "\\item{Lambda}{Value of \\code{lambda}.}","") %>
#'\item{Data$raw.pvalues}{The values of \code{raw.pvalues}}
#'\item{Data$pCDFlist}{The values of \code{pCDFlist}}
#'\item{Data$data.name}{The respective variable names of \code{raw.pvalues} and \code{pCDFlist}}
