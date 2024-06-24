#' @return
#' A \code{DiscreteFDR} S3 class object whose elements are:
#' \item{Rejected}{rejected raw p-values.}
#' \item{Indices}{indices of rejected hypotheses.}
#' \item{Num.rejected}{number of rejections.}
#' \item{Adjusted}{adjusted p-values<%=ifelse(exists("BR") && BR, ""," (only for step-down direction)")%>.}
#' \item{Critical.constants}{critical constants (if requested).}
#' \item{Select$Threshold}{p-value selection \code{threshold} (only present if \code{threshold} \eqn{< 1}).}
#' \item{Select$Effective.Thresholds}{results of each p-value CDF evaluated at the selection threshold (only present if \code{threshold} \eqn{< 1}).}
#' \item{Select$Pvalues}{selected p-values that are \eqn{\leq} selection \code{threshold} (only present if \code{threshold} \eqn{< 1}).}
#' \item{Select$Indices}{indices of p-values \eqn{\leq} selection \code{threshold} (only present if \code{threshold} \eqn{< 1}).}
#' \item{Select$Scaled}{scaled selected p-values, i.e. \eqn{\frac{i-th p-value}{i-th effective threshold}} (only present if \code{threshold} \eqn{< 1}).}
#' \item{Select$Number}{number of selected p-values \eqn{\leq} \code{threshold} (only present if \code{threshold} \eqn{< 1}).}
#' \item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
#' \item{Data$raw.pvalues}{observed p-values.}
#' \item{Data$pCDFlist}{p-value supports.}
#' \item{Data$pCDFlist.indices}{observed p-value indices of the **unique** p-value supports.}
#' \item{Data$FDR.level}{FDR level \code{alpha}.}
#' \item{Data$Data.name}{the respective variable names of the input data.}
#' <%=ifelse(exists("BR") && BR, "\\item{DBR.Tuning}{value of the tuning parameter \\code{lambda}.}","") %> 
