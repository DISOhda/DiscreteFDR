#' @return
#' A \code{DiscreteFDR} S3 class object whose elements are:
#' \item{Rejected}{rejected raw p-values.}
#' \item{Indices}{indices of rejected hypotheses.}
#' \item{Num.rejected}{number of rejections.}
#' \item{Adjusted}{adjusted p-values<%=ifelse(exists("BR") && BR, ""," (only for step-down direction)")%>.}
#' \item{Critical.constants}{critical values (only exists if computations where performed with `ret.crit.consts = TRUE`).}
#' \item{Select$Threshold}{p-value selection `threshold` (only exists if `threshold < 1`).}
#' \item{Select$Effective.Thresholds}{results of each p-value CDF evaluated at the selection threshold (only exists if `threshold < 1`).}
#' \item{Select$Pvalues}{selected p-values that are \eqn{\leq} selection `threshold` (only exists if `threshold < 1`).}
#' \item{Select$Indices}{indices of p-values \eqn{\leq} selection `threshold` (only exists if `threshold < 1`).}
#' \item{Select$Scaled}{scaled selected p-values (only exists if `threshold < 1`).}
#' \item{Select$Number}{number of selected p-values \eqn{\leq} `threshold` (only exists if `threshold < 1`).}
#' \item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
#' \item{Data$raw.pvalues}{observed p-values.}
#' \item{Data$pCDFlist}{list of the p-value supports.}
#' \item{Data$FDR.level}{FDR level `alpha`.}
#' \item{Data$Data.name}{the respective variable names of the input data.}
#' <%=ifelse(exists("BR") && BR, "\\item{DBR.Tuning}{value of the tuning parameter `lambda`.}","") %> 
