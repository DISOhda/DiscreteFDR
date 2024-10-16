#' @return
#' A \code{DiscreteFDR} S3 class object whose elements are:
#' \item{Rejected}{rejected raw \eqn{p}-values.}
#' \item{Indices}{indices of rejected hypotheses.}
#' \item{Num.rejected}{number of rejections.}
#' \item{Adjusted}{adjusted \eqn{p}-values<%=ifelse(exists("BR") && BR, ""," (only for step-down direction)")%>.}
#' \item{Critical.constants}{critical values (only exists if computations where performed with `ret.crit.consts = TRUE`).}
#' \item{Select}{list with data related to \eqn{p}-value selection; only exists if `select.threshold < 1`.}
#' \item{Select$Threshold}{\eqn{p}-value selection threshold (`select.threshold`).}
#' \item{Select$Effective.Thresholds}{results of each \eqn{p}-value CDF evaluated at the selection threshold.}
#' \item{Select$Pvalues}{selected \eqn{p}-values that are \eqn{\leq} selection threshold.}
#' \item{Select$Indices}{indices of \eqn{p}-values \eqn{\leq} selection threshold.}
#' \item{Select$Scaled}{scaled selected \eqn{p}-values.}
#' \item{Select$Number}{number of selected \eqn{p}-values \eqn{\leq} selection threshold.}
#' \item{Data}{list with input data.}
#' \item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
#' \item{Data$Raw.pvalues}{observed \eqn{p}-values.}
#' \item{Data$pCDFlist}{list of the \eqn{p}-value supports.}
#' \item{Data$FDR.level}{FDR level `alpha`.}
#' \item{Data$Data.name}{the respective variable names of the input data.}
#' <%=ifelse(exists("BR") && BR, "\\item{DBR.Tuning}{value of the tuning parameter `lambda`.}","") %> 
