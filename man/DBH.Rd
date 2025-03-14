% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DBH_fun.R
\name{DBH}
\alias{DBH}
\alias{DBH.default}
\alias{DBH.DiscreteTestResults}
\title{Wrapper Functions for the Discrete Benjamini-Hochberg Procedure}
\usage{
DBH(test.results, ...)

\method{DBH}{default}(
  test.results,
  pCDFlist,
  alpha = 0.05,
  direction = "su",
  ret.crit.consts = FALSE,
  select.threshold = 1,
  pCDFlist.indices = NULL,
  ...
)

\method{DBH}{DiscreteTestResults}(
  test.results,
  alpha = 0.05,
  direction = "su",
  ret.crit.consts = FALSE,
  select.threshold = 1,
  ...
)
}
\arguments{
\item{test.results}{either a numeric vector with \eqn{p}-values or an R6 object of class \code{\link[DiscreteTests]{DiscreteTestResults}} from package \link[DiscreteTests]{DiscreteTests} for which a discrete FDR procedure is to be performed.}

\item{...}{further arguments to be passed to or from other methods. They are ignored here.}

\item{pCDFlist}{list of the supports of the CDFs of the \eqn{p}-values; each list item must be a numeric vector, which is sorted in increasing order and whose last element equals 1.}

\item{alpha}{single real number strictly between 0 and 1 indicating the target FDR level.}

\item{direction}{single character string specifying whether to perform a step-up (\code{"su"}; the default) or step-down procedure (\code{"sd"}).}

\item{ret.crit.consts}{single boolean specifying whether critical constants are to be computed.}

\item{select.threshold}{single real number strictly between 0 and 1 indicating the largest raw \eqn{p}-value to be considered, i.e. only \eqn{p}-values below this threshold are considered and the procedures are adjusted in order to take this selection effect into account; if \code{threshold = 1} (the default), all raw \eqn{p}-values are selected.}

\item{pCDFlist.indices}{list of numeric vectors containing the test indices that indicate to which raw \eqn{p}-value each \strong{unique} support in \code{pCDFlist} belongs; ignored if the lengths of \code{test.results} and \code{pCDFlist} are equal.}
}
\value{
A \code{DiscreteFDR} S3 class object whose elements are:
\item{Rejected}{rejected raw \eqn{p}-values.}
\item{Indices}{indices of rejected hypotheses.}
\item{Num.rejected}{number of rejections.}
\item{Adjusted}{adjusted \eqn{p}-values (only for step-down direction).}
\item{Critical.constants}{critical values (only exists if computations where performed with \code{ret.crit.consts = TRUE}).}
\item{Data}{list with input data.}
\item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
\item{Data$Raw.pvalues}{observed \eqn{p}-values.}
\item{Data$pCDFlist}{list of the \eqn{p}-value supports.}
\item{Data$FDR.level}{FDR level \code{alpha}.}

\item{Data$Data.name}{the respective variable names of the input data.}
\item{Select}{list with data related to \eqn{p}-value selection; only exists if \code{select.threshold < 1}.}
\item{Select$Threshold}{\eqn{p}-value selection threshold (\code{select.threshold}).}
\item{Select$Effective.Thresholds}{results of each \eqn{p}-value CDF evaluated at the selection threshold.}
\item{Select$Pvalues}{selected \eqn{p}-values that are \eqn{\leq} selection threshold.}
\item{Select$Indices}{indices of \eqn{p}-values \eqn{\leq} selection threshold.}
\item{Select$Scaled}{scaled selected \eqn{p}-values.}
\item{Select$Number}{number of selected \eqn{p}-values \eqn{\leq} selection threshold.}
}
\description{
\code{DBH()} is a wrapper function of \code{\link[=discrete.BH]{discrete.BH()}} for computing [HSU] and
[HSD]. It simply passes its arguments to \code{\link[=discrete.BH]{discrete.BH()}} with fixed
\code{adaptive = FALSE}.
}
\details{
Computing critical constants (\code{ret.crit.consts = TRUE}) requires considerably
more execution time, especially if the number of unique supports is large.
We recommend that users should only have them calculated when they need them,
e.g. for illustrating the rejection set in a plot or other theoretical
reasons.
}
\examples{
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
df

# Compute p-values and their supports of Fisher's exact test
test.result <- generate.pvalues(df, "fisher")
raw.pvalues <- test.result$get_pvalues()
pCDFlist <- test.result$get_pvalue_supports()

# DBH (step-up) without critical values; using test results object
DBH.su.fast <- DBH(test.result)
summary(DBH.su.fast)

# DBH (step-down) without critical values; using extracted p-values 
# and supports
DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
summary(DBH.sd.fast)

# DBH (step-up) with critical values; using extracted p-values and supports
DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
summary(DBH.su.crit)

# DBH (step-down) with critical values; using test results object
DBH.sd.crit <- DBH(test.result, direction = "sd", ret.crit.consts = TRUE)
summary(DBH.sd.crit)

}
\references{
Döhler, S., Durand, G., & Roquain, E. (2018). New FDR bounds for discrete
and heterogeneous tests. \emph{Electronic Journal of Statistics}, \emph{12}(1),
pp. 1867-1900. \doi{10.1214/18-EJS1441}
}
\seealso{
\code{\link[=discrete.BH]{discrete.BH()}}, \code{\link[=ADBH]{ADBH()}}, \code{\link[=DBR]{DBR()}}, \code{\link[=DBY]{DBY()}}
}
