% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DBR_fun.R
\name{DBR}
\alias{DBR}
\title{[HBR-\eqn{\lambda}] procedure}
\usage{
DBR(
  raw.pvalues,
  pCDFlist,
  alpha = 0.05,
  lambda = NULL,
  ret.crit.consts = FALSE
)
}
\arguments{
\item{raw.pvalues}{vector of the raw observed p-values, as provided by the end user and before matching with their nearest neighbor in the CDFs supports.}

\item{pCDFlist}{a list of the supports of the CDFs of the p-values. Each support is represented by a vector that must be in increasing order.}

\item{alpha}{the target FDR level, a number strictly between 0 and 1.}

\item{lambda}{a number strictly between 0 and 1. If \code{lambda=NULL} (by default), then \code{lambda} is chosen equal to \code{alpha}.}

\item{ret.crit.consts}{a boolean. If \code{TRUE}, critical constants are computed and returned (this is computationally intensive).}
}
\value{
A \code{DiscreteFDR} S3 class object whose elements are:
\item{Rejected}{Rejected raw p-values}
\item{Indices}{Indices of rejected hypotheses}
\item{Num.rejected}{Number of rejections}
\item{Adjusted}{Adjusted p-values }
\item{Critical.constants}{Critical constants (if requested)}
\item{Method}{Character string describing the used algorithm, e.g. 'Discrete Benjamini-Hochberg procedure (step-up)'}
\item{Signif.level}{Significance level \code{alpha}}
\item{Lambda}{Value of \code{lambda}.}
\item{Data$raw.pvalues}{The values of \code{raw.pvalues}}
\item{Data$pCDFlist}{The values of \code{pCDFlist}}
\item{Data$data.name}{The respective variable names of \code{raw.pvalues} and \code{pCDFlist}}
}
\description{
Apply the [HBR-\eqn{\lambda}] procedure, with or without computing the
critical constants, to a set of p-values and their discrete support.
}
\details{
[DBR-\eqn{\lambda}] is the discrete version of the
[Blanchard-Roquain-\eqn{\lambda}] procedure (see References). The authors
of the latter suggest to take \code{lambda = alpha} (see their Proposition 17),
which explains the choice of the default value here.
}
\section{References}{

G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control
under independence and dependence. \emph{Journal of Machine Learning Research},
\emph{10}, pp. 2837-2871.
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

#Construction of the p-values and their support
df.formatted <- fisher.pvalues.support(counts = df, input = "noassoc")
raw.pvalues <- df.formatted$raw
pCDFlist <- df.formatted$support

DBR.fast <- DBR(raw.pvalues, pCDFlist)
summary(DBR.fast)
DBR.crit <- DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
summary(DBR.crit)

}
\seealso{
\link{discrete.BH}
}
