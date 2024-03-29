% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_funs.R
\name{print.DiscreteFDR}
\alias{print.DiscreteFDR}
\title{Printing DiscreteFDR results}
\usage{
\method{print}{DiscreteFDR}(x, ...)
}
\arguments{
\item{x}{an object of class "\code{DiscreteFDR}".}

\item{...}{further arguments to be passed to or from other methods.
They are ignored in this function.}
}
\value{
The respective input object is invisibly returned via \code{invisible(x)}.
}
\description{
Prints the results of discrete FDR analysis, stored in a \code{DiscreteFDR}
S3 class object.
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

DBH.su.crit <- DBH(raw.pvalues, pCDFlist, direction = "su",
                   ret.crit.consts = TRUE)
print(DBH.su.crit)

}
