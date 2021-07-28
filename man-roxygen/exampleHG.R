#'@examples
#'
#'data(amnesia)
#'#We only keep the first 100 lines to keep the computations fast.
#'#We also drop the first column to keep only columns of counts, in the Heller & Gur (2011) setting.
#'amnesia <- amnesia[1:100,2:3]
#'
#'#Construction of the p-values and their support
#'amnesia.formatted <- fisher.pvalues.support(counts = amnesia, input = "HG2011")
#'raw.pvalues <- amnesia.formatted$raw
#'pCDFlist <- amnesia.formatted$support
