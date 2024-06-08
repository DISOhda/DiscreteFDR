## ----format-am----------------------------------------------------------------
data(amnesia)
amnesia.formatted <- fisher.pvalues.support(amnesia[, 2:3], input = "HG2011")
raw.pvalues <- amnesia.formatted$raw
pCDFlist <- amnesia.formatted$support

