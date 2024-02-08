## ----toy-example-data, results='asis'-----------------------------------------
library(knitr)
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
kable(df, caption = "Toy Example")

## ----toy-example-5------------------------------------------------------------
library("DiscreteFDRTest")
DBH.sd.fast <- fast.Discrete(df, alternative = "two.sided", direction = "sd")
print(DBH.sd.fast)

summary(DBH.sd.fast)

## ----toy-example-summary------------------------------------------------------
DBH.sd.fast.summary <- summary(DBH.sd.fast)
DBH.sd.fast.summary$Table

## ----tow-example-2------------------------------------------------------------
p <- fisher.pvalues.support(df, alternative = "two.sided")
raw.pvalues <- p$raw

# or:
raw.pvalues2 <- DBH.sd.fast$Data$raw.pvalues

all(raw.pvalues == raw.pvalues2)

p.adjust(raw.pvalues, method = "BH")

## ----toy-example-crit---------------------------------------------------------
p <- fisher.pvalues.support(df, alternative = "two.sided")
raw.pvalues <- p$raw
pCDFlist <- p$support

DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, 0.05, "sd", TRUE)
crit.vals.BH.disc <- DBH.sd.crit$Critical.values
crit.vals.BH.cont <- 1:9 * 0.05/9
cbind(sort(raw.pvalues), crit.vals.BH.disc, crit.vals.BH.cont)

## ----toy-example-plot-1-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     legend = "topleft", cex = 1.3)

## ----toy-example-plot-2-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     cex = 1.3, ylim = c(0, 0.25), main = "Comparison of discrete and continuous BH procedures")
points(crit.vals.BH.cont, pch = 19, cex = 1.3, lwd = 2)
legend("topright", c("Rejected", "Accepted", "Critical Values (disc.)", "Critical Values (cont.)"),
       col = c("red", "blue", "green", "black"), pch = c(4, 2, 19, 19), lwd = 2, lty = 0)

## ----toy-example-3------------------------------------------------------------
p$support[c(1,5)]

## ----toy-example-4------------------------------------------------------------
pCDFlist <- p$support
stepf <- lapply(pCDFlist, function(x) stepfun(x, c(0, x)))
par(mfcol = c(1, 3), mai = c(1, 0.5, 0.3, 0.1))
plot(stepf[[1]], xlim = c(0, 1), ylim = c(0, 1), do.points = FALSE, lwd = 1, lty = 1, ylab = "F(x)", 
     main = "(a)")
for(i in (2:9)){
  plot(stepf[[i]], add = TRUE, do.points = FALSE, lwd = 1, col = i)
}
segments(0, 0, 1, 1, col = "grey", lty = 2)

#   Plot xi
support <- sort(unique(unlist(pCDFlist)))
components <- lapply(stepf, function(s){s(support) / (1 - s(support))}) 
xi.values <- 1/9 * Reduce('+', components)
xi <- stepfun(support, c(0, xi.values))
plot(xi, xlim = c(0, 0.10), ylim = c(0, 0.10), do.points = FALSE, ylab = expression(xi), main = "(b)")
segments(0, 0, 0.1, 0.1, col = "grey", lty = 2)

#   Plot discrete critical values as well a BH constants and raw p-values
DBH.sd <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
plot(DBH.sd, col = c("black", "black", "red"), pch = c(4, 4, 19), type.crit = 'p', ylim = c(0, 0.15),
     cex = 1.3, main = "(c)", ylab = "Critical Values")
points(1:9, 0.05 * (1:9) / 9, col = "green", pch = 19, cex = 1.3)

mtext("Figure 1", 1, outer = TRUE, line = -2)

