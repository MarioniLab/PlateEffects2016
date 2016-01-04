# This checks the behaviour of limma, limma+cor and limma/sum on normally distributed data.
# The idea is to show that we still get roughly the same performance, so the main cause of
# the decrease in power is probably due to count discreteness. Similar behaviour is observed
# in one-way layouts and with an additive model.

set.seed(0)

# Data generation with a DE gene & plate effect.
plate <- factor(rep(1:6, each=50))
treatment <- factor(rep(1:2, each=150))
ngenes <- 10000
nde <- 1000
is.de <- seq_len(nde)

gvar <- 0.5
dummy <- matrix(rnorm(ngenes * length(treatment), sd=sqrt(gvar)), ncol=length(plate))
dummy[is.de,] <- t(t(dummy[is.de,]) + c(-1, 1)[treatment])
plate.effect <- matrix(rnorm(ngenes * 20, sd=1), ncol=20)
dummy <- dummy+plate.effect[,as.integer(plate)]

# Testing what happens with and without a blocking model.

for (type in c("oneway", "blocked")) { 
if (type=="oneway") { 
    design <- model.matrix(~treatment)
    refdesign <- model.matrix(~factor(rep(1:2, each=3)))
} else {
    batch <- factor(rep(rep(1:3, 2), each=50))
    design <- model.matrix(~batch + treatment)
    refdesign <- model.matrix(~factor(rep(1:3, 2)) + factor(rep(1:2, each=3)))
}

# Analyzing with various methods.

require(limma)
fit <- lmFit(dummy, design)
fit <- eBayes(fit, robust=TRUE)
res <- topTable(fit, coef=ncol(design), sort.by="none", n=Inf)

cor <- duplicateCorrelation(dummy, design, block=plate)
fit.2 <- lmFit(dummy, design, block=plate, correlation=cor$consensus)
fit.2 <- eBayes(fit.2, robust=TRUE)
res.2 <- topTable(fit.2, coef=ncol(design), sort.by="none", n=Inf)

reffit <- lmFit(avearrays(dummy, plate), refdesign)
reffit <- eBayes(reffit, robust=TRUE)
refres <- topTable(reffit, coef=ncol(refdesign), sort.by="none", n=Inf)

# Computing FP and TP rates.

hits <- 10^((-50):0)
tp <- sapply(hits, function(x) { sum(res$P.Value[is.de] <= x) })/nde
fp <- sapply(hits, function(x) { sum(res$P.Value[-is.de] <= x) })/(ngenes-nde)
tp2 <- sapply(hits, function(x) { sum(res.2$P.Value[is.de] <= x) })/nde
fp2 <- sapply(hits, function(x) { sum(res.2$P.Value[-is.de] <= x) })/(ngenes-nde)
tpr <- sapply(hits, function(x) { sum(refres$P.Value[is.de] <= x) })/nde
fpr <- sapply(hits, function(x) { sum(refres$P.Value[-is.de] <= x) })/(ngenes-nde)

thresholds <- c(0.005, 0.01, 0.05)
pch <- c(16, 17, 18)
xtp <-  sapply(thresholds, function(x) { sum(res$P.Value[is.de] <= x) })/nde
xfp <-  sapply(thresholds, function(x) { sum(res$P.Value[-is.de] <= x) })/(ngenes-nde)
xtp2 <- sapply(thresholds, function(x) { sum(res.2$P.Value[is.de] <= x) })/nde
xfp2 <- sapply(thresholds, function(x) { sum(res.2$P.Value[-is.de] <= x) })/(ngenes-nde)
xtpr <- sapply(thresholds, function(x) { sum(refres$P.Value[is.de] <= x) })/nde
xfpr <- sapply(thresholds, function(x) { sum(refres$P.Value[-is.de] <= x) })/(ngenes-nde)

cat("Current type is", type, "\n")
cat("Default voom:\n")
print(data.frame(Threshold=thresholds, TP=xtp, FP=xfp))
cat("voom + correlations:\n")
print(data.frame(Threshold=thresholds, TP=xtp2, FP=xfp2))
cat("voom on averaged arrays:\n")
print(data.frame(Threshold=thresholds, TP=xtpr, FP=xfpr))
cat("\n")

# Making a plot.

setEPS()
postscript(paste0(type, ".eps"), width=10, height=6)
par(mfrow=c(1, 2))
plot(fp, tp, xlim=c(0, 1), ylim=c(0, 1), main=type)
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")
points(xfp, xtp, col="black", pch=pch)
points(xfp2, xtp2, col="grey50", pch=pch)
points(xfpr, xtpr, col="red", pch=pch)

plot(fp, tp, xlim=c(0, 0.1), ylim=c(0, 1), main=type)
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")
points(xfp, xtp, col="black", pch=pch)
points(xfp2, xtp2, col="grey50", pch=pch)
points(xfpr, xtpr, col="red", pch=pch)
abline(v=thresholds, col="dodgerblue", lty=2)
dev.off()
}

# Session information

sessionInfo()
