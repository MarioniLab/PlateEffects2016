# This checks the behaviour of limma, limma+cor and limma/sum on normally distributed data.
# The idea is to show that we still get roughly the same performance.

set.seed(0)

# Data generation with a DE gene & batch effect.
batch <- factor(rep(1:6, each=50))
treatment <- factor(rep(1:2, each=150))

gvar <- 0.5
dummy <- matrix(rnorm(10000 * length(treatment), sd=sqrt(gvar)), ncol=length(batch))
dummy[1:1000,] <- t(t(dummy[1:1000,]) + c(-1, 1)[treatment])

batch.effect <- matrix(rnorm(10000 * 20, sd=1), ncol=20)
dummy <- dummy+batch.effect[,as.integer(batch)]

# Analyzing with various methods.

require(limma)
design <- model.matrix(~treatment)
fit <- lmFit(dummy, design)
fit <- eBayes(fit, robust=TRUE)
res <- topTable(fit, coef=2, sort.by="none", n=Inf)

cor <- duplicateCorrelation(dummy, design, block=batch)
fit.2 <- lmFit(dummy, design, block=batch, correlation=cor$consensus)
fit.2 <- eBayes(fit.2, robust=TRUE)
res.2 <- topTable(fit.2, coef=2, sort.by="none", n=Inf)

refdesign <- model.matrix(~c(0,0,0,1,1,1))
reffit <- lmFit(avearrays(dummy, batch), refdesign)
reffit <- eBayes(reffit, robust=TRUE)
refres <- topTable(reffit, coef=ncol(refdesign), sort.by="none", n=Inf)

# Making a plot.

hits <- 10^((-50):0)
tp <- sapply(hits, function(x) { sum(res$P.Value[1:1000] <= x) })
fp <- sapply(hits, function(x) { sum(res$P.Value[-(1:1000)] <= x) })
tp2 <- sapply(hits, function(x) { sum(res.2$P.Value[1:1000] <= x) })
fp2 <- sapply(hits, function(x) { sum(res.2$P.Value[-(1:1000)] <= x) })
tpr <- sapply(hits, function(x) { sum(refres$P.Value[1:1000] <= x) })
fpr <- sapply(hits, function(x) { sum(refres$P.Value[-(1:1000)] <= x) })

plot(fp, tp, xlim=c(0, 10000), ylim=c(0, 1000))
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")

plot(fp, tp, xlim=c(0, 1000), ylim=c(0, 1000))
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")

