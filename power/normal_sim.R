# This checks the behaviour of limma, limma+cor and limma/sum on normally distributed data.
# The idea is to show that we still get roughly the same performance.

set.seed(0)

# Data generation with a DE gene & plate effect.
plate <- factor(rep(1:6, each=50))
treatment <- factor(rep(1:2, each=150))

gvar <- 0.5
dummy <- matrix(rnorm(10000 * length(treatment), sd=sqrt(gvar)), ncol=length(plate))
dummy[1:1000,] <- t(t(dummy[1:1000,]) + c(-1, 1)[treatment])

plate.effect <- matrix(rnorm(10000 * 20, sd=1), ncol=20)
dummy <- dummy+plate.effect[,as.integer(plate)]

pdf("normal_results.pdf")
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

# Making a plot.

hits <- 10^((-50):0)
tp <- sapply(hits, function(x) { sum(res$P.Value[1:1000] <= x) })
fp <- sapply(hits, function(x) { sum(res$P.Value[-(1:1000)] <= x) })
tp2 <- sapply(hits, function(x) { sum(res.2$P.Value[1:1000] <= x) })
fp2 <- sapply(hits, function(x) { sum(res.2$P.Value[-(1:1000)] <= x) })
tpr <- sapply(hits, function(x) { sum(refres$P.Value[1:1000] <= x) })
fpr <- sapply(hits, function(x) { sum(refres$P.Value[-(1:1000)] <= x) })

par(mfrow=c(1, 2))
plot(fp, tp, xlim=c(0, 10000), ylim=c(0, 1000), main=type)
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")

plot(fp, tp, xlim=c(0, 1000), ylim=c(0, 1000), main=type)
points(fp2, tp2, col="grey50")
points(fpr, tpr, col="red")

}
dev.off()
