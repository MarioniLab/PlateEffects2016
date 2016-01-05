set.seed(0)
require(limma)

# Batch-blocked:
for (x in 1:3) {
    if (x==1L) { 
        treatment <- factor(rep(1:2, each=30))
        batch <- factor(rep(1:6, each=10))
    } else if (x==2L) { 
        treatment <- factor(rep(1:2, each=150))
        batch <- factor(rep(1:6, each=50))
    } else {
        treatment <- factor(rep(1:2, each=100))
        batch <- factor(rep(1:20, each=10))
    }

    gvar <- 1
    design <- model.matrix(~treatment)

    for (sd in c(1:3)) { 
        # Data generation with a batch effect.
        dummy <- matrix(rnorm(10000 * length(treatment), sd=sqrt(gvar)), ncol=length(batch))
        batch.effect <- matrix(rnorm(10000 * 20, sd=sd), ncol=20)
        dummy <- dummy+batch.effect[,as.integer(batch)]
        truth <- (sd^2)/(sd^2+gvar)

        cor <- duplicateCorrelation(dummy, design, block=batch)
        fit <- lmFit(dummy, design, block=batch, correlation=cor$consensus)
        fit <- eBayes(fit, robust=TRUE)
        res <- topTable(fit, coef=ncol(design), sort.by="none", n=Inf)

        # Using the true value.
        alt.fit <- lmFit(dummy, design, block=batch, correlation=truth)
        alt.fit <- eBayes(alt.fit, robust=TRUE)
        alt.res <- topTable(alt.fit, coef=ncol(design), sort.by="none", n=Inf)

        cat(sprintf("Scenario is %s, sd is %i", x, sd), "\n")
        cat(sprintf("Consensus is %.5f, truth is %.5f", cor$consensus, truth), "\n")
        cat(sprintf("Type I error at 1%% for consensus is %.5f", sum(res$P.Value<=0.01)/10000), "\n")
        cat(sprintf("Type I error at 1%% for truth is %.5f", sum(alt.res$P.Value<=0.01)/10000), "\n")
    }
}

