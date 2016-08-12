# This script aims to show that ComBat, despite its EB treatment of 
# the batch effect, is not a panacea for confounding designs.

plate.ids <- factor(rep(1:6, each=50))
conditions <- factor(rep(1:2, each=150))

require(limma)
design <- model.matrix(~conditions)

# We can't use it with the 'design', as it's confounding and breaks it.
require(sva)
x <- matrix(0, ncol=length(plate.ids), nrow=10)
try(y <- ComBat(x, plate.ids, mod=design)) # Confounding, breaks.
try(y <- ComBat(x, plate.ids, mod=design[,-1])) 

for (n in c(100, 200, 500)) { 
    x <- matrix(rnorm(600000), ncol=length(plate.ids))
    first.group <- conditions==1L
    de.genes <- seq_len(n)
    x[de.genes,first.group] <- x[de.genes,first.group] + 2
    
    y <- ComBat(x, plate.ids)
    cat("## Number of DE genes:", n, "\n")
    print(summary(rowMeans(y[de.genes,first.group])-rowMeans(y[de.genes,!first.group]))) # Differences squeezed inwards (but not totally lost, due to EB shrinkage).
    print(summary(rowMeans(x[de.genes,first.group])-rowMeans(x[de.genes,!first.group]))) # Differences retained.
    
    fit.y <- lmFit(y, design)
    fit.x <- lmFit(x, design)
    cat("## Differences between variances\n")
    print(summary(abs(fit.x$sigma-fit.y$sigma)/fit.x$sigma)) # Similar variances, so smaller differences = loss of power.
    cat("\n")
}

# In general, the more DE genes there are, the larger the estimate of the variance of the batch effect (as there's no robustness).
# This results in less shrinkage, which leads to a batch effect estimate closer to the DE effect, and reduction in the effect size upon correction.
