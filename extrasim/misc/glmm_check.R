########################################
# This simulation checks how GLMMs perform in a simple case, to try to explain
# the loss of type I error control in our complex simulations.

library(lme4)
ngenes <- 10000
ngroups <- 2
set.seed(10000)

for (scenario in 1:6) {
    if (scenario<=2L) {
        nbatch <- 2
        if (scenario==1L) { 
            nreps <- 2
        } else { 
            nreps <- 10
        }
        mean.vals <- matrix(exp(5+1/2), nrow=ngenes, ncol=nbatch*2)
    } else {
        if (scenario==3L) {
            nbatch <- 2
            nreps <- 2
        } else if (scenario==4L) {
            nbatch <- 10
            nreps <- 2
        } else if (scenario==5L) {
            nbatch <- 2
            nreps <- 10
        } else {
            nbatch <- 10
            nreps <- 10
        }
        mean.vals <- matrix(exp(rnorm(ngenes*nbatch*2, mean=5)), nrow=ngenes)
    }

    Condition <- factor(rep(LETTERS[seq_len(ngroups)], each=nbatch*nreps))
    Batch <- factor(rep(seq_len(nbatch*ngroups), each=nreps))
    Data <- data.frame(Batch=Batch, Condition=Condition)
    output <- rep(NA, ngenes)
    var.val <- rep(NA, ngenes)

    for (i in seq_along(output)) {
        Data$Y <- rpois(length(Batch), lambda=mean.vals[i,Batch])
        try({ 
            le <- glmer(Y ~ Condition + (1|Batch), data=Data, family=poisson)
            le0 <- glmer(Y ~ (1|Batch), data=Data, family=poisson)
            output[i] <- anova(le, le0)["le","Pr(>Chisq)"] # Avoid using Wald Z-tests, using LRT instead.
            var.val[i] <- unlist(VarCorr(le))
        })
    }
    cat("Scenario is", scenario, '\n')
    print(mean(output <= 0.01, na.rm=TRUE))
    print(mean(var.val, na.rm=TRUE))
}

# There are several take-home points from this simulation:
# - Type I error control is maintained (and in fact, conservative)  when there is no batch effect.
#   This is probably because underestimation of the variance of the random effect is impossible when the true variance is zero!
#   It also indicates that the asymptotic approximations used to compute the p-value are not inherently liberal.
# - With a non-zero variance for the batch term, type I error control is lost.
#   This is associated with an underestimation of the variance of the random effect (it should be 1), which probably causes the liberalness.
#   Intiuitively, the idea is that the standard error of the coefficients would be underestimated.
# - As the number of batches increases, liberalness is mitigated and the variance approaches the correct value.
#   This reflects the increased amount of information available to estimate the variance.
# - The number of replicates per batch has little effect.

