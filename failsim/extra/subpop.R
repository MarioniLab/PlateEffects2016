# This demonstrates that the presence of subpopulations will not introduce
# dependencies between plates, as long as the subpopulation structure is the
# same.  Okay, we're using correlations as a measure of dependence, but it's
# probably good enough if we're not dealing with pathological scenarios.
# Check out 'plot(rowSums(y1), rowSums(y2), log="xy")' for independence.

set.seed(10001)

# Simplest case - a common plate effect:

subpop <- c(1,2,3)
popsize <- c(20, 30, 50)
ngenes <- 10000
each.cell <- rep(subpop, popsize)

for (it in seq_len(10)) { 
    a <- matrix(each.cell*10, byrow=TRUE, nrow=ngenes, ncol=sum(popsize))
    a1 <- a*exp(rnorm(ngenes))
    y1 <- matrix(rnbinom(length(a), mu=a1, size=20), nrow=ngenes)
    a2 <- a*exp(rnorm(ngenes))
    y2 <- matrix(rnbinom(length(a), mu=a2, size=20), nrow=ngenes)
    print(cor(rowSums(y1), rowSums(y2)))
}

# More complex case - subpopulation-specific plate effects:

subpop <- c(1,2,3)
popsize <- c(20, 30, 50)
ngenes <- 10000
each.cell <- rep(subpop, popsize)

for (it in seq_len(10)) { 
    a <- matrix(each.cell * 10, byrow=TRUE, nrow=ngenes, ncol=sum(popsize))
    plate.effect <- matrix(exp(rnorm(ngenes*length(popsize))), nrow=ngenes)
    a1 <- a*plate.effect[,each.cell]
    y1 <- matrix(rnbinom(length(a), mu=a1, size=20), nrow=ngenes)
    plate.effect <- matrix(exp(rnorm(ngenes*length(popsize))), nrow=ngenes)
    a2 <- a*plate.effect[,each.cell]
    y2 <- matrix(rnbinom(length(a), mu=a2, size=20), nrow=ngenes)
    print(cor(rowSums(y1), rowSums(y2)))
}

# Even more complicated - subpopulation-specific effects that are dependent within each plate:

subpop <- c(1,2)
popsize <- c(50, 50)
ngenes <- 10000
each.cell <- rep(subpop, popsize)

for (it in seq_len(10)) { 
    a <- matrix(each.cell * 10, byrow=TRUE, nrow=ngenes, ncol=sum(popsize))
    effect <- exp(rnorm(ngenes))
    plate.effect <- cbind(effect, 1/effect)
    a1 <- a*plate.effect[,each.cell]
    y1 <- matrix(rnbinom(length(a), mu=a1, size=20), nrow=ngenes)
    effect <- exp(rnorm(ngenes))
    plate.effect <- cbind(effect, 1/effect)
    a2 <- a*plate.effect[,each.cell]
    y2 <- matrix(rnbinom(length(a), mu=a2, size=20), nrow=ngenes)
    print(cor(rowSums(y1), rowSums(y2)))
}

# Printing the session information.

sessionInfo()
