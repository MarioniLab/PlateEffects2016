##################################################
# Setting up the simulation parameters.

set.seed(10000)

conditions <- rep(c("A", "B"), each=2)
nplates <- length(conditions)
ngenes <- 10000

all.means <- exp(runif(ngenes, 3, 8))
disp <- 0.5 + 100/all.means
plate.var <- 0.5

plate.grouping <- factor(conditions)
alphas <- c(0.01, 0.002, 0.05)

fc <- 3
first.set <- 1:1000
second.set <- first.set + 1000
combined.set <- c(first.set, second.set)

##################################################
# Odds and ends.

suppressPackageStartupMessages(require(edgeR))

compute.roc <- function(x) {
    x[is.na(x)] <- 1
    fn <- x[combined.set]
    tn <- x[-combined.set]
    return(findInterval(sort(fn), sort(tn))/length(tn))    
}

##################################################
# Running across all options.

for (pv in c(0, 0.5)) { 
    total.collected <- list()               
    for (scenario in 1:4) { 
        for (flip in c(TRUE, FALSE)) {
            edgeR.collected <- list()

            if (scenario==1L) {
                lib.sizes <- lapply(c(100, 100, 10, 10), FUN=rep, x=1)
            } else if (scenario==2L) {
                lib.sizes <- rep(list(rep(1, 50)), nplates)
                lib.sizes[[1]][1:20] <- lib.sizes[[2]][1:20] <- 5
            } else if (scenario==3L) {
                lib.sizes <- rep(list(rep(1, 50)), nplates)
                lib.sizes[[1]][1:10] <- lib.sizes[[2]][1:10] <- 10
            } else if (scenario==4L) { 
                lib.sizes <- rep(list(rep(1, 50)), nplates)
                lib.sizes[[1]][1:5] <- lib.sizes[[2]][1:5] <- 20
            }

            if (flip) { 
                tmp <- lib.sizes[[2]]
                lib.sizes[[2]] <- lib.sizes[[4]]
                lib.sizes[[4]] <- tmp
            }
 
            for (it in 1:10) {
                mean.per.plate <- all.means * exp(matrix(rnorm(nplates*ngenes, -pv/2, sqrt(pv)), ncol=nplates))

                # Adding balanced 2-fold DE to the first 2000 genes. Don't worry about DE being concentrated at high-abundances,
                # the fact that we have a spread of abundances for non-DE genes should mitigate any effects.
                mean.per.plate[first.set,conditions=="A"] <- mean.per.plate[first.set,conditions=="A"]*fc
                mean.per.plate[second.set,conditions=="B"] <- mean.per.plate[second.set,conditions=="B"]*fc

                collected <- list()
                for (x in seq_len(nplates)) {
                    mean.per.cell <- outer(mean.per.plate[,x], lib.sizes[[x]], "*")
                    collected[[x]] <- matrix(rnbinom(length(mean.per.cell), mu=mean.per.cell, size=1/disp), nrow=ngenes)
                }

                for (mode in c("raw", "sum")) {
                    if (mode=="raw") {
                        counts <- do.call(cbind, collected)
                        cell.grouping <- factor(rep(conditions, lapply(collected, ncol)))
                        design <- model.matrix(~cell.grouping)
                    } else {
                        counts <- do.call(cbind, lapply(collected, rowSums))
                        design <- model.matrix(~plate.grouping)
                    }
                    
                    # edgeR QL:
                    y <- DGEList(counts)
                    y <- calcNormFactors(y)
                    y <- estimateDisp(y, design)
                    fit <- glmQLFit(y, design, robust=TRUE)
                    res <- glmQLFTest(fit)
                    
                    edgeR.collected[[mode]][[it]] <- compute.roc(res$table$PValue)
                }
            }
            total.collected[[as.character(scenario * (-1)^flip)]] <- edgeR.collected
        }
    }
            
    saveRDS(file=ifelse(pv<1e-8, "without.rds", "with.rds"), total.collected)
}

##################################################
# End.

sessionInfo()
