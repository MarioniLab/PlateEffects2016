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
design <- model.matrix(~plate.grouping)
alphas <- c(0.01, 0.002, 0.05)

##################################################
# Generating the result file.

log.file <- "temp.txt"
if (file.exists(log.file)) { stop("existing file for results") }

save.fun <- function(label, scenario, pv, mode, x, log.file) {
    discard <- is.na(x)
    totes <- sum(!discard)
    for (alpha in alphas) {
        write.table(data.frame(label, scenario, pv, mode, alpha, sum(x<=alpha & !discard)/totes), file=log.file, 
                    sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

suppressPackageStartupMessages(require(edgeR))

sumPoisFun <- function(counts, mean, disp) {
    mean <- rep(mean, length.out=length(counts))
    ri <- 1 + disp*mean
    vi <- mean*ri
    ro <- rep(1, length(counts))
    vo <- mean
    upper <- (counts >= mean)
    p1 <- p2 <- q1 <- q2 <- counts
    if (any(upper)) {
        p1[upper] <- pnorm(counts[upper], mean=mean[upper],
            sd=sqrt(vi[upper]), lower.tail=FALSE, log.p=TRUE)
        p2[upper] <- pgamma(counts[upper], shape=mean[upper]/ri[upper],
            scale=ri[upper], lower.tail=FALSE, log.p=TRUE)
        q1[upper] <- qnorm(p1[upper], mean=mean[upper], sd=sqrt(vo[upper]),
            lower.tail=FALSE, log.p=TRUE)
        q2[upper] <- qgamma(p2[upper], shape=mean[upper]/ro[upper],
            scale=ro[upper], lower.tail=FALSE, log.p=TRUE)
    }
    lower <- !upper
    if (any(lower)) {
        p1[lower] <- pnorm(counts[lower], mean=mean[lower],
            sd=sqrt(vi[lower]), lower.tail=TRUE, log.p=TRUE)
        p2[lower] <- pgamma(counts[lower], shape=mean[lower]/ri[lower],
            scale=ri[lower], lower.tail=TRUE, log.p=TRUE)
        q1[lower] <- qnorm(p1[lower], mean=mean[lower], sd=sqrt(vo[lower]),
            lower.tail=TRUE, log.p=TRUE)
        q2[lower] <- qgamma(p2[lower], shape=mean[lower]/ro[lower],
            scale=ro[lower], lower.tail=TRUE, log.p=TRUE)
    }
    out <- (q1 + q2)/2 # Approximated quantiles here
    out <- rowSums(out)
    out[out < 0] <- 0
    storage.mode(out) <- "integer"
    out
}

estimateParam <- function(count.list) {
    combo <- rep(seq_along(count.list), sapply(count.list, FUN=ncol))
    super.design <- model.matrix(~factor(combo))
    all.counts <- do.call(cbind, count.list)
    y <- DGEList(all.counts)
    y <- estimateDisp(y, super.design, prior.df=0)
    out <- glmFit(y, super.design, dispersion=y$tagwise.dispersion)

    # Splitting for convenience.
    nrows <- nrow(all.counts)
    fitted <- split(t(out$fitted.values), combo)
    fitted <- lapply(fitted, FUN=function(u) matrix(u, nrow=nrows, byrow=TRUE))
    return(list(fitted=fitted, dispersion=y$tagwise.dispersion))
}

##################################################
# Running across all options.

for (pv in c(0, 0.5)) { 
    for (scenario in 1:4) { 
        for (flip in c(TRUE, FALSE)) {

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
                collected <- list()
                for (x in seq_len(nplates)) {
                    mean.per.cell <- outer(mean.per.plate[,x], lib.sizes[[x]], "*")
                    collected[[x]] <- matrix(rnbinom(length(mean.per.cell), mu=mean.per.cell, size=1/disp), nrow=ngenes)
                }

                for (mode in c("quantile", "sum")) {
                    if (mode=="quantile") {
                        estimates <- estimateParam(collected)
                        counts <- matrix(0L, nrow=ngenes, ncol=nplates)  
                        for (i in seq_len(nplates)) {
                            counts[,i] <- suppressWarnings(sumPoisFun(collected[[i]], estimates$fitted[[i]], estimates$dispersion))
                        }
                        # Removing NA counts. 
                        counts <- counts[rowSums(is.na(counts))==0L,]
                    } else {
                        counts <- do.call(cbind, lapply(collected, rowSums))
                    }

                    # edgeR QL:
                    y <- DGEList(counts)
                    y <- calcNormFactors(y)
                    y <- estimateDisp(y, design, prior.df=0)
                    fit <- glmQLFit(y, design, robust=TRUE)
                    res <- glmQLFTest(fit)
                    save.fun("edgeR (QL)", (-1L)^flip*scenario, pv, mode, res$table$PValue, log.file)
                }
            }
        }
    }
}

##################################################
# End.

sessionInfo()
