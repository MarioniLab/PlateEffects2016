##################################################
# Setting up the simulation parameters.

set.seed(10000)

conditions <- rep(c("A", "B"), each=3)
nplates <- length(conditions)
ngenes <- 10000
ncells <- 50

all.means <- exp(runif(ngenes, 3, 8))
disp <- 0.5 + 100/all.means
plate.var <- 0.5
fc <- 3

cell.grouping <- factor(rep(conditions, each=ncells))
design.all <- model.matrix(~cell.grouping)
plate.grouping <- factor(conditions)
design.sum <- model.matrix(~plate.grouping)

alphas <- c(1e-5, 1e-4, 1e-3, 0.01, 0.1)
alphas <- sort(c(alphas, alphas*3))
first.set <- 1:1000
second.set <- first.set + 1000
combined.set <- c(first.set, second.set)
plate.of.origin <- rep(seq_len(nplates), each=ncells)

##################################################

compute.roc <- function(x) {
    x[is.na(x)] <- 1
    fn <- x[combined.set]
    tn <- x[-combined.set]
    return(findInterval(sort(fn), sort(tn))/length(tn))    
}

##################################################
# Running across all options.

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(DESeq2))

for (pv in c(0, plate.var)) {
    edgeR.collected <- DESeq2.collected <- voom.collected <- list(sum=list(), raw=list())

    for (it in 1:10) {
        mean.per.plate <- all.means * exp(matrix(rnorm(nplates*ngenes, -pv/2, sqrt(pv)), ncol=nplates))
        
        # Adding balanced 2-fold DE to the first 2000 genes.
        mean.per.plate[first.set,conditions=="A"] <- mean.per.plate[first.set,conditions=="A"]*fc
        mean.per.plate[second.set,conditions=="B"] <- mean.per.plate[second.set,conditions=="B"]*fc

        all.counts <- matrix(0L, nrow=ngenes, ncol=ncells*nplates)
        for (x in seq_len(nplates)) {
            all.counts[,seq_len(ncells)+ncells*(x-1L)] <- rnbinom(ngenes*ncells, mu=mean.per.plate[,x], size=1/disp)
        }

        for (mode in c("raw", "sum")) {
            if (mode=="raw") {
                counts <- all.counts
                design <- design.all
                grp <- cell.grouping
            } else {
                counts <- sumTechReps(all.counts, plate.of.origin)
                design <- design.sum
                grp <- plate.grouping
            }

            # edgeR QL:
            y <- DGEList(counts)
            y <- calcNormFactors(y)
            y <- estimateDisp(y, design, prior.df=0)
            fit <- glmQLFit(y, design, robust=TRUE)
            res <- glmQLFTest(fit)
            edgeR.collected[[mode]][[it]] <- compute.roc(res$table$PValue)

            # DESeq2
            suppressMessages(dds <- DESeqDataSetFromMatrix(counts, colData=DataFrame(grp=grp), design = ~grp))
            suppressMessages(dds <- DESeq(dds))
            res <- results(dds, c("grp", "A", "B"))
            DESeq2.collected[[mode]][[it]] <- compute.roc(res$pvalue)

            # voom without correlations
            v.all <- voom(y, design)
            fit <- lmFit(v.all, design)
            fit <- eBayes(fit, robust=TRUE)
            res <- topTable(fit, n=Inf, sort.by="none", coef=2)
            voom.collected[[mode]][[it]] <- compute.roc(res$P.Value)
        }
    }

    saveRDS(file=ifelse(pv<1e-8, "without.rds", "with.rds"), list(edgeR=edgeR.collected, DESeq2=DESeq2.collected, voom=voom.collected))
}

##################################################
# End.

sessionInfo()
