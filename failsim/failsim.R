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

cell.grouping <- factor(rep(conditions, each=ncells))
design.all <- model.matrix(~cell.grouping)
plate.grouping <- factor(conditions)
design.sum <- model.matrix(~plate.grouping)
alphas <- c(0.01, 0.002, 0.05)

plate.of.origin <- rep(seq_len(nplates), each=ncells)

##################################################
# Generating the result file.

log.raw <- "temp_raw.txt"
if (file.exists(log.raw)) { stop("existing file for raw results") }
log.sum <- "temp_sum.txt"
if (file.exists(log.sum)) { stop("existing file for summed results") }

save.fun <- function(label, pv, x, log.file) {
    discard <- is.na(x)
    totes <- sum(!discard)
    for (alpha in alphas) {
        # Protect against NAs (ID'd as outliers by DESeq2).
        write.table(data.frame(label, pv, alpha, sum(x<=alpha & !discard)/totes), file=log.file, 
                    sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

##################################################
# Running across all options.

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(DESeq2))
suppressPackageStartupMessages(require(monocle))
suppressWarnings(suppressPackageStartupMessages(require(MAST)))

for (pv in c(0, plate.var)) {
    for (it in 1:10) {
        mean.per.plate <- all.means * exp(matrix(rnorm(nplates*ngenes, -pv/2, sqrt(pv)), ncol=nplates))
        all.counts <- matrix(0L, nrow=ngenes, ncol=ncells*nplates)
        for (x in seq_len(nplates)) {
            all.counts[,seq_len(ncells)+ncells*(x-1L)] <- rnbinom(ngenes*ncells, mu=mean.per.plate[,x], size=1/disp)
        }

        for (mode in c("raw", "sum")) {
            if (mode=="raw") {
                counts <- all.counts
                design <- design.all
                grp <- cell.grouping
                log.file <- log.raw
            } else {
                counts <- sumTechReps(all.counts, plate.of.origin)
                design <- design.sum
                grp <- plate.grouping
                log.file <- log.sum
            }

            # edgeR QL:
            y <- DGEList(counts)
            y <- calcNormFactors(y)
            y <- estimateDisp(y, design, prior.df=0)
            fit <- glmQLFit(y, design, robust=TRUE)
            res <- glmQLFTest(fit)
            save.fun("edgeR (QL)", pv, res$table$PValue, log.file)

            if (mode=="raw") {
                # edgeR LRT:
                fit <- glmFit(y, design, dispersion=y$tagwise.dispersion)
                res <- glmLRT(fit)
                save.fun("edgeR (LRT)", pv, res$table$PValue, log.file)
            }

            # DESeq2
            suppressMessages(dds <- DESeqDataSetFromMatrix(counts, colData=DataFrame(grp=grp), design = ~grp))
            suppressMessages(dds <- DESeq(dds))
            res <- results(dds, c("grp", "A", "B"))
            save.fun("DESeq2", pv, res$pvalue, log.file)

            # voom without correlations
            v.all <- voom(y, design)
            fit <- lmFit(v.all, design)
            fit <- eBayes(fit, robust=TRUE)
            res <- topTable(fit, n=Inf, sort.by="none", coef=2)
            save.fun("voom", pv, res$P.Value, log.file)

            if (mode=="raw") {
                # voom with correlations
                dc <- duplicateCorrelation(v.all, design, block=plate.of.origin)
                fit <- lmFit(v.all, design, block=plate.of.origin, correlation=dc$consensus)
                fit <- eBayes(fit, robust=TRUE)
                res <- topTable(fit, n=Inf, sort.by="none", coef=2)
                save.fun("voom (cor)", pv, res$P.Value, log.file)
            }

            if (mode=="raw") {
                # MAST
                oldseed <- .Random.seed
                suppressMessages({
                    cpms <- cpm(counts+1, prior.count=0, log=TRUE, lib.size=colSums(counts))
                    sca <- FromMatrix('SingleCellAssay', t(cpms), data.frame(wellKey=seq_along(grp)), data.frame(primerid=seq_len(ngenes)))
                    cData(sca)$cngeneson <- colMeans(counts>0)
                    cData(sca)$condition <- grp
                    fit <- zlm.SingleCellAssay(~ condition + cngeneson, sca, method="bayesglm", ebayes=TRUE, ebayesControl=list(method="MLE", model="H1"))
                    lrt <- lrTest(fit, "condition")
                    save.fun("MAST", pv, lrt[, "hurdle", "Pr(>Chisq)"], log.file)
                })
                .Random.seed <- oldseed # Avoid getting different results because of MAST's randomization methods.
            }

            # SAMstrt doesn't report p-values, just "median FDR" estimates.
            # It also relies on spike-ins, that we don't have here.
        
            if (it==1L && mode=="raw") {
                # Only running monocle in one iteration, as it takes too long.
                cpms <- cpm(counts, prior.count=0)
                pdat <- AnnotatedDataFrame(data=data.frame(grp=grp))
                sampleNames(pdat) <- colnames(cpms)
                HSMM <- new("CellDataSet", exprs=cpms, phenoData=pdat, expressionFamily=negbinomial())
                out <- differentialGeneTest(HSMM, fullModelFormulaStr="expression~grp", cores=6) 
                save.fun("monocle", pv, out$pval, log.file)
            }         
        }
    }
}

##################################################
# End.

sessionInfo()
