####################################################################################
# Need:

counts
sample.formula
sample.data
drop.coefficient 
methods.to.use 

design <- model.matrix(sample.formula, sample.data)
obtained <- list()

####################################################################################
# Computing size factors from the deconvolution method.

if (!exists("normtype")) { normtype <- "libsize" }
normtype <- match.arg(normtype, c("libsize", "deconvolution", "deseq")) 
if (normtype=="libsize") {
    sf <- colSums(counts)
    sf <- sf/mean(sf) 
} else if (normtype=="deseq") { 
    suppressPackageStartupMessages(library(DESeq2))
    sf <- estimateSizeFactorsForMatrix(counts)
} else if (normtype=="deconvolution") {
    suppressPackageStartupMessages(library(scran))
    if (!exists('clusters')) { clusters <- NULL }
    sf <- computeSumFactors(counts, cluster=clusters)
}

eff.lib <- sf*mean(colSums(counts)) 

####################################################################################
# edgeR and variants:

if (any(c("edgeR", "voom", "voomcor", "QLedgeR") %in% methods.to.use)) {
    suppressPackageStartupMessages(library(edgeR))
    y <- DGEList(counts, lib.size=eff.lib) 

    if (any(c("edgeR", "QLedgeR") %in% methods.to.use)) {
        y <- estimateDisp(y, design, prior.df=0)
        
        if ("QLedgeR" %in% methods.to.use) {
            qfit <- glmQLFit(y, design, robust=TRUE)
            qres <- glmQLFTest(qfit, coef=drop.coefficient)
            obtained$QLedgeR <- qres$table$PValue
        }
        if ("edgeR" %in% methods.to.use) {
            fit <- glmFit(y, design, dispersion=y$tagwise.dispersion)
            res <- glmLRT(fit, coef=drop.coefficient)
            obtained$edgeR <- res$table$PValue
        }
    }
}

####################################################################################
# DESeq2:

if ("DESeq2" %in% methods.to.use)  {
    suppressPackageStartupMessages(library(DESeq2))
    suppressMessages(dds <- DESeqDataSetFromMatrix(counts, colData=sample.data, design=sample.formula))
    sizeFactors(dds) <- sf
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds, modelMatrixType="standard")
    
    convec <- integer(ncol(design))
    convec[drop.coefficient] <- 1
    dres <- results(dds, convec)
    obtained$DESeq2 <- dres$pvalue
}

####################################################################################
# voom with or without correlations

if (any(c("voom", "voomcor") %in% methods.to.use)) {
    v.all <- voom(y, design)

    if ("voom" %in% methods.to.use) {
        vfit <- lmFit(v.all, design)
        vfit <- eBayes(vfit, robust=TRUE)
        vres <- topTable(vfit, n=Inf, sort.by="none", coef=drop.coefficient)
        obtained$voom <- vres$P.Value
    } 
    if ("voomcor" %in% methods.to.use) {
        dc <- duplicateCorrelation(v.all, design, block=plate.of.origin)
        vcfit <- lmFit(v.all, design, block=plate.of.origin, correlation=dc$consensus)
        vcfit <- eBayes(vcfit, robust=TRUE)
        vcres <- topTable(vcfit, n=Inf, sort.by="none", coef=drop.coefficient)
        obtained$voomcor <- vcres$P.Value
    }

#    # Voom in conjunction with ComBat (basically detects nothing).
#    # I suspect empirical Bayes doesn't change much here, because the batch effect is precisely estimated with many cells per batch.
#    library(sva)
#    combat.out <- ComBat(dat=my.env$v.all$E, batch=plate.of.origin, mod=model.matrix(~1, my.env$sample.data), par.prior=TRUE)
#    v.all2 <- v.all
#    v.all2$E <- combat.out
#    vfit2 <- lmFit(v.all2, design)
#    vfit2 <- eBayes(vfit2, robust=TRUE)
#    vres2 <- topTable(vfit2, n=Inf, sort.by="none", coef=drop.coefficient)
}

####################################################################################
# MAST:

if ("MAST" %in% methods.to.use) {
    suppressPackageStartupMessages(library(MAST))
    oldseed <- .Random.seed
    suppressMessages({
        lcpms <- log2(t(counts+1)/eff.lib*1e6)
        sca <- FromMatrix('SingleCellAssay', lcpms, data.frame(wellKey=seq_len(nrow(sample.data))), data.frame(primerid=seq_len(nrow(counts))))
        cData(sca) <- cbind(cData(sca), sample.data)
        cData(sca)$cngeneson <- colMeans(counts>0)

        new.formula <- paste(c(as.character(sample.formula), "+ cngeneson"), collapse=" ")
        new.formula <- as.formula(new.formula)
        mfit <- zlm.SingleCellAssay(new.formula, sca, method="bayesglm", ebayes=TRUE, ebayesControl=list(method="MLE", model="H1"))
        mlrt <- lrTest(mfit, Hypothesis(colnames(design)[drop.coefficient], colnames(coef(mfit, "D")))) # cngeneson is added at the end, so it shouldn't affect terms before it.
        obtained$MAST <- mlrt[, "hurdle", "Pr(>Chisq)"]
    })
    .Random.seed <- oldseed # Avoid getting different results because of MAST's randomization methods.
}

####################################################################################
# SAMstrt doesn't report p-values, just "median FDR" estimates.
# It also relies on spike-ins, that we don't have here.

####################################################################################
# monocle:

if ("monocle" %in% methods.to.use) {
    suppressPackageStartupMessages(library(monocle))
    cpms <- t(t(counts)/eff.lib * 1e6)

    # Whole lot of effort to automatically construct a contrast.
    pdat <- AnnotatedDataFrame(as.data.frame(design))
    sampleNames(pdat) <- colnames(cpms)
    varLabels(pdat) <- paste0("X", seq_len(ncol(design)))
    intercept <- colSums(design)==nrow(design)
    if (!any(intercept)) { stop("need an intercept here for monocle") }
    intercept <- which(intercept)

    HSMM <- newCellDataSet(cellData=cpms, phenoData=pdat)
    out <- differentialGeneTest(HSMM, cores=6, 
        fullModelFormulaStr=paste0("expression~", paste(varLabels(pdat)[-intercept], collapse="+")),
        reducedModelFormulaStr=paste0("expression~", ifelse(ncol(pdat)==2L, 1, 
                                      paste(varLabels(pdat)[-c(intercept, drop.coefficient)], collapse="+"))))
    obtained$monocle <- out$pval
}

####################################################################################
# Mixed modelling:

if ("glmer" %in% methods.to.use) {
    my.data <- data.frame(design)
    colnames(my.data) <- paste0("X", seq_len(ncol(design)))
    full.form <- as.formula(paste("Counts ~ 0 + ", paste(colnames(my.data), collapse="+"), "+ (1|Plate) + offset(log(sf))"))
    null.form <- as.formula(paste("Counts ~ 0 +", paste(colnames(my.data)[-drop.coefficient], collapse="+"), "+ (1|Plate) + offset(log(sf))"))
    my.data$Counts <- 0
    my.data$Plate <- factor(plate.of.origin)

    suppressPackageStartupMessages(library(lme4))
    pvalues <- rep(NA, nrow(counts))
    for (g in seq_len(ngenes)) { 
        my.data$Counts <- counts[g,]
        try({
            full.fit <- glmer.nb(full.form, data=my.data)
            null.fit <- glmer.nb(null.form, data=my.data)
            tested <- anova(full.fit, null.fit, test = "Chisq")
            # see ?glmer, http://glmm.wikidot.com/faq; best option other than bootstrapping (which would be bothersome for so many genes).
            # also see http://lme4.r-forge.r-project.org/slides/2010-09-23-Rahway/7GLMM.pdf.
            pvalues[g] <- tested$Pr[2] 
        }, silent=TRUE)
    }
    obtained$glmer <- pvalues
}

####################################################################################
