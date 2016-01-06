#suppressPackageStartupMessages(require(simpaler))
#suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.ensGene))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(DESeq2))

for (data in c("ola", "scialdone")) {
    if (data=="ola") {
        host.dir <- "."
        all.counts <- read.table("ESpresso/counttable_es.csv", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 704)))
        serum <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(all.counts))
        batch <- sub("ola_mES_[^_]+_([^_]+)_.*", "\\1", colnames(all.counts))
        targets <- data.frame(Serum=serum, Batch=batch)

        # Only using data from two batches.
        keep <- targets$Batch %in% c("2", "3")
        all.counts <- all.counts[,keep]
        targets <- targets[keep,]
        targets$Plate <- as.integer(factor(paste0(targets$Serum, targets$Batch)))
        targets[] <- lapply(targets, factor)
        targets$Serum <- factor(targets$Serum, c("lif", "2i", "a2i"))
    } else {
        whee <- new.env()
        load("Scialdone_Gastrulation/aaron_data_gastrulation.RData", env=whee)
        all.counts <- whee$data 
        targets <- data.frame(whee$labels)
        targets[] <- lapply(targets, factor)
    }

    # Quality control on genes (QC on cells already done).
    is.mouse <- grepl("^ENSMUSG", rownames(all.counts))

#    chr.loc <- findChr(rownames(all.counts[is.mouse,]), TxDb.Mmusculus.UCSC.mm10.ensGene)
#    lib.sizes <- colSums(all.counts[is.mouse,])
#    okay.libs <- lib.sizes > 1e5 & colSums(all.counts[is.mouse,][!is.na(chr.loc) & chr.loc=="chrM",])/lib.sizes < 0.1

#    all.counts <- all.counts[,okay.libs]
#    targets <- targets[okay.libs,]

    of.interest <- is.mouse & rowMeans(all.counts) >= 1
    all.counts <- all.counts[of.interest,]

    # Formulating the design matrix.
    if (data=="ola") {
        model.formula <- ~Batch + Serum
        comparisons <- list(c("Serum", "2i", "lif"), c("Serum", "a2i", "lif"))
        coefs <- c(3, 4)
    } else {
        model.formula <- ~Stage
        comparisons <- list(c("Stage", "NP", "PS"))
        coefs <- 2
    }
    refdesign <- model.matrix(model.formula, data=targets)
   
    for (type in c("raw", "sum")) { 
        if (type=="raw") { 
            counts <- all.counts
            design <- refdesign
            deseq.target <- targets
        } else {
            counts <- sumTechReps(all.counts, targets$Plate)
            first.in.each <- !duplicated(targets$Plate)
            design <- refdesign[first.in.each,]
            deseq.target <- targets[first.in.each,]
        }

        # Analyzing with edgeR. 
        y <- DGEList(counts)
        y <- calcNormFactors(y)
        y <- estimateDisp(y, design, prior.df=0)
        fit <- glmQLFit(y, design, robust=TRUE)

        for (i in seq_along(coefs)) { 
            res <- glmQLFTest(fit, coef=coefs[i]) 
            fhandle <- gzfile(paste0("results/edgeR_", data, i, "_", type, ".tsv.gz"), open="wb")
            write.table(topTags(res, n=Inf, sort.by="none")$table, file=fhandle,
                        col.names=NA, sep="\t", quote=FALSE)
            close(fhandle)

        }
        save(y, design, file=paste0("results/objects_", data, "_", type, ".Rda")) # The DGEGLM objects are huge, so we'll go without.     

        # Analyzing with voom.
        v.all <- voom(y, design)
        fit <- lmFit(v.all, design)
        fit <- eBayes(fit, robust=TRUE)

        for (i in seq_along(coefs)) { 
            res <- topTable(fit, n=Inf, sort.by="none", coef=coefs[i])
            fhandle <- gzfile(paste0("results/voom_", data, i, "_", type, ".tsv.gz"), open="wb")
            write.table(res, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
            close(fhandle)
        }

        # Analyzing with DESeq2.
        suppressMessages(dds <- DESeqDataSetFromMatrix(counts, colData=deseq.target, design=model.formula))
        suppressMessages(dds <- DESeq(dds))
        for (i in seq_along(comparisons)) {
            res <- results(dds, comparisons[[i]])
            fhandle <- gzfile(paste0("results/DESeq2_", data, i, "_", type, ".tsv.gz"), open="wb")
            write.table(res, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
            close(fhandle)
        }
    }
}

#######################################################
# Comparing the various types of analyses.

compare.all <- function(qval1, qval2, threshold=0.05) {
    sig1 <- qval1 <= threshold & !is.na(qval1)
    sig2 <- qval2 <= threshold & !is.na(qval2)
    return(list(total1=sum(sig1), total2=sum(sig2), shared=sum(sig1 & sig2)))
}

compare.top <- function(pval1, pval2, top=c(20, 200, 2000)) {
    output <- list()
    o1 <- rank(pval1, ties.method="first")
    o2 <- rank(pval2, ties.method="first") 
    for (x in top) {
        output[[as.character(x)]] <- sum(o1 <= x & o2 <= x)/x
    }
    return(output)
}

log.all <- "results_all.txt"
log.top <- "results_top.txt"
if (file.exists(log.all)) { unlink(log.all) }
if (file.exists(log.top)) { unlink(log.top) }

for (con in c("ola1", "ola2", "scialdone1")) {
    v.sum  <- read.table(paste0("results/voom_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    v.cell <- read.table(paste0("results/voom_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
    all.v <- compare.all(v.cell$adj.P.Val, v.sum$adj.P.Val)
    top.v <- compare.top(v.cell$P.Value, v.sum$P.Value)

    d.sum  <- read.table(paste0("results/DESeq2_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    d.cell <- read.table(paste0("results/DESeq2_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
    all.d <- compare.all(d.cell$padj, d.sum$padj)
    top.d <- compare.top(d.cell$pvalue, d.sum$pvalue)

    e.sum  <- read.table(paste0("results/edgeR_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    e.cell <- read.table(paste0("results/edgeR_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
    all.e <- compare.all(e.cell$FDR, e.sum$FDR)
    top.e <- compare.top(e.cell$PValue, e.sum$PValue)
    
    # Easier to put the separators in, as we need "&" at the front and "\\" at the end.
    write(sprintf("& %s & %i & %i & %i \\\\", c("Single-cell", "Summed", "Both"), unlist(all.d), unlist(all.v), unlist(all.e)), file=log.all, append=TRUE)
    write(file=log.all, "\\hline", append=TRUE)
    write(sprintf("& %s & %.2f & %.2f & %.2f \\\\", names(top.v), unlist(top.d), unlist(top.v), unlist(top.e)), file=log.top, append=TRUE) 
    write(file=log.top, "\\hline", append=TRUE)
}

#######################################################

sessionInfo()

#######################################################
# End.
