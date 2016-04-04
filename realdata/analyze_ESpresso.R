#######################################################
# Loading the ESpresso counts and experimental design.

curdir <- getwd()
setwd("../reference")
source("ESpresso.R")
setwd(curdir)

dir.create("ESpresso", showWarning=FALSE)
of.interest <- rowMeans(all.counts) >= 1
all.counts <- all.counts[of.interest,]
methods.to.use <- c("QLedgeR", "DESeq2", "voom")
coefs <- c(3, 4)

for (type in c("raw", "sum")) { 
    my.env <- new.env()
    my.env$sample.formula <- ~Batch + Serum 
    my.env$drop.coefficient <- coefs[1]
  
    if (type=="raw") { 
        my.env$counts <- all.counts
        my.env$sample.data <- targets
        my.env$normtype <- "deconvolution"  
        my.env$clusters <- by.group 
    } else {
        my.env$counts <- sumTechReps(all.counts, targets$Plate)
        my.env$sample.data <- targets[match(colnames(my.env$counts), targets$Plate),]
        my.env$normtype <- "deseq"  
    }
    
    sys.source("../reference/de_analysis.R", envir=my.env)

    # Analyzing with edgeR. 
    for (i in seq_along(coefs)) { 
        qres <- glmQLFTest(my.env$qfit, coef=coefs[i]) 
        fhandle <- gzfile(paste0("ESpresso/edgeR_", i, "_", type, ".tsv.gz"), open="wb")
        write.table(topTags(qres, n=Inf, sort.by="none")$table, file=fhandle,
                    col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
    }
    save(y, design, file=paste0("ESpresso/objects_", type, ".Rda")) # The DGEGLM objects are huge, so we'll go without.     

    # Analyzing with voom.
    for (i in seq_along(coefs)) { 
        vres <- topTable(my.env$vfit, n=Inf, sort.by="none", coef=coefs[i])
        fhandle <- gzfile(paste0("ESpresso/voom_", i, "_", type, ".tsv.gz"), open="wb")
        write.table(vres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
    }

    # Analyzing with DESeq2.
    for (i in seq_along(coefs)) {
        convec <- integer(ncol(my.env$design))
        convec[coefs[i]] <- 1
        dres <- results(my.env$dds, convec)

        fhandle <- gzfile(paste0("ESpresso/DESeq2_", i, "_", type, ".tsv.gz"), open="wb")
        write.table(dres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
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

for (con in c("1", "2")) {
    v.sum  <- read.table(paste0("ESpresso/voom_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    v.cell <- read.table(paste0("ESpresso/voom_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
    all.v <- compare.all(v.cell$adj.P.Val, v.sum$adj.P.Val)
    top.v <- compare.top(v.cell$P.Value, v.sum$P.Value)

    d.sum  <- read.table(paste0("ESpresso/DESeq2_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    d.cell <- read.table(paste0("ESpresso/DESeq2_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
    all.d <- compare.all(d.cell$padj, d.sum$padj)
    top.d <- compare.top(d.cell$pvalue, d.sum$pvalue)

    e.sum  <- read.table(paste0("ESpresso/edgeR_", con, "_sum.tsv.gz"), header=TRUE, row.names=1)
    e.cell <- read.table(paste0("ESpresso/edgeR_", con, "_raw.tsv.gz"), header=TRUE, row.names=1)
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
