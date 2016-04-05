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

#######################################################
# Running through the directories.

for (curdir in c("ESpresso")) {
    log.all <- file.path(curdir, "results_all.txt")
    log.top <- file.path(curdir, "results_top.txt")
    if (file.exists(log.all)) { unlink(log.all) }
    if (file.exists(log.top)) { unlink(log.top) }

    for (con in c("1", "2")) {
        v.sum  <- read.table(file.path(curdir, paste0("voom_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        v.cell <- read.table(file.path(curdir, paste0("voom_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        all.v <- compare.all(v.cell$adj.P.Val, v.sum$adj.P.Val)
        top.v <- compare.top(v.cell$P.Value, v.sum$P.Value)

        d.sum  <- read.table(file.path(curdir, paste0("DESeq2_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        d.cell <- read.table(file.path(curdir, paste0("DESeq2_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        all.d <- compare.all(d.cell$padj, d.sum$padj)
        top.d <- compare.top(d.cell$pvalue, d.sum$pvalue)

        e.sum  <- read.table(file.path(curdir, paste0("edgeR_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        e.cell <- read.table(file.path(curdir, paste0("edgeR_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        all.e <- compare.all(e.cell$FDR, e.sum$FDR)
        top.e <- compare.top(e.cell$PValue, e.sum$PValue)
        
        # Easier to put the separators in, as we need "&" at the front and "\\" at the end.
        write(sprintf("& %s & %i & %i & %i \\\\", c("Single-cell", "Summed", "Both"), unlist(all.d), unlist(all.v), unlist(all.e)), file=log.all, append=TRUE)
        write(file=log.all, "\\hline", append=TRUE)
        write(sprintf("& %s & %.2f & %.2f & %.2f \\\\", names(top.v), unlist(top.d), unlist(top.v), unlist(top.e)), file=log.top, append=TRUE) 
        write(file=log.top, "\\hline", append=TRUE)
    }
}

#######################################################
# Comparing to bulk, for the datasets that have them.

for (curdir in c("ESpresso")) {
    log.file <- file.path(curdir, "versus_bulk.txt")
    if (file.exists(log.file)) { unlink(log.file) }

    for (con in c("1", "2")) {
        v.sum  <- read.table(file.path(curdir, paste0("voom_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        v.bulk <- read.table(file.path(curdir, paste0("voom_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        top.v <- compare.top(v.bulk$P.Value, v.sum$P.Value)

        d.sum  <- read.table(file.path(curdir, paste0("DESeq2_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        d.bulk <- read.table(file.path(curdir, paste0("DESeq2_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        top.d <- compare.top(d.bulk$pvalue, d.sum$pvalue)

        e.sum  <- read.table(file.path(curdir, paste0("edgeR_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        e.bulk <- read.table(file.path(curdir, paste0("edgeR_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        top.e <- compare.top(e.bulk$PValue, e.sum$PValue)
        
        # Easier to put the separators in, as we need "&" at the front and "\\" at the end.
        write(sprintf("& %s & %.2f & %.2f & %.2f \\\\", names(top.v), unlist(top.d), unlist(top.v), unlist(top.e)), file=log.top, append=TRUE) 
        write(file=log.top, "\\hline", append=TRUE)
    } 
}

#######################################################

sessionInfo()

#######################################################
# End.
