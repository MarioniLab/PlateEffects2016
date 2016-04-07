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
    if (curdir == "ESpresso") {
        coefs <- c("3", "4")
    }

    for (con in coefs) {
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
        write(sprintf("%% %s %s", curdir, con), file=log.all, append=TRUE)
        write(sprintf("& %s & %i & %i & %i \\\\", c("Single-cell", "Summed", "Both"), unlist(all.d), unlist(all.v), unlist(all.e)), file=log.all, append=TRUE)
        write(file=log.all, "\\hline", append=TRUE)
        write(sprintf("& %s & %.2f & %.2f & %.2f \\\\", names(top.v), unlist(top.d), unlist(top.v), unlist(top.e)), file=log.top, append=TRUE) 
        write(file=log.top, "\\hline", append=TRUE)
    }
}

#######################################################
# Comparing to bulk, for the datasets that have them.

getstat <- function(bulk, other, threshold=0.05) {
    both.okay <- !is.na(bulk) & !is.na(other)
    bulksig <- bulk <= threshold 
    othersig <- other <= threshold 
    sum(bulksig & othersig & both.okay)/sum((othersig | bulksig) & both.okay)
}

log.file <- file.path(curdir, "results_bulk.txt")
if (file.exists(log.file)) { unlink(log.file) }

for (curdir in c("ESpresso")) {
     if (curdir == "ESpresso") {
        coefs <- c("3", "4")
    }

    for (con in coefs) {
        v.sum  <- read.table(file.path(curdir, paste0("voom_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        v.cell <- read.table(file.path(curdir, paste0("voom_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        v.bulk <- read.table(file.path(curdir, paste0("voom_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        stat.vc <- getstat(v.bulk$adj.P.Val, v.cell$adj.P.Val)
        stat.vs <- getstat(v.bulk$adj.P.Val, v.sum$adj.P.Val)

        d.sum  <- read.table(file.path(curdir, paste0("DESeq2_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        d.cell <- read.table(file.path(curdir, paste0("DESeq2_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        d.bulk <- read.table(file.path(curdir, paste0("DESeq2_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        stat.dc <- getstat(d.bulk$padj, d.cell$padj)
        stat.ds <- getstat(d.bulk$padj, d.sum$padj)

        e.sum  <- read.table(file.path(curdir, paste0("edgeR_", con, "_sum.tsv.gz")), header=TRUE, row.names=1)
        e.cell <- read.table(file.path(curdir, paste0("edgeR_", con, "_raw.tsv.gz")), header=TRUE, row.names=1)
        e.bulk <- read.table(file.path(curdir, paste0("edgeR_", con, "_bulk.tsv.gz")), header=TRUE, row.names=1)
        stat.ec <- getstat(e.bulk$FDR, e.cell$FDR)
        stat.es <- getstat(e.bulk$FDR, e.sum$FDR)
        
        write(sprintf("%% %s %s", curdir, con), file=log.file, append=TRUE)
        write(sprintf("& %s & %.2f & %.2f & %.2f \\\\", c("Single-cell", "Summed"), c(stat.dc, stat.ds), c(stat.vc, stat.vs), c(stat.ec, stat.es)), file=log.file, append=TRUE)
        write(file=log.file, "\\hline", append=TRUE)
    } 
}

#######################################################
# Generating a bar plot for the scrambled results.

all.scrambles <- list()
all.colors <- list()
leg.contrasts <- list()
leg.colors <- list()

i <- 1
for (curdir in c("ESpresso")) {
    current <- read.table(file.path(curdir, "scrambled.txt"))
    if (curdir=="ESpresso") { 
        contrasted <- c("3"="grey80", "4"="grey40")
        leg.contrasts[[i]] <- c('2i vs. serum', 'a2i vs. serum')
    }
    mode <- paste0(current[,1], current[,3])
    colors <- split(contrasted[as.character(current[,2])], mode)
    fprs <- split(current[,4], mode)

    all.scrambles[[i]] <- do.call(cbind, fprs)
    all.colors[[i]] <- do.call(cbind, colors)
    leg.colors[[i]] <- contrasted
    i <- i + 1
}

scrambled <- do.call(rbind, all.scrambles)
colors <- do.call(rbind, all.colors)

# Renaming.
renamed <- colnames(scrambled)
renamed <- sub("sum$", "\n(sum)", renamed)
renamed <- sub("raw$", "\n(single)", renamed)
renamed <- sub("QLedgeR", "QL edgeR", renamed)

# Plotting.
setEPS()
width <- 9
postscript("all_scrambled.eps", width=width, height=7)
layout(rbind(c(1,2)), width=c(width-2, 2))
par(mar=c(6.2, 5.1, 2.1, 0.1))

barplot(scrambled, names.arg=renamed, beside=TRUE, col=colors, ylab="Proportion of genes rejected", 
        las=2, cex.lab=1.4, cex.axis=1.2, cex.names=1.2)
abline(h=0.01, col="red", lwd=2, lty=2)

par(mar=c(0.2, 0.1, 5.1, 0.1))
plot(0, 0, type="n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
legend("topleft", legend=unlist(leg.contrasts), fill=unlist(leg.colors), cex=1.2)
dev.off()

#######################################################

sessionInfo()

#######################################################
# End.
