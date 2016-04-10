# This makes plots about the variability of the top set of genes for the Kolod data,
# using cells in the first batch for easy visualization (not Scialdone, as there aren't any DE genes).

curdir <- getwd()
setwd("../reference")
source("ESpresso.R")
coefs <- c(3, 4)
setwd(curdir)

for (x in coefs) {
    current <- targets$Serum %in% c(ifelse(x==3, '2i', 'a2i'), "lif") & targets$Batch == "2"
    data <- all.counts[,current]
    grouping <- targets$Serum[current]
    require(edgeR)
    cpms <- cpm(data, log=TRUE)

    method <- "edgeR"
    cur.field <- switch(method, edgeR="PValue", voom="P.Value", DESeq2="pvalue")
    sum.res <- read.table(sprintf("ESpresso/%s_%s_sum.tsv.gz", method, x), header=TRUE, row.names=1)
    sum.res <- sum.res[order(sum.res[[cur.field]]),]
    sum.top <- rownames(sum.res)[1:10]
    raw.res <- read.table(sprintf("ESpresso/%s_%s_raw.tsv.gz", method, x), header=TRUE, row.names=1)
    raw.res <- raw.res[order(raw.res[[cur.field]]),]
    raw.top <- rownames(raw.res)[1:10]
    
    unique.sum <- setdiff(sum.top, raw.top)
    unique.raw <- setdiff(raw.top, sum.top)
    sum.values <- cpms[rownames(cpms) %in% unique.sum,]
    raw.values <- cpms[rownames(cpms) %in% unique.raw,]
    
    require(org.Mm.eg.db)
    sum.names <- select(org.Mm.eg.db, keys=rownames(sum.values), columns="SYMBOL", keytype="ENSEMBL")$SYMBOL
    missing.sum <- is.na(sum.names)
    sum.names[missing.sum] <- rownames(sum.values)[missing.sum]
    raw.names <- select(org.Mm.eg.db, keys=rownames(raw.values), columns="SYMBOL", keytype="ENSEMBL")$SYMBOL
    missing.raw <- is.na(raw.names)
    raw.names[missing.raw] <- rownames(raw.values)[missing.raw]
    
    set.seed(1000)
    pdf(sprintf("ESpresso/top_%s_%s_sum.pdf", method, x), width=10, height=3)
    par(mfrow=c(1, nrow(sum.values)), mar=c(1.1, 2.1, 2.1, 1.1))
    for (i in seq_len(nrow(sum.values))) { 
        plot(jitter(ifelse(grouping=="lif", -0.5, 0.5), amount=0.4), sum.values[i,], col=ifelse(grouping=="lif", "black", "grey"), pch=16,
             xaxt="n", ylab="", xlab="", main=sum.names[i], cex.main=ifelse(missing.sum[i], 0.7, 1))
    }
    dev.off()

#    # Not particularly interested in looking at those *not* in the top for summation.   
#    set.seed(1000)
#    pdf(sprintf("top_%s_%s_raw.pdf", method, comp), width=10, height=3)
#    par(mfrow=c(1, nrow(sum.values)), mar=c(1.1, 2.1, 2.1, 1.1))
#    for (i in seq_len(nrow(sum.values))) { 
#        plot(jitter(ifelse(grouping=="lif", -0.5, 0.5), amount=0.4), raw.values[i,], col=ifelse(grouping=="lif", "black", "grey"), pch=16,
#             xaxt="n", ylab="", xlab="", main=raw.names[i], cex.main=ifelse(missing.raw[i], 0.7, 1))
#    }
#    dev.off()
}
