# This makes plots about the variability of the top set of genes for the Kolod data,
# using cells in the first batch for easy visualization (not Scialdone, as there aren't any DE genes).

curdir <- getwd()
setwd("../reference")
source("ESpresso.R")
coefs <- c(3, 4)
setwd(curdir)

for (x in coefs) {
    current <- targets$Serum %in% c(ifelse(x==3, '2i', 'a2i'), "lif") 
    data <- all.counts[,current]
    grouping <- targets$Serum[current]
    require(edgeR)
    cpms <- cpm(data, log=TRUE)

    # Removing the batch effect, but replacing zeroes with NAs.
    cpms <- removeBatchEffect(cpms, batch=targets$Batch[current])
    cpms[data==0L] <- NA

    method <- "edgeR"
    cur.field <- switch(method, edgeR="PValue", voom="P.Value", DESeq2="pvalue")
    sum.res <- read.table(sprintf("ESpresso/%s_%s_sum.tsv.gz", method, x), header=TRUE, row.names=1)
    sum.res <- sum.res[order(sum.res[[cur.field]]),]
    sum.top <- rownames(sum.res)[1:10]
    sum.values <- cpms[sum.top,]
    
    require(org.Mm.eg.db)
    sum.names <- select(org.Mm.eg.db, keys=rownames(sum.values), columns="SYMBOL", keytype="ENSEMBL")
    sum.names <- sum.names$SYMBOL[match(rownames(sum.values), sum.names$ENSEMBL)]
    missing.sum <- is.na(sum.names)
    sum.names[missing.sum] <- rownames(sum.values)[missing.sum]
    
    set.seed(1000)
    pdf(sprintf("ESpresso/top_%s_%s_sum.pdf", method, x), width=10, height=3)
    par(mfrow=c(1, nrow(sum.values)), mar=c(1.1, 2.1, 2.1, 1.1))
    for (i in seq_len(nrow(sum.values))) { 
        is.lif <- grouping=="lif"
        ybounds <- range(sum.values[i,], na.rm=TRUE)
        plot(0, 0, type="n", xaxt="n", ylab="", xlab="", main=sum.names[i], cex.main=ifelse(missing.sum[i], 0.6, 1), xlim=c(-0.9, 0.9), ylim=ybounds)
        is.zero <- is.na(sum.values[i,])
        lif.prop.zero <- sum(is.zero[is.lif])/sum(is.lif)
        nlif.prop.zero <- sum(is.zero[!is.lif])/sum(!is.lif)
        rect(-1, ybounds[1], -0.05, ybounds[1]-10, col=rgb(1, 0, 0, lif.prop.zero), border=NA)
        rect(0.05, ybounds[1], 1, ybounds[1]-10, col=rgb(1, 0, 0, nlif.prop.zero), border=NA)
        points(jitter(ifelse(is.lif, -0.5, 0.5), amount=0.4), sum.values[i,], col=ifelse(is.lif, "black", "grey"), pch=16)
        box()
    }
    dev.off()
}
