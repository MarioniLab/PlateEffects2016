host.dir <- file.path("../../DataFiles/ESpresso")
all.files <- c("2i_bulk_1.gene.count", "2i_bulk_2.gene.count", "a2i_bulk_2a.gene.count", "a2i_bulk_2b.gene.count", "serum_bulk_1.gene.count", "serum_bulk_2.gene.count")
combined <- lapply(file.path(host.dir, all.files), read.table, row.names=1)

# Checking all the rownames are the same.
refnames <- rownames(combined[[1]])
for (x in seq_along(combined)) {
    stopifnot(identical(rownames(combined[[x]]), refnames))
}

# Merging technical replicates.
biolibs <- sub("2[ab].gene.count", "2.gene.count", all.files)
library(edgeR)
summed <- sumTechReps(do.call(data.frame, combined), biolibs)
rownames(summed) <- refnames

# Retaining all genes used in the single-cell analysis.
tmp.env <- new.env()
load("ESpresso/objects_sum.Rda", envir=tmp.env)
to.keep <- rownames(summed) %in% rownames(tmp.env$tmp.y)
summed <- summed[to.keep,]

# Setting up the analysis.
tdat <- strsplit(sub(".gene.count$", "", colnames(summed)), split="_")
serum <- sapply(tdat, "[", i=1)
batch <- sapply(tdat, "[", i=3)

my.env <- new.env()
my.env$counts <- summed
my.env$sample.formula <- ~ Batch + Serum
my.env$sample.data <- data.frame(Batch=batch, Serum=factor(serum, levels=c("serum", '2i', 'a2i')))

coefs <- c(3,4)
my.env$drop.coefficient <- coefs[1]
my.env$normtype <- "deseq"
methods.to.use <- c("QLedgeR", "DESeq2", "voom")
sys.source("../reference/de_analysis.R", envir=my.env)

# Analyzing with edgeR. 
for (i in coefs) { 
    qres <- glmQLFTest(my.env$qfit, coef=i) 
    fhandle <- gzfile(paste0("ESpresso/edgeR_", i, "_bulk.tsv.gz"), open="wb")
    write.table(topTags(qres, n=Inf, sort.by="none")$table, file=fhandle,
                col.names=NA, sep="\t", quote=FALSE)
    close(fhandle)
}

# Analyzing with voom.
for (i in coefs) { 
    vres <- topTable(my.env$vfit, n=Inf, sort.by="none", coef=i)
    fhandle <- gzfile(paste0("ESpresso/voom_", i, "_bulk.tsv.gz"), open="wb")
    write.table(vres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
    close(fhandle)
}

# Analyzing with DESeq2.
for (i in coefs) {
    convec <- integer(ncol(my.env$design))
    convec[i] <- 1
    dres <- results(my.env$dds, convec)
    
    fhandle <- gzfile(paste0("ESpresso/DESeq2_", i, "_bulk.tsv.gz"), open="wb")
    write.table(dres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
    close(fhandle)
}

