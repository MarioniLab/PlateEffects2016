for (dataset in c("ESpresso")) { 

#######################################################
# Loading the counts and experimental design.

curdir <- getwd()
setwd("../reference")
if (dataset=="ESpresso") { 
    source("ESpresso.R")
    coefs <- c(3, 4)
    sample.formula <- ~Batch + Serum
}
setwd(curdir)

dir.create(dataset, showWarning=FALSE)
of.interest <- rowMeans(all.counts) >= 1
all.counts <- all.counts[of.interest,]
methods.to.use <- c("QLedgeR", "DESeq2", "voom")

#######################################################
# Running through the real data analyses. 

for (type in c("raw", "sum")) { 
    my.env <- new.env()
    my.env$sample.formula <- sample.formula
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
    for (i in coefs) { 
        qres <- glmQLFTest(my.env$qfit, coef=i) 
        fhandle <- gzfile(file.path(dataset, paste0("edgeR_", i, "_", type, ".tsv.gz")), open="wb")
        write.table(topTags(qres, n=Inf, sort.by="none")$table, file=fhandle,
                    col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
    }

    tmp.y <- my.env$y
    tmp.design <- my.env$design
    save(tmp.y, tmp.design, file=file.path(dataset, paste0("objects_", type, ".Rda"))) # The DGEGLM objects are huge, so we'll go without.     

    # Analyzing with voom.
    for (i in coefs) { 
        vres <- topTable(my.env$vfit, n=Inf, sort.by="none", coef=i)
        fhandle <- gzfile(file.path(dataset, paste0("voom_", i, "_", type, ".tsv.gz")), open="wb")
        write.table(vres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
    }

    # Analyzing with DESeq2.
    for (i in coefs) {
        convec <- integer(ncol(my.env$design))
        convec[i] <- 1
        dres <- results(my.env$dds, convec)

        fhandle <- gzfile(file.path(dataset, paste0("DESeq2_", i, "_", type, ".tsv.gz")), open="wb")
        write.table(dres, file=fhandle, col.names=NA, sep="\t", quote=FALSE)
        close(fhandle)
    }
}

#######################################################
}

sessionInfo()

#######################################################
# End.
