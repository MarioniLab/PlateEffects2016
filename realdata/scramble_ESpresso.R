#######################################################
# Loading the ESpresso counts and experimental design.

curdir <- getwd()
setwd("../reference")
source("ESpresso.R")
setwd(curdir)

of.interest <- rowMeans(all.counts) >= 1
all.counts <- all.counts[of.interest,]
methods.to.use <- c("QLedgeR", "DESeq2", "voom")

dir.create("ESpresso", showWarning=FALSE)
outputfile <- file.path("ESpresso/scrambled.txt")
if (file.exists(outputfile)) { unlink(outputfile) }

for (extra.term in c('lif', 'a2i', '2i')) {
    targets$Interaction <- as.integer(extra.term == targets$Serum & targets$Batch == '3')

    for (type in c("raw", "sum")) { 
        my.env <- new.env()
        my.env$sample.formula <- ~Batch + Serum + Interaction
        my.env$drop.coefficient <- 5

        if (type=="raw") { 
            my.env$counts <- all.counts
            my.env$sample.data <- targets
            my.env$normtype <- "deconvolution"
            my.env$clusters <- by.group
        } else {
            my.env$counts <- sumTechReps(all.counts, targets$Plate)
            my.env$sample.data <- targets[match(colnames(my.env$counts), targets$Plate),]
            my.env$norm.type <- "deseq"
        }
   
        sys.source("../reference/de_analysis.R", envir=my.env)
        for (method in names(my.env$obtained)) {
            cur.pval <- my.env$obtained[[method]]
            is.okay <- !is.na(cur.pval)
            write.table(file=outputfile, data.frame(method, extra.term, type, sum(cur.pval <= 0.01 & is.okay)/sum(is.okay)), 
                  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
        }
    }
}

#######################################################

sessionInfo()

#######################################################
# End.
