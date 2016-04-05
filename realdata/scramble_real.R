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

of.interest <- rowMeans(all.counts) >= 1
all.counts <- all.counts[of.interest,]
methods.to.use <- c("QLedgeR", "DESeq2", "voom")

dir.create(dataset, showWarning=FALSE)
outputfile <- file.path(dataset, "scrambled.txt")
if (file.exists(outputfile)) { unlink(outputfile) }

# Running through sum and row, while shuffling targets.

for (scenario in seq_along(coefs)) {
    new.targets <- targets
    drop.coefficient <- coefs[scenario]

    if (dataset=="ESpresso") {
        if (drop.coefficient==3) { shuffle <- '2i' } else { shuffle <- 'a2i' }
        cur.other <- which(shuffle == targets$Serum & targets$Batch == '3')
        cur.lif <- which('lif' == targets$Serum & targets$Batch == '3')
        new.targets$Serum[cur.other] <- 'lif'
        new.targets$Serum[cur.lif] <- shuffle
    }

    for (type in c("raw", "sum")) { 
        my.env <- new.env()
        my.env$sample.formula <- sample.formula
        my.env$drop.coefficient <- drop.coefficient

        if (type=="raw") { 
            my.env$counts <- all.counts
            my.env$sample.data <- new.targets
            my.env$normtype <- "deconvolution"
            my.env$clusters <- by.group
        } else {
            my.env$counts <- sumTechReps(all.counts, targets$Plate)
            my.env$sample.data <- new.targets[match(colnames(my.env$counts), targets$Plate),]
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
