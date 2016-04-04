##################################################

compute.roc <- function(x, is.de) {
    x[is.na(x)] <- 1
    fn <- x[is.de]
    tn <- x[!is.de]
    return(findInterval(sort(fn), sort(tn))/length(tn))    
}

compute.fdr <- function(x, is.de, threshold=0.05) {
    sig <- p.adjust(x, method="BH") <= threshold
    if (any(sig, na.rm=TRUE)) { 
        return(sum(sig & !is.de, na.rm=TRUE)/sum(sig, na.rm=TRUE))
    } else {
        return(0)
    }
}

##################################################
# Setting up the simulation parameters.

set.seed(10000)
conditions <- rep(c("A", "B"), each=3)
ngenes <- 10000
nplates <- length(conditions)

args <- commandArgs(trailingOnly = TRUE)
seed <- 10000
scenario <- 1L
for (i in seq_along(args)){
        eval(parse(text=args[[i]]))
}
scenario <- as.integer(scenario)
print(scenario)
print(seed)

# Input/output files.
print(indir)
COUNT_GEN <- readRDS(file.path(indir, "function.rds"))

print(outdir)
dir.create(outdir, showWarning=FALSE)
log.raw <- file.path(outdir, sprintf("raw_%i.txt", scenario))
log.sum <- file.path(outdir, sprintf("sum_%i.txt", scenario))

##################################################
# Setting up scenario-specific parameters.

mod.shape <- 1 # Defaults.
de.fc <- 3
nde <- 2000
ncell.range <- c(50, 100)

if (scenario==2L) {
    mod.shape <- 1e6 # No plate effect.
} else if (scenario==3L) {
    de.fc <- 6 # larger fold change.
} else if (scenario==4L) {
    nde <- 4000 # More DE genes.
}

##################################################
# Actually running it.

methods.to.use <- c("voom", "QLedgeR", "DESeq2")
collected <- list()
for (method in methods.to.use) { 
    collected[[method]] <- list(raw=list(), sum=list())
}

for (it in 1:10) {
    ncells <- round(runif(nplates, ncell.range[1], ncell.range[2]))
    plate.of.origin <- rep(seq_len(nplates), ncells)
    cell.grouping <- factor(rep(conditions, ncells))
    plate.grouping <- factor(conditions)

    simulated <- COUNT_GEN(plate.of.origin, ngenes, mod.shape=mod.shape, nde=nde, ingroup=which(plate.grouping=="A"), fc=de.fc)
    all.counts <- simulated$counts
    true.de <- simulated$is.de

    for (mode in c("raw", "sum")) {
        my.env <- new.env()
        my.env$normtype <- "libsize"
        if (mode=="raw") {
            my.env$counts <- all.counts
            my.env$sample.data <- data.frame(group=cell.grouping)
            my.env$sample.formula <- ~group
            my.env$drop.coefficient <- 2
            log.file <- log.raw            
        } else {
            my.env$counts <- sumTechReps(all.counts, plate.of.origin)
            my.env$sample.data <- data.frame(group=plate.grouping)
            my.env$sample.formula <- ~group
            my.env$drop.coefficient <- 2
            log.file <- log.sum 
        }

        # Saving the ROC stats.
        sys.source("../reference/de_analysis.R", envir=my.env)
        for (method in names(my.env$obtained)) { 
            collected[[method]][[mode]][[it]] <- compute.roc(my.env$obtained[[method]], true.de)
        }

        # Also saving the FDRs.
        for (method in names(my.env$obtained)) { 
            write.table(data.frame(method, scenario, compute.fdr(my.env$obtained[[method]], true.de)),
                        file=log.file, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
    }
}

# Writing the ROC stats to file.
saveRDS(file=file.path(outdir, paste0("roc_", scenario, ".rds")), collected)

##################################################
# End.

sessionInfo()
