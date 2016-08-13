##################################################
# Generating the result file.

alphas <- c(0.01, 0.002, 0.05)
save.fun <- function(label, pv, x, log.file) {
    discard <- is.na(x)
    totes <- sum(!discard)
    for (alpha in alphas) {
        # Protect against NAs (ID'd as outliers by DESeq2).
        write.table(data.frame(label, pv, alpha, sum(x<=alpha & !discard)/totes), file=log.file, 
                    sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

##################################################
# Setting up the simulation parameters.

ngenes <- 10000

args <- commandArgs(trailingOnly = TRUE)
seed <- 10000
scenario <- 1L
niter <- 10L
do.long <- TRUE
for (i in seq_along(args)){
    eval(parse(text=args[[i]]))
}
scenario <- as.integer(scenario)
print(scenario)
print(seed)
print(niter)
print(do.long)

# Input/output files.
print(indir)
COUNT_GEN <- readRDS(file.path(indir, "function.rds"))

print(outdir)
dir.create(outdir, showWarning=FALSE)
log.raw <- file.path(outdir, sprintf("raw_%i.txt", scenario))
if (file.exists(log.raw)) { unlink(log.raw) }
log.sum <- file.path(outdir, sprintf("sum_%i.txt", scenario))
if (file.exists(log.sum)) { unlink(log.sum) }

##################################################
# Setting up scenario-specific parameters.

mod.shape <- 1 # Defaults.
ncell.range <- c(50, 100)
mod.lib.fun <- function(n) 1
zinb <- FALSE
conditions <- rep(c("A", "B"), each=3)

if (scenario==2L) {
    mod.shape <- 2 # half the plate effect.
} else if (scenario==3L) {
    mod.shape <- 1e6 # no plate effect.
} else if (scenario==4L) {
    ncell.range <- c(20, 100) # increased range of cells.
} else if (scenario==5L) {
    mod.lib.fun <- function(n) 2^rnorm(n, mean=0, sd=0.5) # increased range of library sizes.
} else if (scenario==6L) {
    zinb <- TRUE # zero-inflated NB simulation.
} else if (scenario==7) {
    conditions <- rep(c("A", "B"), each=6)
} else if (scenario==8) {
    mod.shape <- rchisq(ngenes, df=5)/5 # Heteroscedastic plate effects across genes.
}

nplates <- length(conditions)

##################################################
# Actually running it.

set.seed(seed)
for (it in seq_len(niter)) {
    ncells <- round(runif(nplates, ncell.range[1], ncell.range[2]))
    plate.of.origin <- rep(seq_len(nplates), ncells)
    simulated <- COUNT_GEN(plate.of.origin, ngenes, mod.shape=mod.shape, mod.lib=mod.lib.fun(sum(ncells)), zinb=zinb)
    all.counts <- simulated$counts

    cell.grouping <- factor(rep(conditions, ncells))
    plate.grouping <- factor(conditions)

    for (mode in c("raw", "sum")) {
        my.env <- new.env()
        my.env$normtype <- "libsize"

        if (mode=="raw") {
            my.env$counts <- all.counts
            my.env$sample.data <- data.frame(group=cell.grouping)
            my.env$sample.formula <- ~group
            my.env$drop.coefficient <- 2

            log.file <- log.raw            
            methods.to.use <- c("edgeR", "voom", "voomcor", "QLedgeR", "DESeq2", "MAST")
            if (do.long && it==1L) { methods.to.use <- c(methods.to.use, "monocle", "glmer") }
        } else {
            my.env$counts <- sumTechReps(all.counts, plate.of.origin)
            my.env$sample.data <- data.frame(group=plate.grouping)
            my.env$sample.formula <- ~group
            my.env$drop.coefficient <- 2

            log.file <- log.sum
            methods.to.use <- c("voom", "QLedgeR", "DESeq2")
        }

        sys.source("../reference/de_analysis.R", envir=my.env)
        for (method in names(my.env$obtained)) { 
            save.fun(method, scenario, my.env$obtained[[method]], log.file)
        }
    }
}

##################################################
# End.

sessionInfo()
