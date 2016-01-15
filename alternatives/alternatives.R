##################################################
# Setting up the simulation parameters.

set.seed(10000)

conditions <- rep(c("A", "B"), each=3)
nplates <- length(conditions)
ngenes <- 10000
ncells <- 50

all.means <- exp(runif(ngenes, 3, 8))
disp <- 0.5 + 100/all.means
plate.var <- 0.5

cell.grouping <- factor(rep(conditions, each=ncells))
alphas <- c(0.01, 0.002, 0.05)

plate.of.origin <- factor(rep(seq_len(nplates), each=ncells))

##################################################
# Generating the result file.

log.file <- "temp.txt"
if (file.exists(log.file)) { stop("existing file for results") }

save.fun <- function(label, mode, x, log.file) {
    discard <- is.na(x)
    totes <- sum(!discard)
    for (alpha in alphas) {
        # Protect against NAs (ID'd as outliers by DESeq2).
        write.table(data.frame(label, mode, alpha, sum(x<=alpha & !discard)/totes), file=log.file, 
                    sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

##################################################
# Running across all options.

suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(DESeq2))
suppressPackageStartupMessages(require(lme4))

for (it in 1:10) {
    mean.per.plate <- all.means * exp(matrix(rnorm(nplates*ngenes, -plate.var/2, sqrt(plate.var)), ncol=nplates))
    counts <- matrix(0L, nrow=ngenes, ncol=ncells*nplates)
    for (x in seq_len(nplates)) {
        counts[,seq_len(ncells)+ncells*(x-1L)] <- rnbinom(ngenes*ncells, mu=mean.per.plate[,x], size=1/disp)
    }

    for (mode in c("Blocked", "Averaged")) {
        if (mode=="Blocked") {
            targets <- DataFrame(Grouping=cell.grouping, Plate=plate.of.origin)
            design <- model.matrix(~Plate + Grouping, data=targets)
            design <- design[,-4]
            con <- c(0,0,0,0,0,1)
        } else if (mode=="Averaged") { 
            targets <- DataFrame(Plate=plate.of.origin)
            design <- model.matrix(~0 + Plate, data=targets)
            con <- c(1,1,1,-1,-1,-1)/3
        } 

        # edgeR QL:
        y <- DGEList(counts)
        y <- calcNormFactors(y)
        y <- estimateDisp(y, design, prior.df=0)
        fit <- glmQLFit(y, design, robust=TRUE)
        res <- glmQLFTest(fit, contrast=con)
        save.fun("edgeR (QL)", mode, res$table$PValue, log.file)
    }

    if (it==1L) { 
        pvalues <- rep(NA, ngenes)
        my.data <- data.frame(Counts=0L, Group=cell.grouping, Plate=plate.of.origin)
        for (g in seq_len(ngenes)) { 
            my.data$Counts <- counts[g,]
            try({
                full.fit <- glmer.nb(Counts ~ Group + (1|Plate), data=my.data)
                null.fit <- glmer.nb(Counts ~ (1|Plate), data=my.data)
                tested <- anova(full.fit, null.fit, test = "Chisq")
                # see ?glmer, http://glmm.wikidot.com/faq; best option other than bootstrapping (which would be bothersome for so many genes).
                pvalues[g] <- tested$Pr[2] 
            })
        }
        save.fun("lme4", "Mixed", pvalues, log.file)
    }
}

##################################################
#
#set.seed(2000)
#fc <- 1
#first.set <- 1:100 # Takes too long with 1000 genes.
#second.set <- first.set + 100
#
#mean.per.plate <- all.means * exp(matrix(rnorm(nplates*ngenes, -plate.var/2, sqrt(plate.var)), ncol=nplates))
#mean.per.plate[first.set,conditions=="A"] <- mean.per.plate[first.set,conditions=="A"]*fc
#mean.per.plate[second.set,conditions=="B"] <- mean.per.plate[second.set,conditions=="B"]*fc
#
#counts <- matrix(0L, nrow=ngenes, ncol=ncells*nplates)
#for (x in seq_len(nplates)) {
#    counts[,seq_len(ncells)+ncells*(x-1L)] <- rnbinom(ngenes*ncells, mu=mean.per.plate[,x], size=1/disp)
#}
#
## Mixed modelling, *once*, because I have to iterate across genes
#
## Compare to p-values from edgeR with summation.
#
#counts <- sumTechReps(counts, plate.of.origin)
#design <- model.matrix(~conditions)
#y <- DGEList(counts)
#y <- calcNormFactors(y)
#y <- estimateDisp(y, design, prior.df=0)
#fit <- glmQLFit(y, design, robust=TRUE)
#res <- glmQLFTest(fit)
#
#save(pvalues, res, file="MixedvSummed.Rda")
#
##################################################
# End.

sessionInfo()
