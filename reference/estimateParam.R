# Given a design matrix and some counts, we estimate various parameters for simulation.
# - average count.
# - spread of library sizes.
# - dispersion trend for counts (NB).
# - dispersion trend for counts (ZINB).
# - variance trend for log-summed counts.
# - also values from duplicateCorrelation.

all.counts
refdesign
by.plate
picdir
dir.create(picdir, showWarnings=FALSE)

##############################################

# Simple count loading and filtering.
library(edgeR)
by.cell <- DGEList(all.counts)
by.cell <- by.cell[rowMeans(by.cell$counts) >= 1L,]

# Normalizing with the deconvolution methods.
library(scran)
if (!exists("by.group")) { by.group <- quickCluster(by.cell$counts) }
sf <- computeSumFactors(by.cell$counts, cluster=by.group)
by.cell$samples$norm.factors <- sf/by.cell$samples$lib.size
by.cell$samples$norm.factors <- by.cell$samples$norm.factors/exp(mean(log(by.cell$samples$norm.factors))) # Mean-centering for convenience.

pdf(file.path(picdir, "libsizes.pdf"))
hist(log2(sf), xlab=expression(Log[2]~"size factor"), cex.lab=1.4, cex.axis=1.2, ylab="Number of cells", main="", col="grey80")
dev.off()

# Estimating the NB dispersion (assuming sufficient residual d.f. to estimate the dispersion without EB shrinkage).
plateX <- model.matrix(~by.plate)
by.cell <- estimateDisp(by.cell, plateX, prior.df=0, trend='none')

# Estimating the log-overall mean.
centered.off <- getOffset(by.cell)
centered.off <- centered.off - mean(centered.off)
logmeans <- mglmOneGroup(by.cell$counts, offset=centered.off, dispersion=by.cell$tagwise.dispersion)

pdf(file.path(picdir, "avecounts.pdf"))
hist(logmeans/log(2), xlab=expression(Log[2]~"average count"), cex.lab=1.4, cex.axis=1.2, ylab="Number of cells", main="", col="grey80")
dev.off()

pdf(file.path(picdir, "celldisp.pdf"))
par(mar=c(5.1,4.5,2.1,2.1), cex.axis=1.2, cex.lab=1.4)
ldisp <- log2(by.cell$tagwise.dispersion)
med.ldisp <- median(ldisp)
mad.ldisp <- mad(ldisp)
dontshow <- ldisp <= med.ldisp - 5*mad.ldisp | ldisp >= med.ldisp + 5*mad.ldisp
smoothScatter(logmeans[!dontshow]/log(2), ldisp[!dontshow], xlab=expression(Log[2]~"average count"), ylab=expression(Log[2]~"NB dispersion"))
dev.off()

##############################################

# Estimating the plate effect variance.
library(parallel)
library(lme4)
collected <- mclapply(seq_len(nrow(by.cell)), function(i) {
    output <- NA_real_         
    try({ 
        out <- glmer(Counts ~ 0 + refdesign + (1|Plate) + offset(log(sf)), 
                     data=data.frame(Counts=as.integer(all.counts[i,]), Group=by.group, Plate=by.plate), 
                     family=negative.binomial(1/by.cell$tagwise[i]))
        output <- unlist(VarCorr(out))
    })
    return(output)
}, mc.cores=4)

sigma2 <- mean(unlist(collected), na.rm=TRUE)

pdf(file.path(picdir, "platevar.pdf")) 
collected2 <- unlist(collected)
collected2 <- collected2[collected2 <= 3]
hist(collected2, xlab=expression("Estimated"~sigma^2), cex.lab=1.4, cex.axis=1.2, ylab="Number of genes", main="", col="grey80", breaks=20)
abline(v=sigma2, col="red", lwd=2, lty=2)
dev.off()

## Test code, to see whether this strategy generally works:
#true.mean <- exp(rnorm(100, 10, 2))
#assignments <- rep(seq_len(100), each=10)
#actual.count <- rnbinom(1000, mu=true.mean[assignments], size=2)
#library(lme4)
#out <- glmer(Count ~ (1|Assigned), data=data.frame(Count=actual.count, Assigned=assignments), family=negative.binomial(2))
#unlist(VarCorr(out))

##############################################
# Repeating the estimation of the dispersion with ZINB models.

zinb.prop <- rep(-Inf, nrow(by.cell))
zinb.disp <- by.cell$tagwise.dispersion
zinb.mean <- exp(logmeans)

library(pscl)
for (i in which(rowSums(by.cell$counts==0)>0)) { 
    try({
        zfit <- zeroinfl(by.cell$count[i,] ~ 0 + by.plate | 1, dist="negbin", offset=log(sf))
        zinb.mean[i] <- mean(exp(zfit$coefficients$count))
        zinb.prop[i] <- zfit$coefficients$zero
        zinb.disp[i] <- 1/zfit$theta
    })
}
zinb.prop <- exp(zinb.prop)/(1+exp(zinb.prop))

##############################################
# Constructing sampling functions.

COUNT_FUN_GEN <- function(cell.mu, cell.disp, sigma2, relative.size, zi.mu, zi.disp, zi.prop) {
    relative.size <- relative.size/mean(relative.size) # mean centering. 
    sigma2 # evaluating the promise in this environment.
    if (length(unique(c(length(cell.mu), length(cell.disp), 
                        length(zi.mu), length(zi.disp), length(zi.prop))))!=1L) {
        stop("vector lengths must be the same")
    }

    function(plates, ngenes, nde=0, ingroup, fc, zinb=FALSE, mod.lib=1, mod.shape=1) {
        chosen <- sample(length(cell.mu), ngenes, replace=TRUE)
        if (!zinb) {
            cur.mean <- cell.mu[chosen]
            cur.disp <- cell.disp[chosen]
        } else {
            cur.mean <- zi.mu[chosen]
            cur.disp <- zi.disp[chosen]
        }

        plates <- as.factor(plates)
        nplates <- nlevels(plates)
        sigma2 <- sigma2/mod.shape
        plate.means <- cur.mean * matrix(exp(rnorm(ngenes*nplates, mean=-sigma2/2, sd=sqrt(sigma2))), nrow=ngenes, ncol=nplates)

        # Adding DE, if so desired (identifies the plates within each group).
        if (nde > 0) {
            m <- match(ingroup, levels(plates))
            chosen.de <- sample(ngenes, nde)
            is.up <- rep(c(TRUE, FALSE), length.out=nde)
            fc <- sqrt(fc)
            plate.means[chosen.de[is.up],m] <- plate.means[chosen.de[is.up],m] * fc
            plate.means[chosen.de[is.up],-m] <- plate.means[chosen.de[is.up],-m]/fc
            plate.means[chosen.de[!is.up],m] <- plate.means[chosen.de[!is.up],m]/fc
            plate.means[chosen.de[!is.up],-m] <- plate.means[chosen.de[!is.up],-m]*fc
        } else {
            chosen.de <- integer(0)
        }
        is.de <- logical(ngenes)
        is.de[chosen.de] <- TRUE

        # Expanding and adding variable library sizes.
        expanded <- as.integer(plates)
        cell.means <- plate.means[,expanded]
        ncells <- length(plates)
        chosen.lib.mult <- sample(relative.size, ncells, replace=TRUE) * mod.lib
        cell.means <- t(t(cell.means)*chosen.lib.mult)

        # Simulating counts.
        counts <- matrix(rnbinom(ncells*ngenes, mu=cell.means, size=1/cur.disp), ncol=ncells, nrow=ngenes)
        if (zinb) {
            is.zero <- matrix(rbinom(ncells*ngenes, 1, zi.prop[chosen]), ncol=ncells, nrow=ngenes)==1L
            counts[is.zero] <- 0
        }
        return(list(counts=counts, is.de=is.de))
    }
}

COUNT_GEN <- COUNT_FUN_GEN(cell.mu=exp(logmeans), cell.disp=by.cell$tagwise.dispersion, sigma2=sigma2,
                           relative.size=by.cell$samples$lib.size, zi.mu=zinb.mean, zi.disp=zinb.disp, zi.prop=zinb.prop)
saveRDS(COUNT_GEN, file=file.path(picdir, "function.rds"))

