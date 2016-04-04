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
hist(log10(by.cell$samples$lib.size), xlab=expression(Log[10]~"library size"), cex.lab=1.4, cex.axis=1.2, ylab="Number of cells", main="", col="grey80")
dev.off()

# Estimating the NB dispersion (assuming sufficient residual d.f. to estimate the dispersion without EB shrinkage).
plateX <- model.matrix(~by.plate)
by.cell <- estimateDisp(by.cell, plateX, prior.df=0, trend='none')

# Estimating the log-overall mean.
centered.off <- getOffset(by.cell)
centered.off <- centered.off - mean(centered.off)
logmeans <- mglmOneGroup(by.cell$counts, offset=centered.off, dispersion=by.cell$tagwise.dispersion)

pdf(file.path(picdir, "avecounts.pdf"))
hist(logmeans/log(10), xlab=expression(Log[10]~"average count"), cex.lab=1.4, cex.axis=1.2, ylab="Number of cells", main="", col="grey80")
dev.off()

# Refitting a mean-dispersion trend to the logs.
ldisp <- log10(by.cell$tagwise.dispersion)
lfit <- loessFit(x=logmeans, y=ldisp, span=0.1)

pdf(file.path(picdir, "celldisp.pdf"))
par(mar=c(5.1,4.5,2.1,2.1), cex.axis=1.2, cex.lab=1.4)
med.ldisp <- median(ldisp)
mad.ldisp <- mad(ldisp)
dontshow <- ldisp <= med.ldisp - 5*mad.ldisp | ldisp >= med.ldisp + 5*mad.ldisp
smoothScatter(logmeans[!dontshow]/log(10), ldisp[!dontshow], xlab=expression(Log[10]~"average count"), ylab=expression(Log[10]~"NB dispersion"))
o <- order(logmeans)
lines(logmeans[o]/log(10), lfit$fitted[o], col="red", lwd=2)
dev.off()

# Estimating the dispersion for the summed counts.
summed <- sumTechReps(by.cell$counts, by.plate)
first.in.each <- !duplicated(by.plate)
design <- refdesign[first.in.each,]

by.sum <- DGEList(summed)
by.sum <- calcNormFactors(by.sum)
by.sum <- estimateDisp(by.sum, design)

pdf(file.path(picdir, "sumdisp.pdf"))
ldisp <- log10(by.sum$tagwise.dispersion)
lfit <- loessFit(x=logmeans, y=ldisp, span=0.1)
par(mar=c(5.1,4.5,2.1,2.1), cex.axis=1.2, cex.lab=1.4)
smoothScatter(logmeans/log(10), ldisp, xlab=expression(Log[10]~"average count"), ylab=expression(Log[10]~"NB dispersion"))
o <- order(logmeans)
lines(logmeans[o]/log(10), lfit$fitted[o], col="red", lwd=2)
dev.off()

# Plotting the absolute differences, which represents the reciprocal of the Gamma shape parameter for a NB-Gamma mixture model.
pdf(file.path(picdir, "rate.pdf"))
leftover <- by.sum$trended - 10^lfit$fitted/mean(table(by.plate))
par(mar=c(5.1,4.5,2.1,2.1), cex.axis=1.2, cex.lab=1.4)
plot(logmeans/log(10), log10(leftover), xlab=expression(Log[10]~"average count"), ylab=expression("\u2013"~Log[10]~"Gamma shape"))
dev.off()

##############################################
# Repeating the estimation of the dispersion with ZINB models.

zinb.prop <- numeric(nrow(by.cell))
zinb.disp <- by.cell$tagwise.dispersion
zinb.mean <- exp(logmeans)

library(pscl)
for (i in which(rowSums(by.cell$counts==0)>0)) { 
    try({
        zfit <- zeroinfl(by.cell$count[i,] ~ 0 + by.plate | 1, dist="negbin")
        zinb.mean[i] <- mean(exp(zfit$coefficients$count))
        zinb.prop[i] <- zfit$coefficients$zero
        zinb.disp[i] <- 1/zfit$theta
    })
}
zinb.prop <- exp(zinb.prop)/(1+exp(zinb.prop))

##############################################
# Constructing sampling functions.

COUNT_FUN_GEN <- function(cell.mu, cell.disp, plate.shape, relative.size, zi.mu, zi.disp, zi.prop) {
    relative.size <- relative.size/mean(relative.size) # mean centering. 
    if (length(unique(c(length(cell.mu), length(cell.disp), length(plate.shape), 
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
        cur.shape <- plate.shape[chosen] * mod.shape
        plate.means <- matrix(rgamma(ngenes*nplates, shape=cur.shape, scale=cur.mean/cur.shape), nrow=ngenes, ncol=nplates)

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

COUNT_GEN <- COUNT_FUN_GEN(cell.mu=exp(logmeans), cell.disp=by.cell$tagwise.dispersion, plate.shape=1/leftover,
                           relative.size=by.cell$samples$lib.size, zi.mu=zinb.mean, zi.disp=zinb.disp, zi.prop=zinb.prop)
saveRDS(COUNT_GEN, file=file.path(picdir, "function.rds"))

