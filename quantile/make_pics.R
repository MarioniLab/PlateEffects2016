# This makes some pretty pictures.

stuff <- read.table("results.txt", sep="\t", colClasses="character")
threshold <- "0.01"
plot.order <- as.character(c(1:4, -(1:4)))
xpos <- rep(1:8*1.2, each=2) + rep(c(-0.2, 0.2), 8) # + rep(c(0, 0.4), each=8)

for (pv in c("with", "without")) { 
    chosen <- stuff[,3]==ifelse(pv=="with", "0.5", "0") & stuff[,5]==threshold 
    chosen.sum <- chosen & stuff[,4]=="sum"
    out.sum <- split(log10(as.numeric(stuff[chosen.sum,6])), stuff[chosen.sum,2])
    chosen.quantile <- chosen & stuff[,4]=="quantile"
    out.quantile <- split(log10(as.numeric(stuff[chosen.quantile,6])), stuff[chosen.quantile,2])

    sum.means <- sapply(out.sum, FUN=mean)
    sum.se <- sqrt(sapply(out.sum, FUN=var)/lengths(out.sum))
    quantile.means <- sapply(out.quantile, FUN=mean)
    quantile.se <- sqrt(sapply(out.quantile, FUN=var)/lengths(out.quantile))
    
    stopifnot(identical(names(sum.means), names(quantile.means)))
    m <- match(plot.order, names(sum.means))
    all.means <- as.vector(rbind(sum.means, quantile.means)[,m])
    all.se <- as.vector(rbind(sum.se, quantile.se)[,m])

    is.quant <- rep(c(FALSE, TRUE), 8)
    cols <- ifelse(is.quant, "grey50", "black")

    setEPS()
    postscript(sprintf("%s.eps", pv))
    xbounds <- range(xpos)
    ybounds <- c(-3, -1)
    plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="Observed type I error rate", ylim=ybounds, xlim=xbounds, cex.lab=1.4)

    # Adding the alternating shading.
    xmid <- (xpos[1:8*2-1] + xpos[1:8*2])/2
    for (shade in seq_along(xmid)) {
        if (shade%%2==1L) {
            rect(xmid[shade]-0.6, -4, xmid[shade]+0.6, 0, col=rgb(240, 250, 250, max=255), border=NA)
        }
    }
    box()
    points(xpos, all.means, pch=ifelse(is.quant, 18, 16), cex=1.5, col=cols)

    # Adding the axes.
    axis(1, (xpos[1:8*2]+xpos[1:8*2-1])/2, labels=c(1:4, 1:4), cex.axis=1.2)
    h <- -3.3
    segments(xpos[9], h, xpos[16], h, xpd=TRUE, lwd=2)
    text((xpos[9]+xpos[16])/2, h-0.1, "Swapped", xpd=TRUE, cex=1.2) 

    yticks <- seq(ybounds[1], ybounds[2], by=1)
    axis(2, yticks, labels=10^yticks, cex.axis=1.2)

    # Adding the error bars.
    segments(xpos, all.means, xpos, all.means+all.se, col=cols)
    segments(xpos-0.1, all.means+all.se, xpos+0.1, all.means+all.se, lwd=2, col=cols)

    # Adding the threshold.
    abline(h=log10(as.numeric(0.01)), col="red", lwd=2, lty=2)

    dev.off()
}

