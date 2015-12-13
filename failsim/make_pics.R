# This makes some pretty pictures.

for (type in c("raw", "sum")) { 
    stuff <- read.table(sprintf("results_%s.txt", type), sep="\t", colClasses="character")
    threshold <- "0.01"

    for (pv in c("with", "without")) { 
        chosen <- stuff[,2]==ifelse(pv=="with", "0.5", "0") & stuff[,3]==threshold
        out <- split(log10(as.numeric(stuff[chosen,4])), stuff[chosen,1])

        all.means <- sapply(out, FUN=mean)
        all.se <- sqrt(sapply(out, FUN=var)/lengths(out))

        setEPS()
        postscript(sprintf("%s_%s.eps", pv, type))
        par(mar=c(6.2, 4.1, 4.1, 2.1))
        xbounds <- c(0.5, length(all.means)+0.5)
        ybounds <- c(-3, 0)
        
        plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="Observed type I error rate", ylim=ybounds, xlim=xbounds, cex.lab=1.4,
                main=sprintf("%s plate effects", paste0(toupper(substring(pv, 1,1)), substring(pv, 2))), cex.main=1.4)
        xticks <- seq_along(out)
        abline(v=xticks, col="grey", lwd=1.5, lty=3)
        points(xticks, all.means, pch=16, cex=2)

        # Adding the axes.
        axis(1, xticks, labels=FALSE)
        text(x=xticks, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=names(all.means), srt=45, adj=1, xpd=TRUE, cex=1.4)

        yticks <- seq(ybounds[1], ybounds[2], by=1)
        axis(2, yticks, labels=10^yticks, cex.axis=1.2)

        # Adding the error bars.
        segments(xticks, all.means, xticks, all.means+all.se)
        segments(xticks-0.2, all.means+all.se, xticks+0.2, all.means+all.se, lwd=2)

        # Adding the threshold.
        abline(h=log10(as.numeric(0.01)), col="grey50", lwd=2, lty=2)


        dev.off()
    }
}
