cur.dir <- "ESpresso"
out.dir <- file.path(cur.dir, "results_failsim")

#########################################################################################
# Defining a function to extract mean and standard errors of the type I error rate.

threshold <- "0.01"
extractor <- function(fname) {
    out <- read.table(file.path(out.dir, fname), sep="\t", stringsAsFactors=FALSE)
    chosen <- out[,3]==threshold
    out <- out[chosen,]
    out[,4] <- pmax(out[,4], 1e-8) # for visualization purposes, to avoid undefined logs.
    out <- split(log10(out[,4]), out[,1])
    all.means <- sapply(out, FUN=mean)
    all.se <- sqrt(sapply(out, FUN=var)/lengths(out))
    return(list(mean=all.means, se=all.se))
}

renamed <- c(DESeq2="DESeq2", edgeR="edgeR", glmer="GLMM", MAST="MAST", monocle="Monocle", QLedgeR="QL edgeR", voom="voom", voomcor="voom\n+ cor")

#########################################################################################
# First, making the centrepiece for the raw results.

for (mode in 1:4) {
    devfun <- function(fname, ...) {
        setEPS()
        postscript(fname, ...)
    } 

    if (mode==1L) {
        # Main figure for raw results.
        raw1 <- extractor("raw_1.txt")
        raw3 <- extractor("raw_3.txt")
        
        stopifnot(identical(names(raw1$mean), names(raw3$mean)))
        all.means <- rbind(raw1$mean, raw3$mean)
        all.se <- rbind(raw1$se, raw3$se)
        pch <- 21
        inner.color <- c("black", "white")

        modes <- c('Plate effect', 'No plate effect')
        out.pic <- file.path(out.dir, "main_raw.eps")
        width <- 8
        height <- 7
    } else if (mode==2L) {
        # Supplementaries for raw results.
        raw2 <- extractor("raw_2.txt")
        raw4 <- extractor("raw_4.txt")
        raw5 <- extractor("raw_5.txt")
        raw6 <- extractor("raw_6.txt")
        raw7 <- extractor("raw_7.txt")
        
        stopifnot(identical(names(raw2$mean), names(raw4$mean)))
        stopifnot(identical(names(raw2$mean), names(raw5$mean)))
        stopifnot(identical(names(raw2$mean), names(raw6$mean)))
        stopifnot(identical(names(raw2$mean), names(raw7$mean)))
        all.means <- rbind(raw2$mean, raw4$mean, raw5$mean, raw6$mean, raw7$mean)
        all.se <- rbind(raw2$se, raw4$se, raw5$se, raw6$se, raw7$se)
        pch <- c(22, 22, 24, 24, 23)
        inner.color <- c("black", "white", "black", "white", "black")

        modes <- c("Half effect", "Variable cells", "Variable sizes", "Zero inflation", "More plates")
        out.pic <- file.path(out.dir, "supp_raw.pdf")
        devfun <- pdf
        width <- 11
        height <- 7
    } else if (mode==3L) {
        # Main figure for summed results.
        raw1 <- extractor("raw_1.txt")
        sum1 <- extractor("sum_1.txt")
        sum3 <- extractor("sum_3.txt")

        raw1$mean <- raw1$mean[match(names(sum1$mean), names(raw1$mean))]
        raw1$se <- raw1$se[match(names(sum1$se), names(raw1$se))]
        stopifnot(identical(names(raw1$mean), names(sum1$mean)))
        stopifnot(identical(names(raw1$mean), names(sum3$mean)))

        all.means <- rbind(sum1$mean, raw1$mean, sum3$mean)
        all.se <- rbind(sum1$se, raw1$se, sum3$se)
        pch <- c(21, 21, 21)
        inner.color <- c("black", "grey80", "white")

        modes <- c("Plate effect", "Unsummed", "No plate effect")
        out.pic <- file.path(out.dir, "main_sum.eps")
        width <- 5
        height <- 7
    } else if (mode==4L) {
        # Supplementaries for summed results.
        sum2 <- extractor("sum_2.txt")
        sum4 <- extractor("sum_4.txt")
        sum5 <- extractor("sum_5.txt")
        sum6 <- extractor("sum_6.txt")
        sum7 <- extractor("sum_7.txt")
 
        stopifnot(identical(names(sum2$mean), names(sum4$mean)))
        stopifnot(identical(names(sum2$mean), names(sum5$mean)))
        stopifnot(identical(names(sum2$mean), names(sum6$mean)))
        stopifnot(identical(names(sum2$mean), names(sum7$mean)))
        all.means <- rbind(sum2$mean, sum4$mean, sum5$mean, sum6$mean, sum7$mean)
        all.se <- rbind(sum2$se, sum4$se, sum5$se, sum6$se, sum7$se)
        pch <- c(22, 22, 24, 24, 23)
        inner.color <- c("black", "white", "black", "white", "black")

        modes <- c("Half effect", "Variable cells", "Variable sizes", "Zero inflation", "More plates")
        out.pic <- file.path(out.dir, "supp_sum.pdf")
        devfun <- pdf
        width <- 7
        height <- 7
    }

    coords <- matrix(1, nrow=nrow(all.means), ncol=ncol(all.means))
    coords[1,] <- coords[1,] + 1
    coords[] <- cumsum(coords)

    devfun(out.pic, width=width, height=height)
    layout(rbind(c(1,2)), width=c(width-1.5, 1.5))
    par(mar=c(6.2, 4.1, 2.1, 1.1))

    xbounds <- c(min(coords)-1, max(coords)+1)
    ybounds <- c(-3, 0)
    plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="Observed type I error rate", ylim=ybounds, xlim=xbounds, cex.lab=1.4, cex.main=1.4, bty="n", xaxs="i")

    # Adding the background bars.
    bar.cols <- rep(c("white", "lightblue"), length.out=ncol(coords))
    for (i in seq_len(ncol(coords))) { 
        rect(coords[1,i]-1, ybounds[1]-1, xbounds[2]+1, ybounds[2]+1, border=NA, col=bar.cols[i])
    }
    box()

    # Adding the axes.
    xticks <- colMeans(coords)
    names(xticks) <- renamed[names(raw1$mean)]
    axis(1, xticks, names(xticks), las=2, cex.axis=1.2)

    yticks <- seq(ybounds[1], ybounds[2], by=1)
    axis(2, yticks, labels=10^yticks, cex.axis=1.2)

    # Adding the error bars.
    all.combo <- all.means + all.se
    segments(coords, all.means, coords, all.combo)
    segments(coords-0.2, all.combo, coords+0.2, all.combo, lwd=2)

    # Adding the points.
    points(coords, all.means, pch=pch, col="black", bg=inner.color)

    # Adding the threshold.
    abline(h=log10(as.numeric(0.01)), col="red", lwd=2, lty=2)

    # Adding a legend.
    par(mar=c(0.2, 0.1, 5.1, 0.1))
    plot(0, 0, type="n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
    legend("topleft", legend=modes, pch=pch, col="black", pt.bg=inner.color)

    dev.off()
}

#########################################################################################
# End.

