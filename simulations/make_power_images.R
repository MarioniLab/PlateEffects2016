cur.dir <- "ESpresso"
out.dir <- file.path(cur.dir, "results_power")

#########################################################################################
# First, making the ROC plots.

colors <- c(QLedgeR="red", DESeq2="black", voom="grey")
for (mode in 1:4) {
    cur.roc <- readRDS(file.path(out.dir, paste0("roc_", mode, ".rds")))
    raw.roc <- sum.roc <- list()
    for (method in names(cur.roc)) {
        raw.roc[[method]] <- colMeans(do.call(rbind, cur.roc[[method]]$raw))
        sum.roc[[method]] <- colMeans(do.call(rbind, cur.roc[[method]]$sum))
    }
    
    for (zoom in c(FALSE, TRUE)) {
        if (!zoom) {
            xup <- 1
            yup <- 1
            fname <- file.path(out.dir, paste0("ROC_", mode, ".eps"))
        } else {
            if (mode!=1) { break }
            xup <- 0.1
            yup <- 0.7
            fname <- file.path(out.dir, paste0("ROC_zoom_", mode, ".eps"))
        }

        setEPS()
        postscript(fname)
        par(mar=c(6.1, 5.1, 2.1, 1.1))
        plot(0,0, type="n", xlab="False positive rate", ylab="True positive rate", cex.axis=1.2, cex.lab=1.4, xlim=c(0, xup), ylim=c(0, yup))
        for (method in names(raw.roc)) { 
            current <- raw.roc[[method]]
            lines(current, seq_along(current)/length(current), col=colors[method], lwd=2)
            current <- sum.roc[[method]]
            lines(current, seq_along(current)/length(current), col=colors[method], lwd=2, lty=2)
        }
        
        if (mode <= 2L && !zoom) {
            # Adding a legend.
            all.names <- c(names(raw.roc), names(sum.roc))
            renames <- paste0(all.names, " (", rep(c("single", "sum"), c(length(raw.roc), length(sum.roc))), ")")
            legend("bottomright", lwd=2, lty=rep(1:2, each=length(raw.roc)), legend=renames, col=colors[all.names], cex=1.2)
        }
        dev.off()
    }
}

#########################################################################################
# Now making the FDR plot (basically ripped out of failsim image maker).

extractor <- function(fname) {
    out <- read.table(file.path(out.dir, fname), sep="\t", stringsAsFactors=FALSE)
    out <- split(log10(out[,3]), out[,1])
    all.means <- sapply(out, FUN=mean)
    all.se <- sqrt(sapply(out, FUN=var)/lengths(out))
    return(list(mean=all.means, se=all.se))
}

raw1 <- extractor("raw_1.txt")
raw2 <- extractor("raw_2.txt")
raw3 <- extractor("raw_3.txt")
raw4 <- extractor("raw_4.txt")
sum1 <- extractor("sum_1.txt")
sum2 <- extractor("sum_2.txt")
sum3 <- extractor("sum_3.txt")
sum4 <- extractor("sum_4.txt")
        
stopifnot(identical(names(raw1$mean), names(raw2$mean)))
stopifnot(identical(names(raw1$mean), names(raw3$mean)))
stopifnot(identical(names(raw1$mean), names(raw4$mean)))
stopifnot(identical(names(sum1$mean), names(sum2$mean)))
stopifnot(identical(names(sum1$mean), names(sum3$mean)))
stopifnot(identical(names(sum1$mean), names(sum4$mean)))

matchnames <- c(DESeq2="DESeq2", voom="voom", QLedgeR="QL edgeR")
renamed <- c(paste0(matchnames[names(raw1$mean)], "\n(single)"), paste0(matchnames[names(sum1$mean)], "\n(sum)"))
all.means <- cbind(
                   rbind(raw1$mean, raw2$mean, raw3$mean, raw4$mean),
                   rbind(sum1$mean, sum2$mean, sum3$mean, sum4$mean)
                   )
all.se <- cbind(
                rbind(raw1$se, raw2$se, raw3$se, raw4$se), 
                rbind(sum1$se, sum2$se, sum3$se, sum4$se)
                )

# Setting up the plot parameters.

pch <- c(23, 23, 25, 25)
inner.color <- c("black", "white", "black", "white")
modes <- c("Plate effect", "No plate effect", "Stronger DE", "More DE")

coords <- matrix(1, nrow=nrow(all.means), ncol=ncol(all.means))
coords[1,] <- coords[1,] + 1
coords[] <- cumsum(coords)

setEPS()
width <- 9
postscript(file.path(out.dir, "FDR.eps"), width=width, height=7)
layout(rbind(c(1,2)), width=c(width-1.5, 1.5))
par(mar=c(6.2, 4.1, 2.1, 1.1))

xbounds <- c(min(coords)-1, max(coords)+1)
ybounds <- c(-2, 0)
plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="Observed FDR", ylim=ybounds, xlim=xbounds, cex.lab=1.4, cex.main=1.4, bty="n", xaxs="i")

# Adding the background bars.
bar.cols <- rep(c("white", "lightblue"), length.out=ncol(coords))
for (i in seq_len(ncol(coords))) { 
    rect(coords[1,i]-1, ybounds[1]-1, xbounds[2]+1, ybounds[2]+1, border=NA, col=bar.cols[i])
}
box()

# Adding the axes.
xticks <- colMeans(coords)
axis(1, xticks, renamed, las=2, cex.axis=1.2)

yticks <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
axis(2, log10(yticks), labels=yticks, cex.axis=1.2)

# Adding the error bars.
all.combo <- all.means + all.se
segments(coords, all.means, coords, all.combo)
segments(coords-0.2, all.combo, coords+0.2, all.combo, lwd=2)

# Adding the points.
points(coords, all.means, pch=pch, col="black", bg=inner.color)

# Adding the threshold.
abline(h=log10(as.numeric(0.05)), col="red", lwd=2, lty=2)

# Adding a legend.
par(mar=c(0.2, 0.1, 5.1, 0.1))
plot(0, 0, type="n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
legend("topleft", legend=modes, pch=pch, col="black", pt.bg=inner.color)

dev.off()

#########################################################################################
# End.

