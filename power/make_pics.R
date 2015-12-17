# This makes RDS curves using all the stuff in the thing.

for (x in c("with", "without")) {
    current <- readRDS(sprintf("%s.rds", x))

    for (zoom in c(TRUE, FALSE)) {
        setEPS()
        postscript(sprintf("%s%s.eps", x, ifelse(zoom, "_zoom", "")))
        plot(0,0, type="n", xlab="False positive rate", cex.lab=1.4, cex.axis=1.2,
             ylab="True positive rate", xlim=c(0, ifelse(zoom, 0.1, 1)), ylim=c(0, ifelse(zoom, 0.6, 1)))
        collected.names <- collected.col <- collected.lty <- list()

        for (method in names(current)) {
            curmeth <- current[[method]]
            if (method=="edgeR") { col <- "red" } 
            else if (method=="voom") { col <- "grey" }
            else { col <- "black" }

            for (modes in names(curmeth)) { 
                all.hits <- curmeth[[modes]]
                averaged <- colMeans(do.call(rbind, all.hits))

                cur.lty <- ifelse(modes=="raw", 1, 2)
                lines(averaged, seq_along(averaged)/length(averaged), lwd=2, lty=cur.lty, col=col)

                collected.names <- c(collected.names, paste(method, ifelse(modes=="raw", "", "(sum)")))
                collected.col <- c(collected.col, col)
                collected.lty <- c(collected.lty, cur.lty)
            }
        }

        if (!zoom) {
            legend(1, 0, legend=unlist(collected.names), lwd=2, col=unlist(collected.col), lty=unlist(collected.lty), xjust=1, yjust=0)
        }
        dev.off()
    }
}
