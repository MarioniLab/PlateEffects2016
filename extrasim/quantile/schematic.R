# This generates the schematic requried to explain the simulation design.

pdf("sim_schematic.pdf", width=12, height=6)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-2, 58), ylim=c(0, 25))

top.plate <- 12
bottom.plate <- 0
plate.size <- 10
shift <- 13

text(-1, c(bottom.plate+plate.size/2), expression(g[2]), pos=2, cex=1.5)
text(-1, c(top.plate+plate.size/2), expression(g[1]), pos=2, cex=1.5)

for (scenario in 1:4) { 
    left.edge <- (scenario - 1L) * shift
    rect(left.edge, bottom.plate, left.edge+plate.size, bottom.plate+plate.size)  
    rect(left.edge, top.plate, left.edge+plate.size, top.plate+plate.size)  
    text(left.edge + plate.size/2, 23, scenario, pos=3, cex=1.5)

    if (scenario==1L) {
        top.cutoff <- 100
        bottom.cutoff <- 10
        ncolored <- 0
        colval <- "white" 
    } else {
        top.cutoff <- 50
        bottom.cutoff <- 50

        if (scenario==2L) {
            ncolored <- 20
            colval <- "yellow"
        } else if (scenario==3L) {
            ncolored <- 10
            colval <- "orange"
        } else if (scenario==4L) {
            ncolored <- 5
            colval <- "red"
        }
    }

    for (plate.type in c("top", "bottom")) { 
        if (plate.type=="top") {
            plate.y <- top.plate
            cutoff <- top.cutoff
        } else {
            plate.y <- bottom.plate
            cutoff <- bottom.cutoff 
        }

        for (y in seq_len(plate.size)) {
            remaining <- cutoff - (y-1)*plate.size
            if (remaining <= 0) {  break }
            to.add <- pmin(remaining, plate.size)

            if (plate.type=="top") {
                remaining.col <- ncolored - (y-1)*plate.size
                all.cols <- rep("grey80", plate.size)
                if (remaining.col > 0) { 
                    all.cols[1:pmin(remaining.col, plate.size)] <- colval
                }
            }
            
            symbols(left.edge + seq_len(to.add)-0.5, plate.y+rep(y, to.add)-0.5, circles=rep(0.3, to.add), 
                    bg=all.cols, add=TRUE, inches=FALSE)
        }
    }
}


legend(51, top.plate+2, legend=c("1x", "5x", "10x", "20x"), fill=c("grey80", "yellow", "orange", "red"), cex=1.5) 
dev.off()
