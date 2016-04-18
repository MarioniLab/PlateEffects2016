#################################################

library(org.Mm.eg.db)

add.gene.names <- function(stuff) {
    entrez <- select(org.Mm.eg.db, keys=rownames(stuff), keytype="ENSEMBL", column=c("ENTREZID", "SYMBOL"))
    entrez <- entrez[match(rownames(stuff), entrez$ENSEMBL),]
    cbind(entrez, stuff)
}

# Pulled out of Supplementary figure S5 of the mESC paper.
of.interest <- c("Prdm1", "Bmp4", "Ccnf", "Cdh1", "Dnmt3a", "Dnmt3l", "Dppa4", "Eed", "Esrrb", "Etv5", "Ezh2", "Fbxo15",
                 "Fgf4", "Gdf3", "Il6st", "Idh1", "Jarid2", "Kdm6a", "Klf2", "Klf4", "Lifr", "Kdm1a", "Nanog", "Nodal",
                 "Nr0b1", "Nr5a2", "Pou5f1", "Otx2", "Pfkm", "Pfkp", "Prdm14", "Rarg", "Sall1", "Sall4", "Smarcd1", "Sox2",
                 "Dppa3", "Suz12", "Tdgf1", "Tead1", "Tert", "Tet1", "Tet2", "Tfap2c", "Tfcp2l1", "Utf1", "Zfp42")

contrasted <- c("3"="red", "4"="blue")

#################################################

sum3 <- read.table("ESpresso/edgeR_3_sum.tsv.gz", header=TRUE)
raw3 <- read.table("ESpresso/edgeR_3_raw.tsv.gz", header=TRUE)
sum3 <- add.gene.names(sum3)
raw3 <- add.gene.names(raw3)

ribo.sum <- match(of.interest, sum3$SYMBOL)
ribo.raw <- match(of.interest, raw3$SYMBOL)
r1 <- rank(sum3$PValue)[ribo.sum]
r2 <- rank(raw3$PValue)[ribo.raw]

x3 <- pmin(r1, r2)
y3 <- r2 - r1
mat3 <- cbind(c(sum(r1 <= 250), sum(r1 > 250)), c(sum(r2 <= 250), sum(r2 > 250)))

sum4 <- read.table("ESpresso/edgeR_4_sum.tsv.gz", header=TRUE)
raw4 <- read.table("ESpresso/edgeR_4_raw.tsv.gz", header=TRUE)
sum4 <- add.gene.names(sum4)
raw4 <- add.gene.names(raw4)

ribo.sum <- match(of.interest, sum4$SYMBOL)
ribo.raw <- match(of.interest, raw4$SYMBOL)

r1 <- rank(sum4$PValue)[ribo.sum]
r2 <- rank(raw4$PValue)[ribo.raw]
x4 <- pmin(r1, r2)
y4 <- r2 - r1
mat4 <- cbind(c(sum(r1 <= 250), sum(r1 > 250)), c(sum(r2 <= 250), sum(r2 > 250)))

fisher.test(mat3 + mat4)  # pooling to get more power.

#################################################
# Plotting.

pdf("ESpresso/real_ranks.pdf")
par(mar=c(5.1, 5.1, 2.1, 2.1))
plot(0, 0, xlim=c(0, 250), ylim=c(-500, 2200), xlab="Minimum rank", ylab="Rank difference (single - sum)", type="n", cex.lab=1.4, cex.axis=1.2)
rect(-100, -100, 500, 100, col="grey90", border=NA)
box()

points(x3, y3, pch=16, col=contrasted["3"])
text(x3, y3, of.interest, col=contrasted["3"], pos=4, cex=0.8)

points(x4, y4, pch=16, col=contrasted["4"])
chosen <- of.interest %in% c("Tfcp2l1", "Nr0b1", "Klf2", "Nanog", "Pfkp")
text(x4, y4, of.interest, col=contrasted["4"], pos=ifelse(chosen, 3, 4), cex=0.8)

legend("topleft", legend=c("2i vs. serum", "a2i vs. serum"), pch=16, col=contrasted)
dev.off()


