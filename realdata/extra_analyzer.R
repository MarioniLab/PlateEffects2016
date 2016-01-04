# Quick DE analysis script for GSE59114 on single-cell RNA-seq in aging mice.
# Note there's a couple of assumptions here; they only report TPMs, so we'll
# assume library sizes of 1 million (this seems close to their actual library 
# sizes on GEO, so we'll go with it). Also, they don't indicate which plate
# each cell came from; we'll assume those numbered 96 and below come from one
# plate, and those numbered 97 and above come from another plate.

require(edgeR)
ltpms <- read.csv("GSE59114/GSE59114_C57BL6_GEO_all.csv", header=TRUE, skip=1)
ltpms <- ltpms[,-(1:2)]
raw.counts <- 2^ltpms -  1
keep <- rowMeans(raw.counts) > 1
raw.counts <- raw.counts[keep,]

# As single cell analysis.
groupings <- factor(sub("_[0-9]+$", "", colnames(raw.counts)))
design <- model.matrix(~0 + groupings)
colnames(design) <- levels(groupings)

y <- DGEList(raw.counts)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

res1 <- glmQLFTest(fit, contrast=makeContrasts(old_MPP - young_MPP, levels=design))
summary(decideTestsDGE(res1))
res2 <- glmQLFTest(fit, contrast=makeContrasts(young_ST_HSC - young_LT_HSC, levels=design))
summary(decideTestsDGE(res2))
res3 <- glmQLFTest(fit, contrast=makeContrasts(old_ST_HSC - old_LT_HSC, levels=design))
summary(decideTestsDGE(res3))

# After summation
platedness <- as.integer(sub(".*_([0-9]+)$", "\\1", colnames(raw.counts))) <= 96
same.plate <- paste0(groupings, "_", platedness)
sum.counts <- sumTechReps(raw.counts, same.plate)
sum.groups <- factor(sub("_(TRUE|FALSE)$", "", colnames(sum.counts)))
sum.design <- model.matrix(~0 + sum.groups)
colnames(sum.design) <- levels(sum.groups)

sum.y <- DGEList(sum.counts)
sum.y <- calcNormFactors(sum.y)
sum.y <- estimateDisp(sum.y, sum.design)
sum.fit <- glmQLFit(sum.y, sum.design, robust=TRUE)

sum.res1 <- glmQLFTest(sum.fit, contrast=makeContrasts(old_MPP - young_MPP, levels=sum.design))
summary(decideTestsDGE(sum.res1))
sum.res2 <- glmQLFTest(sum.fit, contrast=makeContrasts(young_ST_HSC - young_LT_HSC, levels=sum.design))
summary(decideTestsDGE(sum.res2))
sum.res3 <- glmQLFTest(sum.fit, contrast=makeContrasts(old_ST_HSC - old_LT_HSC, levels=sum.design))
summary(decideTestsDGE(sum.res3))

# Reporting session information.

sessionInfo()
