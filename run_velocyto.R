library(velocyto.R)
library(Matrix)
ldat <- read.loom.matrices("velocyto/combined_4samples_v2.loom")
print(str(ldat))
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
dev.off()
emat <- ldat$spliced
nmat <- ldat$unspliced

cell.colors <- data.frame("Rep1_6T8PI:wgEncodeCshlLongRnaSeqGm12878CellTotalAlnRep1" = "gm12878", "Rep2_E9TS0:wgEncodeCshlLongRnaSeqGm12878CellTotalAlnRep2" = "gm12878", "K562_Rep1_withtags_TXR6K:K562_Rep1" = "k562", "K562_Rep2_withtags_OWVYP:K562_Rep2" = "k562")
names(cell.colors) <- colnames(emat)
print(colnames(emat) %in% names(cell.colors))
print(colnames(emat))
print(names(cell.colors))

# filter expression matrices based on some minimum max-cluster averages
#emat <- filter.genes.by.cluster.expression(as.matrix(emat),cell.colors,min.max.cluster.average = 5)
#nmat <- filter.genes.by.cluster.expression(as.matrix(nmat),cell.colors,min.max.cluster.average = 1)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))
fit.quantile <- 0.05;

pdf("velocity.pdf")
rvel1 <- gene.relative.velocity.estimates(as.matrix(emat),as.matrix(nmat),deltaT=1,deltaT2 = 1,kCells = 1, fit.quantile=fit.quantile)
save(rvel1, file = "rvel1.robj")
pca.velocity.plot(rvel1, nPcs=3, plot.cols=2, ,cex=1.2,pcount=0.1,pc.multipliers=c(1,-1, 1))
dev.off()
#rvel.qf <- gene.relative.velocity.estimates(as.matrix(emat), as.matrix(nmat), deltaT=1, kCells = 0, fit.quantile = fit.quantile)
