# Before you running this script, the only one thing that you should prepare is the file which
# contains the expression data. Better the TPM form.
# Here, the expression file is:
# 'COADREAD_trans_GeneExp_TPM_Normal.csv': which was obtained from Counts form (with gencode file) 
# or from FPKM (without gencode file) by using the python script [TCGADownloader.script.calc*.*]
# Note: gencode file - 'gencode.v22.annotation.used4FPKM.csv', which was extracted from the GTF file, 
# but contained less but nessesary info. See more detail in [TCGADownloader.script.gencode.getFeatureLen.awk]

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/TCGADownloader/script/')

library(RUVSeq)

normal <- read.csv(file = './data/COADREAD_trans_GeneExp_TPM_Normal.csv', header = TRUE, row.names = 1)
# 每个基因的表达值大于 5 的样本个数大于 25
filter <-apply(normal, 1, function(x) length(x[x>1])>=25)
filtered <- normal[filter,]
# filtered <- normal
genes <-rownames(filtered)
rawcoln <- colnames(filtered)

# By checking the RLE and PCA behaviors of raw data, discard some samples
# filtered$TCGA.AA.3489.11A.01R.1839.07 <- NULL
# filtered$TCGA.AZ.6603.11A.02R.1839.07 <- NULL
# filtered$TCGA.AF.2692.11A.01R.A32Z.07 <- NULL
# filtered$TCGA.AG.3725.11A.01R.1736.07 <- NULL
# filtered$TCGA.AG.3742.11A.01R.1660.07 <- NULL
# filtered$TCGA.AA.3496.11A.01R.1839.07 <- NULL
# filtered$TCGA.AA.3713.11A.01R.1723.07 <- NULL
# filtered$TCGA.AZ.6605.11A.01R.1839.07 <- NULL
# filtered$TCGA.AG.3732.11A.01R.1660.07 <- NULL

x <- c()
coln <- c()
for (item in colnames(filtered)) {
  tmp <- strsplit(x = item, split = '[.]')
  x <- c(x, tmp[[1]][5])
  coln <- c(coln, paste(tmp[[1]][2], tmp[[1]][3], sep = '_'))
}
colnames(filtered) <- coln
x <- factor(x)

# By checking the RLE and PCA behaviors of raw data, discard some samples
selectedSample <- c('AZ_6600', 'A6_2685', 'AA_3534', 'AA_3660', 'A6_2683', 'AA_3522', 'AA_3527',
                    'AA_3663', 'AA_3525', 'A6_2671', 'A6_2682', 'AZ_6598', 'AA_3520', 'AF_2691',
                    'AF_2689', 'A6_5662')
x <- x[which(colnames(filtered) %in% selectedSample)]
filtered <- filtered[,selectedSample]

colors <-RColorBrewer::brewer.pal(3, "Set2")
# RLE and PCA of raw data
set <-newSeqExpressionSet(as.matrix(filtered), phenoData =data.frame(x, row.names=colnames(filtered)))
png(filename = 'RAW_RLE_TPM_rm.png', width = 1000, height = 500)
plotRLE(set, outline=FALSE, col=colors[x], cex.lab=0.2, las=2, srt=45)
dev.off()
png(filename = 'RAW_PCA_TPM_rm.png', width = 500, height = 500)
plotPCA(set, k=2, col=colors[x], cex=1)
dev.off()

write.csv(x = filtered, file = 'COADREAD_trans_GeneExp_TPM_Normal_selected.csv')

# Skip the following steps
# RLE and PCA after UQ-normalization
UQset <- betweenLaneNormalization(set, which="upper", round=FALSE)
png(filename = 'UQ_RLE_TPM_rm.png', width = 1000, height = 500)
plotRLE(UQset, outline=FALSE, col=colors[x], cex.lab=0.5, las=2, srt=45)
dev.off()
png(filename = 'UQ_PCA_TPM_rm.png', width = 500, height = 500)
plotPCA(UQset, k=2, col=colors[x], cex=1)
dev.off()

norm_data <- UQset@assayData$normalizedCounts
write.csv(x = norm_data, file = 'COADREAD_n_norm.csv')

# RUVs (skiped)
differences <-makeGroups(x)
for (j in 1:5) {
  RUVset <-RUVs(set, genes, k=j, differences)
  par(mfrow = c(5, 2), pin=c(0.5,0.5), mar=c(4.1, 3.9, 3.2, 1.1))
  plotRLE(RUVset, outline=FALSE, col=colors[x])
  plotPCA(RUVset, k=3, col=colors[x], cex=1.2)
}
## k = 3 is better (maybe)
