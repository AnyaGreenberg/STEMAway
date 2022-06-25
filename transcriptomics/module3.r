setwd('~/Documents/transcriptomics/')

library(tidyverse)
library(affy)
library(arrayQualityMetrics)
library(affyPLM)
library(pheatmap)

gse <- ReadAffy(compress=TRUE, celfile.path='data/GSE8671_RAW')


## ARRAY QUALITY METRICS
arrayQualityMetrics(gse, 'quality_control/array_quality_metrics/', force=TRUE, do.logtransform=TRUE)


## AFFYPLM
pset <- fitPLM(gse, background=TRUE, normalize=TRUE)
par(mai=c(2,1,0.5,0.5))
RLE(pset, main='RLE', ylab='Probe Intensities', las=2)

par(mai=c(2,1,0.5,0.5))
NUSE(pset, main='NUSE', ylab='Probe Intensities', las=2)


## CORRECTION/NORMALIZATION
norm <- rma(gse)
rma <- exprs(norm)
colnames(rma) <- gsub('.CEL.gz' ,'', colnames(rma))
write.table(rma, 'data/GSE8671_rma.csv', col.names=TRUE, row.names=TRUE, quote=FALSE, sep=',')


## BOXPLOT
boxplot(gse, main='Raw', ylab='Probe Intensities', las=2, col=c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
boxplot(norm, main='RMA', ylab='Probe Intensities', las=2, col=c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))


## PCA
meta <- read.csv('data/GSE8671_metadata.csv', sep='\t')
rownames(meta) <- meta$Sample
group <- as.factor(meta$Tissue)

pca_raw <- prcomp(exprs(gse), scale=F, center=F)
pca <- as.data.frame(pca_raw$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c('hotpink1', 'navy'))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle('PCA - Raw')+ theme_bw()

pca_norm <- prcomp(norm, scale=F, center=F)
pca <- as.data.frame(pca_norm$rotation)
ggplot(pca, aes(x=PC1, y=PC2, color=group))+
  geom_point()+ stat_ellipse()+
  scale_colour_manual(values=c('hotpink1', 'navy'))+
  geom_text(aes(label=row.names(meta)),hjust=0, vjust=0)+
  ggtitle('PCA - RMA')+ theme_bw()


## HEATMAP
group <- data.frame(Tissue=meta$Tissue)
row.names(group) <- paste(meta$Sample, '.CEL.gz', sep='')

dismat_raw <- 1-cor(exprs(gse))
pheatmap(dismat_raw, annotation_col=group, annotation_row=group, main='Heatmap - Raw', labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(ADENOMA='hotpink1', NORMAL='navy')))

dismat_norm <- 1-cor(rma)
pheatmap(dismat_norm, annotation_col=group, annotation_row=group, main='Heatmap - RMA', labels_col=meta$Tissue, 
         annotation_colors=list(Tissue=c(ADENOMA='hotpink1', NORMAL='navy')))
