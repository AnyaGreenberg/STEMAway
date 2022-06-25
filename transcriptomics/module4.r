setwd('~/Documents/transcriptomics/')

library(tidyverse)
library(hgu133plus2.db)
library(limma)
library(EnhancedVolcano)
library(pheatmap)

rma <- read.table('data/GSE8671_rma.csv', sep=',')
rma <- rma[,-c(13,43)]

## ANNOTAIONS
symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(rma), columns=c('SYMBOL', 'ENTREZID'))
symbols <- symbols[!duplicated(symbols$PROBEID),]

all(row.names(rma) == symbols$PROBEID)

rma$SYMBOL <- symbols$SYMBOL

rma <- na.omit(rma)

rma <- rma[!duplicated(rma$SYMBOL),]
rownames(rma) <- rma$SYMBOL


## GENE FILTERING
means <- rowMeans(rma[,1:62])
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=TRUE))

filt <- rma[which(means >= perc2),1:62]
write.table(filt, 'data/GSE8671_filtered.csv', col.names=TRUE, row.names=TRUE, quote=FALSE, sep=',')


## LIMMA
meta <- read.csv('data/GSE8671_metadata.csv', sep='\t')
meta <- meta[-c(13,43),]
Tissue <- factor(meta$Tissue, levels=c('ADENOMA', 'NORMAL'), order=FALSE)
row.names(meta) <- meta$Sample
design <- model.matrix(~0+Tissue, meta)

lm <- lmFit(filt, design)
contrast <- makeContrasts(TissueADENOMA-TissueNORMAL, levels=design)
fit <- contrasts.fit(lm, contrast)
efit <- eBayes(fit)

tT <- topTable(efit, p.value=0.05, adjust.method='fdr', sort.by='logFC', genelist=row.names(filt), number=length(row.names(efit)))


## Volcano plot
EnhancedVolcano(tT, lab=tT$ID, x='logFC', y='adj.P.Val', 
                pCutoff=1e-5, FCcutoff=2,
                pointSize=1, legendLabSize=10, labSize=3, 
                title='GSE8671 DEG', subtitle='RMA Normalization')


## HEATMAP
up <- tT %>%
  filter(adj.P.Val <= 1e-5 & logFC >= 1)
up <- up[order(up$logFC, decreasing=TRUE),]

down <- tT %>%
  filter(adj.P.Val <= 1e-5, logFC <= -1)
down <- down[order(down$logFC, decreasing=FALSE),]

hgenes <- rbind(head(up,25), head(down,25))
hdf <- filt[row.names(filt) %in% hgenes$ID,]
group <- data.frame(Location=meta$Location)
row.names(group) <- meta$Sample
pheatmap(hdf, annotation_col=group, cluster_rows=T, main="Top DEG") 
