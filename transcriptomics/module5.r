setwd('~/Documents/transcriptomics/')

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(ggnewscale) #for cnetplot

# ENRICHMENT (UP-REGULATED GENES)
## GENE LIST
### vector of logFC values, named with gene symbol, sorted in decreasing order - filter for most differentially expressed
tT <- read.table('data/dge.csv', sep=',', header=TRUE) %>%
  filter(logFC > 1.5)

deg_logfc <- tT$logFC
names(deg_logfc) <- tT$ID
deg_logfc <- sort(deg_logfc, decreasing=TRUE)

### need ENTREZID for some functions
deg_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=names(deg_logfc), columns=c('ENTREZID'), keytype='SYMBOL')
any(is.na(deg_entrez$ENTREZID))

## GENE ONTOLOGY
go_cc <- enrichGO(deg_entrez$ENTREZID, org.Hs.eg.db, ont='CC', readable=TRUE)
go_bp <- enrichGO(deg_entrez$ENTREZID, org.Hs.eg.db, ont='BP', readable=TRUE)
go_mf <- enrichGO(deg_entrez$ENTREZID, org.Hs.eg.db, ont='MF', readable=TRUE)
barplot(go_cc, title='Cellular components')
barplot(go_bp, title='Biological processes')
barplot(go_mf, title='Molecular function')

## KEGG PATHWAYS
kegg <- enrichKEGG(deg_entrez$ENTREZID)
dotplot(kegg, title='Enriched KEGG pathways')

# NETWORK
## GENE-CONCEPT NETWORK
gcnet <- setReadable(kegg, org.Hs.eg.db, keyType='ENTREZID') #map ENTREZID to gene symbol
cnetplot(gcnet, foldChange=deg_logfc, categorySize='p.adjust', colorEdge=TRUE)

## TRANSCRIPTION FACTOR NETWORK
msig <- msigdbr(species='Homo sapiens', category='C3')
c3 <- msig %>% select(gs_name, entrez_gene)

tf <- enricher(deg_entrez$ENTREZID, TERM2GENE=c3)
tfnet <- setReadable(tf, org.Hs.eg.db, keyType='ENTREZID')
cnetplot(tfnet, foldChange=deg_logfc, categorySize='p.adjust', colorEdge=TRUE)

## EXTERNAL TOOLS
write.table(names(deg_logfc), 'data/deg.up.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)


# OVER-REPRESENTATION
## GSEA
### GENE LIST
### vector of logFC values, named with ENTREZID, sorted in decreasing order - include all genes
tT <- read.table('data/dge.csv', sep=',', header=TRUE)
gsea_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=tT$ID, columns=c('ENTREZID'), keytype='SYMBOL')
gsea_entrez <- gsea_entrez[!duplicated(gsea_entrez$SYMBOL),]
gsea_entrez <- gsea_entrez[!duplicated(gsea_entrez$ENTREZID),]
gsea_entrez <- gsea_entrez %>% rename(ID=SYMBOL)
gsea_genes <- merge(gsea_entrez, tT, by='ID')

genelist <- gsea_genes$logFC
names(genelist) <- gsea_genes$ENTREZID
genelist <- sort(genelist, decreasing=TRUE)

## ANALYSIS
h <- msig %>% select(gs_name, entrez_gene)
gsea <- GSEA(genelist, TERM2GENE=h, eps=0)
par(mai=c(1,1,1,2))
gseaplot2(gsea, geneSetID=1:5, pvalue_table=TRUE)
