setwd('~/Documents/transcriptomics/')

## METADATA
library(tidyverse)

### STEP 1: create file with just the "header" of the series matrix file
df <- read.table('data/GSE8671_series_matrix.head.txt', header=FALSE, sep='\t')

### STEP 2: transpose data
dft <- t(df)

### STEP 3: format data (optional)
dft[1,] <- gsub('!', '', dft[1,])

colnames(dft) <- dft[1,]

dft <- dft[-1,]
dft <- as.data.frame(dft)

### STEP 4: select covariates
meta <- dft[,c(2,10,11)]

### STEP 5: format data (optional)
colnames(meta) <- c('Sample', 'Tissue', 'Location')
meta <- meta %>%
  mutate(Tissue = ifelse(grepl('adenoma', Tissue), 'ADENOMA', 
                         ifelse(grepl('normal', Tissue), 'NORMAL', 'UNK'))) %>%
  mutate(Location = gsub('LOCATION: ', '', toupper(Location)))

table(meta$Tissue)
table(meta$Tissue, meta$Location)

write.table(meta, 'data/GSE8671_metadata.csv', quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
