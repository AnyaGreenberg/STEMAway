## METADATA
library(tidyverse)

### STEP 1: create file with just the "header" of the series matrix file
df <- read.table('data/GSE19804_series_matrix.head.txt', header=FALSE)

### STEP 2: transpose data
dft <- t(df)

### STEP 3: format data (optional)
dft[1,] <- gsub('!', '', dft[1,])

colnames(dft) <- dft[1,]

dft <- dft[-1,]
dft <- as.data.frame(dft)

### STEP 4: select covariates
meta <- dft[,c(2,10,12,13)]

### STEP 5: format data (optional)
colnames(meta) <- c('Sample', 'Tissue', 'Age', 'Stage')
meta <- meta %>%
  mutate(Tissue = ifelse(grepl('cancer', Tissue), 'CANCER', 
                         ifelse(grepl('normal', Tissue), 'NORMAL', 'UNK'))) %>%
  mutate(Age = gsub('age: ', '', Age)) %>%
  mutate(Stage = gsub('stage: ', '', Stage)) %>%
  mutate(Stage = gsub('n/a', 'NORMAL', Stage))

table(meta$Tissue)
table(meta$Tissue, meta$Stage)

write.table(meta, 'data/metadata.csv', quote=FALSE, row.names=FALSE, col.names=TRUE)

## ARRAY DATA
library(affy)

gse <- ReadAffy(compress=TRUE, celfile.path='data/GSE19804_RAW/')
