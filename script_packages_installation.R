## set working directory-----
getwd()
setwd("E:/Bioinformatics/How to install packages")

## Load libraries-----
!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
library(dplyr)

# read in the data ---------
dat <- read.csv(file = "E:\\Bioinformatics\\GSE183947_fpkm.csv")

dim(dat)

# get metadata --------
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

gse

metadata <- pData(phenoData(gse[[1]]))

head(metadata)

# select, mutate, rename ------------
metadata.subset <- select(metadata, c(1,10,11,17))

or

# we can use piping operator----
metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))
  head()
  
# we have to make a single dataset to avaoid masssy dataset----
metadata.modified <-metadata %>%
    select(1,10,11,17) %>%
    rename(tissue = characteristics_ch1) %>%
    rename(metastasis = characteristics_ch1.1) %>%
    mutate(tissue = gsub("tissue: ", "", tissue)) %>%
    mutate(metastasis = gsub("metastasis: ", "", metastasis))

# looking at gene expression data ---------
  head(dat)
  
# reshaping data - from wide to long-------- 
  dat %>%
  rename(gene = X) %>%
    gather(key = 'samples', value = 'FPKM', -gene)
    head()

# We have to save above long dataset into a single dataset-------- 
dat.long <- dat %>%
      rename(gene = X) %>%
      gather(key = 'samples', value = 'FPKM', -gene)

# join dataframes = dat.long + metadata.modified---
dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))
head()

# join dataframes = dat.long + metadata.modified save this like---
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# explore data ------
# filter, group_by, summarize and arrange 
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue.x) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
arrange(mean_FPKM)
