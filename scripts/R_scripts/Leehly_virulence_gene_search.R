#!/usr/bin/env Rscript

#####################################################################################################
#FS19C Roary - search for isolates with no LEE or hemolysin virulence genes
#Kathy Mou

#Purpose: Identify isolates that don't possess any of the LEE or hemolysin virulence genes

library(tidyverse)
library(pheatmap)


#Read file (generated from Roary) into matrix
PA <- read_tsv('gene_presence_absence.Rtab') %>%
  column_to_rownames(var = 'Gene') %>% 
  as.matrix()
colnames(PA)
rownames(PA)

#Heatmap of presence/absence of select virulence genes (doi: 10.1128/iai.68.11.6115-6126.2000)
keep_rows <- rownames(PA) %>% grep(pattern="^bfp|^cesT|^eae|^esc|^esp|^hly|^ler|^sepZ|^stx|^tir|^tagA")
PA2 <- PA[keep_rows,]
rownames(PA2)
colnames(PA2)
pheatmap(PA2, cellheight = 6) #only isolate 79 does not possess these virulence genes.

#Read file (generated from Roary) into dataframe
gene <- read.csv('gene_presence_absence.csv', sep=",")
sort(unique(gene$Non.unique.Gene.name))

#Sort out virulence genes from Non.unique.Gene.name column
nonun <- gene$Non.unique.Gene.name %>%  grep(pattern="^bfp|^cesT|^eae|^esc|^esp|^hly|^ler|^sepZ|^stx|^tir|^tagA")
nonun2 <- gene[nonun,]
virgenes <- nonun2

#Sort out virulence genes from Gene column
vir <- gene$Gene %>%  grep(pattern="^bfp|^cesT|^eae|^esc|^esp|^hly|^ler|^sepZ|^stx|^tir|^tagA")
vir2 <- gene[vir,]
virgenes <- rbind(virgenes, vir2)

write.csv(virgenes, "Virulence_genes_list.csv")
