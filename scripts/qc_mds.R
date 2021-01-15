#####################################################################################################
#FS19C QC - MDS of fastANI and Mash output
#Kathy Mou

#Purpose: Convert ANI pairwise genome-genome similarity calculations from Mash and fastANI to a distance matrix and visualize via MDS, heatmap.

#Load library packages
library(ggplot2)
library(tidyverse)
library(ggrepel) #for geom_text_repel()
install.packages("pheatmap")
library(pheatmap)
install.packages("gtools")
library(gtools) #for mixedorder() to sort reference_id numerically, which has mixed numeric and characters

sessionInfo()
#R version 4.0.2 (2020-06-22)


#####################################################################################################

#In mash and fastani files, removed extra path names (/home/..., _pol.fasta and .fasta) from elements in reference_id and query_id in excel before importing distances.tab to R so that distance calculation would work
#Otherwise, it produces errors saying duplicate 'row.names' are not allowed, non-unique values when setting 'row.names'


### FastANI ###
#Import fastani file, remove orthologous_matches, total_seq_fragments columns
fast_tab <- read_tsv('./Files/fs19cfastanioutput2.out.tab', col_names = c("reference_id", "query_id", "ani_dist", "orthologous_matches", "total_seq_fragments")) %>%
  select(-orthologous_matches) %>% 
  select(-total_seq_fragments)

#Make distance matrix
fast_dist <- fast_tab %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id') 

#Heatmap of fastani distance matrix as another way of looking at clustering
#Adapted code from: https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
fast_heatmap <- pheatmap(fast_dist, scale = "row", main = "fastANI distance matrix heatmap")
fast_heatmap
ggsave("FS19C_fastani_heatmap.tiff", plot=fast_heatmap, width = 13, height = 14, dpi = 500, units =c("in"))

#Pairwise correlation between samples (columns)
cols.cor <- cor(fast_dist, use = "pairwise.complete.obs", method='pearson')
#pairwise correlation between samples (rows)
rows.cor <- cor(t(fast_dist), use = 'pairwise.complete.obs', method='pearson')

#Plot heatmap of pairwise correlations
fast_corr_heatmap <- pheatmap(fast_dist, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor),
                              clustering_distance_rows = as.dist(1 - rows.cor),
                              main = "fastANI pairwise correlation heatmap")
fast_corr_heatmap
ggsave("FS19C_fastani_correlation_heatmap.tiff", plot=fast_corr_heatmap, width = 11, height = 14, dpi = 500, units =c("in"))

#Calculate distance for MDS
fast_dist2 <- as.dist(fast_dist) #distance matrix computation that computes distances between rows of a data matrix
head(fast_dist2)

#Generate MDS
fast_mds <- cmdscale(fast_dist2) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')
#cmdscale = classic MDS of a data matrix

#Plot MDS
plot_fast_mds <-
  fast_mds %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("fastANI MDS")
plot_fast_mds
ggsave("FS19C_fastaniMDS.tiff", plot=plot_fast_mds, width = 9, height = 8, dpi = 500, units =c("in"))



### Mash ###
#Import mash file, add ani_dist column, remove ANI / pvalue / matching_hashes columns
mash_tab <- read_tsv('./Files/distances.tab', col_names = c("reference_id", "query_id", "ANI", "pvalue", "matching_hashes")) %>%
  mutate(ani_dist=1-ANI) %>%
  select(-ANI) %>% 
  select(-pvalue) %>% 
  select(-matching_hashes)
mash_tab <- mash_tab[mixedorder(as.character(mash_tab$reference_id)),] #I reordered reference_id to make samples list in numerical order versus ordering like 1, 10-19, 2, 20-29, etc.
#I hoped this could fix FS19C_mash_correlation_heatmap.tiff so that it could list samples on the right of heatmap in the same order as FS19C_fastani_correlation_heatmap.tiff. Did not work :(

#Make distance matrix
mash_dist <- mash_tab %>%
  arrange(reference_id) %>% 
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id')

#Heatmap of mash distance matrix as another way of looking at clustering
mash_heatmap <- pheatmap(mash_dist, scale = "row", main = "mash distance matrix heatmap")
mash_heatmap
ggsave("FS19C_mash_heatmap.tiff", plot=mash_heatmap, width = 13, height = 14, dpi = 500, units =c("in"))

#Pairwise correlation between samples (columns)
cols.cor.mash <- cor(mash_dist, use = "pairwise.complete.obs", method='pearson')
#pairwise correlation between samples (rows)
rows.cor.mash <- cor(t(mash_dist), use = 'pairwise.complete.obs', method='pearson')

#Plot heatmap of pairwise correlations
mash_corr_heatmap <- pheatmap(mash_dist, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor),
                              clustering_distance_rows = as.dist(1 - rows.cor),
                              main = "mash pairwise correlation heatmap")
mash_corr_heatmap
ggsave("FS19C_mash_correlation_heatmap.tiff", plot=mash_corr_heatmap, width = 11, height = 14, dpi = 500, units =c("in"))

  
#Calculate distance for MDS
mash_dist2 <- as.dist(mash_dist) #distance matrix computation that computes distances between rows of a data matrix
head(mash_dist2)

#Generate MDS
mash_mds <- cmdscale(mash_dist2) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')
#cmdscale = classic MDS of a data matrix

#Plot MDS
plot_mash_mds <-
  mash_mds %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
plot_mash_mds
ggsave("FS19C_mashMDS.tiff", plot=plot_mash_mds, width = 9, height = 8, dpi = 500, units =c("in")) #this was the original mashMDS (reference_id not sorted, it listed entries as 1, 10-19, 2, 2-29, etc.)
ggsave("FS19C_mashMDS2.tiff", plot=plot_mash_mds, width = 9, height = 8, dpi = 500, units =c("in")) #this is a different mashMDS after I used mixedorder on reference_id





#Jules sample code
ani_tab <- read_tsv('./Files/ANIm_percentage_identity.tab') %>%
  pivot_longer(cols = -X1, names_to='genome', values_to='ANI') %>%
  mutate(ani_dist=1-ANI) %>%
  select(-ANI)

ani_dist <- ani_tab %>%
  pivot_wider(names_from = genome, values_from=ani_dist) %>%
  column_to_rownames(var='X1') %>%
  as.dist()

ani_mds <- cmdscale(ani_dist) %>% as.data.frame() %>%
  rownames_to_column(var='genome')

p_mds <-
  ani_mds %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=genome), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2')

p_mds