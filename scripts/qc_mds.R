#####################################################################################################
#FS19C QC - MDS of fastANI and Mash output
#Kathy Mou

#Purpose: Convert ANI pairwise genome-genome similarity calculations of fastANI and mash output to a distance matrix and visualize distance matrices with MDS and heatmaps.

#Load library packages
library(ggplot2)
library(tidyverse)
library(ggrepel) #for geom_text_repel()
library(pheatmap)
library(gtools) #for mixedorder() to sort reference_id numerically, which has mixed numeric and characters

sessionInfo()
#R version 4.0.2 (2020-06-22)

#####################################################################################################

#In mash and fastani files, I removed extra path names (/home/..., _pol.fasta and .fasta) from elements in reference_id and query_id 
#in excel before importing distances.tab to R so that distance calculation would work
#Otherwise, it produces errors saying duplicate 'row.names' are not allowed, non-unique values when setting 'row.names'
#fastANI gives you ANI similarity in percentage, mash already did the distance calculation for you (hence the small decimal numbers) 
#so you don't need to do distance calculation (1-ANI).
#You will need to convert fastANI ANI percentage values to decimal and then subtract that from 1 to get distance calculation.
#I did a little experiment to see if reordering reference_id of mash file (distance.tab) with mixedorder() would change MDS and 
#heatmap structure compared to if I didn't reorder reference_id. In conclusion, no difference, so don't need to reorder reference_id

################################### FastANI ##################################################################
######## Import fastani file, remove orthologous_matches, total_seq_fragments columns ########
fast_tab <- read_tsv('./Files/fs19cfastanioutput2.out.tab', col_names = c("reference_id", "query_id", "ani", "orthologous_matches", "total_seq_fragments")) %>%
  select(-orthologous_matches) %>% 
  select(-total_seq_fragments) %>% 
  mutate(ani_dist= 1-(ani/100)) %>% #add new ani_dist column to calculate distance
  select(-ani)

######## Make distance matrix ########
fast_dist <- fast_tab %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id') 

######## Heatmap of fastani distance matrix as another way of looking at clustering ########
#Adapted code from: https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
fast_heatmap <- pheatmap(fast_dist, scale = "row", main = "fastANI distance matrix heatmap")
fast_heatmap
ggsave("FS19C_fastani_heatmap.tiff", plot=fast_heatmap, width = 13, height = 14, dpi = 500, units =c("in"))

######## Pairwise correlation between samples ########
#by columns
cols.cor <- cor(fast_dist, use = "pairwise.complete.obs", method='pearson')
#by rows
rows.cor <- cor(t(fast_dist), use = 'pairwise.complete.obs', method='pearson')

######## Plot heatmap of pairwise correlations ########
fast_corr_heatmap <- pheatmap(fast_dist, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor),
                              clustering_distance_rows = as.dist(1 - rows.cor),
                              main = "fastANI pairwise correlation heatmap")
fast_corr_heatmap
ggsave("FS19C_fastani_correlation_heatmap.tiff", plot=fast_corr_heatmap, width = 11, height = 14, dpi = 500, units =c("in"))

######## Calculate distance for MDS ########
fast_dist2 <- as.dist(fast_dist) #distance matrix computation that computes distances between rows of a data matrix
head(fast_dist2)

######## Generate MDS ########
fast_mds <- cmdscale(fast_dist2) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')
#cmdscale = classic MDS of a data matrix

######## Plot MDS ########
plot_fast_mds <-
  fast_mds %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("fastANI MDS")
plot_fast_mds
ggsave("FS19C_fastaniMDS.tiff", plot=plot_fast_mds, width = 9, height = 8, dpi = 500, units =c("in"))

############################################ Mash #########################################################
######## Import mash file, add ani_dist column, remove pvalue, matching_hashes columns ########
mash_tab <- read_tsv('./Files/distances.tab', col_names = c("reference_id", "query_id", "ani_dist", "pvalue", "matching_hashes")) %>%
  select(-pvalue) %>% 
  select(-matching_hashes)
mash_tab2 <- mash_tab[mixedorder(as.character(mash_tab$reference_id)),] #I reordered reference_id to make samples list in numerical order versus 
#ordering like 1, 10-19, 2, 20-29, etc.
#I hoped this could fix FS19C_mash_correlation_heatmap.tiff so that it could list samples on the right of heatmap in the same order as 
#FS19C_fastani_correlation_heatmap.tiff. Did not work :(

######## Make distance matrix ########
#mash_dist with mash_tab
mash_dist <- mash_tab %>%
  arrange(reference_id) %>% 
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id')

#mash_dist2 with mash_tab2
mash_dist2 <- mash_tab2 %>%
  arrange(reference_id) %>% 
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id')

######## Heatmap of mash distance matrix as another way of looking at clustering ########
#using mash_dist
mash_heatmap <- pheatmap(mash_dist, scale = "row", main = "mash distance matrix heatmap")
mash_heatmap
ggsave("FS19C_mash_heatmap.tiff", plot=mash_heatmap, width = 13, height = 14, dpi = 500, units =c("in"))

#using mash_dist2
mash_heatmap2 <- pheatmap(mash_dist2, scale = "row", main = "mash distance matrix heatmap")
mash_heatmap2 #looks the same as mash_heatmap

######## Pairwise correlation between samples (columns) ########
#using mash_dist
cols.cor.mash <- cor(mash_dist, use = "pairwise.complete.obs", method='pearson')
#pairwise correlation between samples (rows)
rows.cor.mash <- cor(t(mash_dist), use = 'pairwise.complete.obs', method='pearson')

#using mash_dist2
cols.cor.mash2 <- cor(mash_dist2, use = "pairwise.complete.obs", method='pearson')
#pairwise correlation between samples (rows)
rows.cor.mash2 <- cor(t(mash_dist2), use = 'pairwise.complete.obs', method='pearson')

######## Plot heatmap of pairwise correlations ########
#using mash_dist
mash_corr_heatmap <- pheatmap(mash_dist, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor.mash),
                              clustering_distance_rows = as.dist(1 - rows.cor.mash),
                              main = "mash pairwise correlation heatmap")
mash_corr_heatmap
ggsave("FS19C_mash_correlation_heatmap.tiff", plot=mash_corr_heatmap, width = 11, height = 14, dpi = 500, units =c("in"))

#using mash_dist2
mash_corr_heatmap2 <- pheatmap(mash_dist2, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor.mash2),
                              clustering_distance_rows = as.dist(1 - rows.cor.mash2),
                              main = "mash pairwise correlation heatmap2")
mash_corr_heatmap2 #looks the same as mash_corr_heatmap

######## Calculate distance for MDS ########
#using mash_dist
mash_distA <- as.dist(mash_dist) #distance matrix computation that computes distances between rows of a data matrix
head(mash_distA)

#using mash_dist2
mash_distB <- as.dist(mash_dist2) #distance matrix computation that computes distances between rows of a data matrix
head(mash_distB)

######## Generate MDS ########
#using mash_distA
mash_mdsA <- cmdscale(mash_distA) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')
#cmdscale = classic MDS of a data matrix

#using mash_distB
mash_mdsB <- cmdscale(mash_distB) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')

######## Plot MDS ########
#using mash_mdsA
plot_mash_mdsA <-
  mash_mdsA %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
plot_mash_mdsA
ggsave("FS19C_mashMDS.tiff", plot=plot_mash_mdsA, width = 9, height = 8, dpi = 500, units =c("in"))

#using mash_mdsB
plot_mash_mdsB <-
  mash_mdsB %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
plot_mash_mdsB #exact same plot as plot_mash_mdsA
ggsave("FS19C_mashMDS2.tiff", plot=plot_mash_mdsB, width = 9, height = 8, dpi = 500, units =c("in")) #this is a different mashMDS after I used mixedorder on reference_id





###################################### Jules sample code ###############################################################
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