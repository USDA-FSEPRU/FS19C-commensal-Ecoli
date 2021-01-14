#####################################################################################################
#FS19C QC - MDS of fastANI and Mash output
#Kathy Mou

#Purpose: Convert ANI pairwise genome-genome similarity calculations from Mash and fastANI to a distance matrix and visualize via MDS.

#Load library packages
library(ggplot2)
library(tidyverse)
library(ggrepel) #for geom_text_repel

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

#Calculate distance
fast_dist <- fast_tab %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id') %>%
  as.dist() #distance matrix computation that computes distances between rows of a data matrix

#Generate MDS
fast_mds <- cmdscale(fast_dist) %>% as.data.frame() %>%
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



### Mash ###
#Import mash file, add ani_dist column, remove ANI / pvalue / matching_hashes columns
mash_tab <- read_tsv('./Files/distances.tab', col_names = c("reference_id", "query_id", "ANI", "pvalue", "matching_hashes")) %>%
  mutate(ani_dist=1-ANI) %>%
  select(-ANI) %>% 
  select(-pvalue) %>% 
  select(-matching_hashes)

#Calculate distance
mash_dist <- mash_tab %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id') %>%
  as.dist() #distance matrix computation that computes distances between rows of a data matrix

#Generate MDS
mash_mds <- cmdscale(mash_dist) %>% as.data.frame() %>%
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