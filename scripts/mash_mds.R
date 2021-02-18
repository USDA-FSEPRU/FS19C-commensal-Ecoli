#!/usr/bin/env Rscript

#####################################################################################################
#FS19C QC - MDS of Mash output
#Kathy Mou

#Purpose: Convert mash output to a distance matrix and visualize distance matrices with MDS and heatmaps.

#Load library packages
library(ggplot2)
library(tidyverse)
library(ggrepel) #for geom_text_repel()
library(gtools) #for mixedorder() to sort reference_id numerically, which has mixed numeric and characters

sessionInfo()
#R version 4.0.2 (2020-06-22)

#####################################################################################################

#Mash already did the distance calculation for you (hence the small decimal numbers) so you don't need to do distance calculation (1-ANI).

############################################ Mash #########################################################
######## Import mash file, add ani_dist column, remove pvalue, matching_hashes columns ########
mash_tab <- read_tsv('./Files/distances_thirdrun.tab', col_names = c("reference_id", "query_id", "ani_dist", "pvalue", "matching_hashes")) %>%
  select(-pvalue) %>%
  select(-matching_hashes)

######## Make distance matrix ########
#mash_dist with mash_tab
mash_dist <- mash_tab %>%
  arrange(reference_id) %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id')

######## Calculate distance for MDS ########
#using mash_dist
mash_distA <- as.dist(mash_dist) #distance matrix computation that computes distances between rows of a data matrix
head(mash_distA)

######## Generate MDS ########
#using mash_distA
mash_mdsA <- cmdscale(mash_distA) %>% as.data.frame() %>%
  rownames_to_column(var='reference_id')
#cmdscale = classic MDS of a data matrix

######## Plot MDS ########
#using mash_mdsA
mash_mdsA$label <- mash_mdsA$reference_id
mash_mdsA$label[1:95] <- "" 
plot_mash_mdsA <-
  mash_mdsA %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=label), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
plot_mash_mdsA
ggsave("FS19C_mashMDS_thirdrun_onlyrefgenomes.tiff", plot=plot_mash_mdsA, width = 9, height = 8, dpi = 500, units =c("in"))

mash_mdsA$label <- mash_mdsA$reference_id
mash_mdsA$label[1:95] <- "" 
plot_mash_mdsA <-
  mash_mdsA %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
plot_mash_mdsA
ggsave("FS19C_mashMDS_thirdrun_all.tiff", plot=plot_mash_mdsA, width = 9, height = 8, dpi = 500, units =c("in"))
