#####################################################################################################
#FS19C QC - MDS
#Kathy Mou

#Purpose: Convert ANI pairwise genome-genome similarity calculations from Mash and fastANI to a distance matrix and visualize via MDS.

#Load library packages
library(ggplot2)
library(tidyverse)

sessionInfo()
#R version 4.0.2 (2020-06-22)


#####################################################################################################

#Import file and reorganize column
mash <- read_tsv('./Files/distances.tab') %>%
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