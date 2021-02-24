#!/usr/bin/env Rscript

#####################################################################################################
#FS19C Roary
#Kathy Mou

#Purpose: Analyze pan-genome results from roary. 

#Load library packages
library(ggplot2)
library(tidyverse)
library(splitstackshape)
library(pheatmap)

sessionInfo()
#R version 4.0.2 (2020-06-22)

#####################################################################################################

#Code taken from: https://github.com/sanger-pathogens/Roary/blob/master/bin/create_pan_genome_plots.R
#It takes in the *.Rtab files and produces graphs on 
#how the pan genome varies as genomes are added (in random orders).

#Number of new genes
mydata = read.table('number_of_new_genes.Rtab')
boxplot(mydata, data=mydata, main="Number of new genes",
         xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

#Number of conserved genes
mydata = read.table("number_of_conserved_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of conserved genes",
          xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

#Number of genes in pan-genome
mydata = read.table("number_of_genes_in_pan_genome.Rtab")
boxplot(mydata, data=mydata, main="No. of genes in the pan-genome",
          xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

#Number of unique genes
mydata = read.table("number_of_unique_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of unique genes",
         xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

#Blast identity frequency
mydata = read.table("blast_identity_frequency.Rtab")
plot(mydata,main="Number of blastp hits with \ndifferent percentage \n identity",  xlab="Blast percentage identity", ylab="No. blast results")

#Number of conserved genes vs total genes
conserved = colMeans(read.table("number_of_conserved_genes.Rtab"))
total = colMeans(read.table("number_of_genes_in_pan_genome.Rtab"))
genes = data.frame( genes_to_genomes = c(conserved,total),
                    genomes = c(c(1:length(conserved)),c(1:length(conserved))),
                    Key = c(rep("Conserved genes",length(conserved)), rep("Total genes",length(total))) )
ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
theme_classic() +
ylim(c(1,max(total)))+
xlim(c(1,length(total)))+
xlab("No. of genomes") +
ylab("No. of genes")+ theme_bw(base_size = 16) +  theme(legend.justification=c(0,1),legend.position=c(0,1))
#ggsave(filename="conserved_vs_total_genes.png", scale=1)

#Number of unique genes vs new genes
unique_genes = colMeans(read.table("number_of_unique_genes.Rtab"))
new_genes = colMeans(read.table("number_of_new_genes.Rtab"))
genes = data.frame( genes_to_genomes = c(unique_genes,new_genes),
                    genomes = c(c(1:length(unique_genes)),c(1:length(unique_genes))),
                    Key = c(rep("Unique genes",length(unique_genes)), rep("New genes",length(new_genes))) )
ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
theme_classic() +
ylim(c(1,max(unique_genes)))+
xlim(c(1,length(unique_genes)))+
xlab("No. of genomes") +
ylab("No. of genes")+ theme_bw(base_size = 16) +  theme(legend.justification=c(1,1),legend.position=c(1,1))
#ggsave(filename="unique_vs_new_genes.png", scale=1)


#####################################################################################################
#Heatmap of gene presence and absence of sugar catabolism genes

# Previous issues with as.matrix()
# 1. I read `./Files/gene_presence_absence.Rtab` with read.table  %>%  column_to_rownames and it couldn't find column "Gene" in data.
# 1b. Also, when I assigned `Gene` as rownames, and removed extra column, as.matrix() wouldn't convert the dataframe to a matrix. It stayed as a dataframe.
# 2. I ran read.csv and read_tsv  %>%  column_to_rownames. When I ran as.matrix() separately, it wouldn't convert to matrix and stayed as a dataframe.
# 3. So it looks like I need to read.csv or read_tsv %>% column_rownames %>% as.matrix() all at the same time and then object stays as a matrix

#read_tsv
PA <- read_tsv('./Files/gene_presence_absence.Rtab') %>%
        column_to_rownames(var = 'Gene') %>% 
        as.matrix()
colnames(PA)
rownames(PA)
pa <- c("1-428RN3A_pol", "10-434FEN3_pol", "11-434FEN3_pol", 
        "12-435FEN3_pol", "13-435FEN3_pol", "14-437FEN5_pol", 
        "15-437FEN5_pol", "57-436REC_pol",  "58-436RED_pol",  
        "59-437REC_pol",  "6-437REN3B_pol", "60-437RED_pol",
        "61-438REC_pol",  "62-438RED_pol", "Ecoli_HS", "Ecoli_K-12_MG1655", 
        "Ecoli_NADC6564", "Ecoli_Nissle1917", "Ecoli_O157H7_EDL933", "Ecoli_TW14588")
keep_rows <- rownames(PA) %>% grep("^ara|^ed|^eut|^fuc|^gal|^man|^nag|^nan|^rbs|^suc|^uxa", .) #narrow down gene list to sugar catabolism genes
PA2 <- PA[keep_rows,pa]
rownames(PA2)
colnames(PA2)
class(PA2) #matrix and array
pheatmap(PA2, cellheight = 6) #increased cell height so easier to read gene names on right side of heatmap
#export as 500 width, 1300 height
pheatmap(PA2)

#Test with read.csv and get same results as read_tsv
PA <-read.csv("./Files/gene_presence_absence.Rtab",  sep = "\t", quote = "") %>% 
        column_to_rownames(var="Gene") %>% 
        as.matrix()
class(PA)
pa <- c("1-428RN3A_pol", "10-434FEN3_pol", "11-434FEN3_pol", 
            "12-435FEN3_pol", "13-435FEN3_pol", "14-437FEN5_pol", 
            "15-437FEN5_pol", "57-436REC_pol",  "58-436RED_pol",  
            "59-437REC_pol",  "6-437REN3B_pol", "60-437RED_pol",
            "61-438REC_pol",  "62-438RED_pol", "EDL933", 
            "MG1655", "NADC6564", "Nissle1917", "TW14588")

#sugar genes
keep_rows <- rownames(PA) %>% grep("^ara|^ed|^eut|^fuc|^gal|^lac|^man|^nag|^nan|^rbs|^suc|^uxa", .) #narrow down gene list to sugar catabolism genes
PA2 <- PA[keep_rows,pa]
rownames(PA2)
colnames(PA2)
pheatmap(PA2, cellheight = 6)

#stx genes
stx <- rownames(PA) %>% grep("^stx", .)
PA3 <- PA[stx,pa]
pheatmap(PA3, cellheight = 6)

##### Jules Zone #####
library(tidyverse)

PA <- read_tsv('./Files/gene_presence_absence.Rtab') %>%
        column_to_rownames(var = 'Gene') %>% 
        as.matrix()


PA2 <- PA[keep_rows,pa]
rownames(PA2)
colnames(PA2)

pheatmap(PA2)
