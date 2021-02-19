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
presenceabsence <-read.table("./Files/gene_presence_absence.Rtab", header = FALSE, sep = "\t",
                    quote = "")
colnames(presenceabsence) <- presenceabsence[1,] #set first row as column name
presenceabsence <- presenceabsence[-1,] # remove extra row
presenceabsence[1:5,1:5] # 1 = present, 0 = absent
colnames(presenceabsence)
rownames(presenceabsence) <- presenceabsence[,1] #set first column as row name
presenceabsence <- presenceabsence[,-1]
pa <- c("1-428RN3A_pol.fasta%.fasta", "10-434FEN3_pol.fasta%.fasta", "11-434FEN3_pol.fasta%.fasta", 
        "12-435FEN3_pol.fasta%.fasta", "13-435FEN3_pol.fasta%.fasta", "14-437FEN5_pol.fasta%.fasta", 
        "15-437FEN5_pol.fasta%.fasta", "57-436REC_pol.fasta%.fasta",  "58-436RED_pol.fasta%.fasta",  
        "59-437REC_pol.fasta%.fasta",  "6-437REN3B_pol.fasta%.fasta", "60-437RED_pol.fasta%.fasta",
        "61-438REC_pol.fasta%.fasta",  "62-438RED_pol.fasta%.fasta", "EDL933.fasta%.fasta", 
        "MG1655.fasta%.fasta", "NADC6564.fasta%.fasta", "Nissle1917.fasta%.fasta", "TW14588.fasta%.fasta")
pa1 <- presenceabsence[pa] #select only the isolates from pa
keep_rows <- rownames(pa1) %>% grep("^ara|^ed|^eut|^fuc|^gal|^man|^nag|^nan|^rbs|^suc|^uxa", .) #narrow down gene list to sugar catabolism genes
pa2 <- pa1[keep_rows,]
pheatmap(pa2) #Error during wrapup: 'x' must be numeric. 
#Tried a few different solutions from various forums, got an error message that list object cannot be coerced to type double. Checked out this site:
#https://www.programmingr.com/r-error-messages/list-object-cannot-be-coerced-to-type-double/
pa2_num <- unlist(pa2) #convert list to single vector
pa2_numeric <- lapply(pa2_num, as.numeric) #convert to numeric
as.matrix(pa2_numeric) #convert to matrix
pheatmap(pa2_numeric,mean="pheatmap default") #Error during wrapup: must have n >= 2 objects to cluster





##### Jules Zone #####
library(tidyverse)

PA <- read_tsv('./Files/gene_presence_absence.Rtab') %>%
        column_to_rownames(var = 'Gene') %>% 
        as.matrix()


PA2 <- PA[keep_rows,pa]
rownames(PA2)
colnames(PA2)

pheatmap(PA2)
