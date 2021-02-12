#!/usr/bin/env Rscript

#####################################################################################################
#FS19C Roary
#Kathy Mou

#Purpose: Analyze pan-genome results from roary. It takes in the *.Rtab files and produces graphs on 
#how the pan genome varies as genomes are added (in random orders).
#Code taken from: https://github.com/sanger-pathogens/Roary/blob/master/bin/create_pan_genome_plots.R

#Load library packages
library(ggplot2)
library(tidyverse)

sessionInfo()
#R version 4.0.2 (2020-06-22)

#####################################################################################################

mydata = read.table('number_of_new_genes.Rtab')
boxplot(mydata, data=mydata, main="Number of new genes",
         xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_conserved_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of conserved genes",
          xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_genes_in_pan_genome.Rtab")
boxplot(mydata, data=mydata, main="No. of genes in the pan-genome",
          xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_unique_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of unique genes",
         xlab="No. of genomes", ylab="No. of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("blast_identity_frequency.Rtab")
plot(mydata,main="Number of blastp hits with \ndifferent percentage \n identity",  xlab="Blast percentage identity", ylab="No. blast results")

##################

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

######################

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

######################
presenceabsence <-read.table("gene_presence_absence.Rtab", header = TRUE, sep = "\t",
                    quote = "")
presenceabsence[1:5,1:5] # 1 = present, 0 = absent
rownames(presenceabsence) <- presenceabsence$Gene # set Gene as rownames
presenceabsence <- presenceabsence[,-1] # remove extra Gene column

#pull out eut-only genes
presenceabsence_eut <-presenceabsence[grep('^eut', rownames(presenceabsence)),]
names(which(colSums(presenceabsence_eut == 1) > 0))
names(which(rowSums(presenceabsence_eut == 1) > 0))


#pull out aaeA genes as a test
presenceabsence_aae <-presenceabsence[grep('^aae', rownames(presenceabsence)),]
names(which(colSums(presenceabsence_aaeA == 1) > 0))
names(which(rowSums(presenceabsence_aaeA == 1) > 0))
