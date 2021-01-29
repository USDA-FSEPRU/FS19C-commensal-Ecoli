# Methods

Written summary of methods performed in this repo. Includes lab notes of how methods performed

## Sample Collection
See OneNote FS19C lab notebook entries xxx
* FS19C Samples 1-96 Final Data.xlsx
* Hannah Sorbitol-positive isolates - MALDI, list for sequencing.xlsx
* Sorbitol-negative isolates - agglutination, MALDI, list for sequencing.xlsx
* KathyMou_NovaSeq_Submission_Form_8June2020.xlsx
* FS19C_metadata.xlsx


## Conda environments
### fastanienv
* FastANI, multiqc, prokka

### mashenv
* Mash
  * Directions: http://mash.readthedocs.org

### prokka_env
* prokka

## Sequence Analysis
### Sequence assembly
* Summary:
* Completed on:
* Platform:
* Commands:
1. Make text file replace.sh
```
#usr/bin/bash
while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
```

2. Run replace.sh
```
sh replace.sh
```

3. Rename *.fastq.gz files to fit the '*_1.fastq.gz' or '*_2.fastq.gz' format using the following bash script to create a directory "linked", take the fastq.gz files in current directory, shorten the names from something like 1-H12-96-441FEC_S103_L001_R2_001.fastq.gz to 1-H12-96-441FEC_2.fastq.gz, and copy over the short-named fastq.gz files to "linked" directory while also linking the short-named fastq.gz files to the original long-named fastq.gz files. Saved as bash script "renamefiles.sh" and ran on ceres.
```
mkdir linked
for myfile in ./*.fastq.gz; do
        target=$(echo $myfile|sed -r 's/1\-[A-H][0-9]+\-([0-9]+\-[0-9]+[A-Z]+[0-9]?[A-Z]?\_)[A-Z]+[0-9]+\_[A-Z]+[0-9]+\_R([0-9])\_[0-9]+/\1\2/g')
        echo "$myfile is being renamed to $target"
        ln -sr $myfile linked/$target
done
```


4. Run the follow slurm script (SRAassemblyPipeline.FS19C.SLURM_TEMPLATE)
```
#!/bin/bash
#SBATCH --job-name=REPLACE                              # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 12:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
```

```
set -e
module load java
module load pilon
module load samtools
module load r

#source ~/templates/adapt_polish.sh

#convert to interleaved
reformat.sh in1=REPLACE_1.fastq.gz in2=REPLACE_2.fastq.gz out=REPLACE.fq.gz
ln -s REPLACE.fq.gz REPLACE_temp.fq.gz

#diagnose interleaving
#bbsplitpairs.sh in=20-427FEC.fq.gz out=fixed.20-427FEC.fq.gz outs=singletons.20-427FEC.fq.gz minlen=70 fint
#bbsplitpairs.sh didn't produce any output
100X: repair.sh in=old_20-427FEC.fq.gz out=20-427FEC_temp.fq.gz outs=nonint_singletons.20-427FEC.fq.gz repair
250X: repair.sh in=20-427FEC.fq.gz out=20-427FEC_temp.fq.gz outs=nonint_singletons.20-427FEC.fq.gz repair

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=REPLACE_temp.fq.gz out=REPLACE_trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=/home/kathy.mou/software/bbmap/resources/adapters.fa ftm=5 ordered interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_trimmed.fq.gz REPLACE_temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=REPLACE_temp.fq.gz out=REPLACE_filtered.fq.gz k=31 ref=/home/kathy.mou/software/bbmap/resources/sequencing_artifacts.fa.gz,/home/kathy.mou/software/bbmap/resources/phix174_ill.ref.fa.gz ordered cardinality interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_filtered.fq.gz REPLACE_temp.fq.gz

#subsample to approx 100x
reformat.sh in=REPLACE_temp.fq.gz sbt=600000000 (100X) /  1500000000 (250X) out=REPLACE_subsamp.fq.gz interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_subsamp.fq.gz REPLACE_temp.fq.gz

#change sbt=600000000 (page 19 of https://www.fsis.usda.gov/wps/wcm/connect/e0b56754-662e-4043-91ab-f4e84aff169f/mlg-42.pdf?MOD=AJPERES)

#Error-correct phase 1
bbmerge.sh in=REPLACE_temp.fq.gz out=REPLACE_ecco.fq.gz ecco mix vstrict ordered interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_ecco.fq.gz REPLACE_temp.fq.gz

#Error-correct phase 2
clumpify.sh in=REPLACE_temp.fq.gz out=REPLACE_eccc.fq.gz ecc passes=4 reorder interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_eccc.fq.gz REPLACE_temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
tadpole.sh in=REPLACE_temp.fq.gz out=REPLACE_ecct.fq.gz ecc k=62 ordered tossjunk=t interleaved=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_ecct.fq.gz REPLACE_temp.fq.gz

#Assemble with Spades
spades.py --12 REPLACE_ecct.fq.gz -o REPLACE_spades_out --only-assembler -t 16 -m 32 -k 25,55,95,125

# --- Evaluation ---

#Calculate the coverage distribution, and capture reads that did not make it into the assembly
bbmap.sh in=REPLACE_subsamp.fq.gz ref=REPLACE_spades_out/scaffolds.fasta nodisk covstats=REPLACE_covstats.txt maxindel=200 minid=90 qtrim=10 untrim ambig=all interleaved=t

#Calculate assembly stats
stats.sh in=REPLACE_spades_out/scaffolds.fasta

#Remove bad contigs and generate list of good contigs
Rscript ~/software/good_contig_names.R REPLACE_covstats.txt
#rm REPLACE_covstats.txt

filterbyname.sh names=REPLACE.names in=REPLACE_spades_out/scaffolds.fasta out=REPLACE.fasta include=t
#rm REPLACE.names

# if this all finishes correcly then remove all the intermediates #
#rm -r REPLACE_spades_out/
#rm REPLACE_eccc.fq.gz REPLACE_ecco.fq.gz REPLACE_filtered.fq.gz #REPLACE_trimmed.fq.gz REPLACE_1.fastq.gz REPLACE_2.fastq.gz #REPLACE.fq.gz REPLACE_ecct.fq.gz REPLACE_temp.fq.gz

# run polishing script #
source ~/software/adapt_polish.sh

adapt_polish REPLACE.fasta REPLACE_subsamp.fq.gz 4

#rm REPLACE.SLURM
```


### QC
#### FastQC
* Summary: Ran fastqc on FS19C sequence data to assess sequence quality of individual reads for each sample, and to use output for multiqc (forgot to do this prior to sequence assembly)
* Completed on: 7Jan2021
* Platform: Ceres, slurm
* fastqc.batch slurm script (aka fastqc.slurm):

  ```
  #!/bin/bash
  #SBATCH --job-name=fastqc                              # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 12:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=32G   # memory
  #Enter commands here:
  module load fastqc
  fastqc -t 16 *.fastq.gz  
  mv *fastqc* ./fastqc/
  #End of file
  ```

* Files generated:
  * *fastqc.zip
  * *fastqc.html

#### MultiQC
* Tutorial: https://www.youtube.com/watch?v=qPbIlO_KWN0
* Summary: Ran multiqc with fastqc output of FS19C sequence data to assess quality of all sequences for samples 1-94 (forgot to do this prior to sequence assembly). I did not include samples 95 and 96 in the multiqc run as I realized these samples were ran on a second sequencing run, so their coverage is different than the first sequencing run that had all samples (but I would only consider  samples 1-94 from that first run). In addition, examining fastqc output for samples 95 and 96 was fine enough.
* Completed on: 7Jan2021
* Platform: Ceres, conda environment
* Commands:

  ```
  salloc
  module load miniconda
  source activate fastanienv
  conda install -c bioconda multiqc
  multiqc *.fastqc.zip
  conda deactivate
  ```

1. Why are the plots flat plot (not interactive)?
  * From: https://multiqc.info/docs/
  "Flat plots
  Reports with large numbers of samples may contain flat plots. These are rendered when the MultiQC report is generated using MatPlotLib and are non-interactive (flat) images within the report. The reason for generating these is that large sample numbers can make MultiQC reports very data-intensive and unresponsive (crashing people's browsers in extreme cases). Plotting data in flat images is scalable to any number of samples, however.
  Flat plots in MultiQC have been designed to look as similar to their interactive versions as possible. They are also copied to multiqc_data/multiqc_plots"
2. I noticed the FS19all_multiqc_report.html report showed samples 95 and 96 having huge number of reads, per sequence GC content had several large peaks, and samples 95 and 96 making up the majority of the overrepresented sequences.
3. To try to eliminate the discrepancies due to samples 95 and 96 (these two samples had so many reads as a result of re-sequencing them) and are therefore skewing the report stats, I moved samples 95 and 96 fastqc.zip to Fastqc_Sample95_96 directory, and ran multiqc to generate FS19_1-94_multiqc_report.html and FS19_1-94_data directory.
4. multiqc report of samples 1-94 look a lot better with the sequence count ranges being a lot closer among all samples, sequence quality histograms all in green, per sequence quality scores in green, less than 1% of reads making up overrepresented sequences, and a single peak for per sequence GC content
5. I also looked at the fastqc reports for samples 95 and 96 individually.
  * **95** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9
  * **96** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9
6. Note: Jules says if you're able to get assemblies to work, you can be sure that the sequences' qualities are good

* Files generated:
 * FS19all_multiqc_report.html
 * FS19all_multiqc_data directory
 * FS19_1-94_multiqc_report.html
 * FS19_1-94_multiqc_data directory
 * 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
 * 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R1_001_fastqc.html

#### FastANI
* Summary: ran fastANI on FS19C sequence data by running in conda environment to estimate Average Nucleotide Identity (ANI) using alignment-free approximate sequence mapping. It calculates distance between 2 sequences. Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* FastANI publication: DOI: 10.1038/s41467-018-07641-9
* Completed on: 29Dec2020
* Platform: Ceres, conda
* Commands ran:
1. Preparing conda environment on Ceres
```
salloc
module load miniconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name fastanienv
source activate fastanienv
conda install -c bioconda fastani
```

2. Testing fastani on sample genomes provided by https://github.com/ParBLiSS/FastANI
```
fastANI -q Shigella_flexneri_2a_01.fna -r Escherichia_coli_str_K12_MG1655.fna -o testfastani.out
```
Viewed output testfastani.out and it looked the same as output from https://github.com/ParBLiSS/FastANI. Awesome!

3. Made query list in polished_genomes/ with this command and checked to make sure query list only contains sample names (tells the path to where reference genome files are):
```
ls -dv "$PWD"/* > quertylist.txt
```

4. Made reference list with this (tells the path to where reference genome files are):
```
ls -dv "$PWD"/* > referencelist.txt
```
Reference genomes include:
* E. coli MG1655 aka https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3
  * Saved as Ecoli_K-12_MG1655.fasta
* E. coli HS aka https://www.ncbi.nlm.nih.gov/nuccore/NC_009800.1
  * Saved as Ecoli_HS.fasta
* E. coli Nissle 1917 aka https://www.ncbi.nlm.nih.gov/nuccore/CP007799.1
  * Saved as Ecoli_Nissle1917.fasta
* E. coli O157:H7 str. NADC 6564 aka https://www.ncbi.nlm.nih.gov/nuccore/CP017251
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6408774/#CR28
  * Saved as Ecoli_NADC6564.fasta
* E. coli O157:H7 EDL933 aka https://www.ncbi.nlm.nih.gov/nuccore/CP008957.1
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6408774/
  * Saved as Ecoli_O157H7_EDL933.fasta



5. Moved querylist.txt and referencelist.txt to same directory and ran fastani:
```
fastANI --ql querylist.txt --rl referencelist.txt -o fs19cfastanioutput.out
```

6. Downloaded fs19cfastanioutput.out and made FS19CfastANIoutput.xlsx

7. After discussion with Jules about how I'd run mash, I realized I ran fastani incorrectly. I should be measuring distances of all genomes (references and samples) between each other. So combine all file names in a new querylist2.txt and run fastani with this file listed as query and as reference lists.
```
fastANI --ql querylist2.txt --rl querylist2.txt -o fs19cfastanioutput2.out
```

* Files generated:
  * FS19CfastANIoutput.xlsx
  * FS19CfastANIoutput2.xlsx
  * fs19cfastanioutput.out
  * fs19cfastanioutput2.out

#### Mash
* Summary: ran Mash on FS19C sequence data by running in conda environment to compare results with fastANI. The *sketch* function converts a sequence or collection of sequences into a MinHash sketch. The *dist* function compares two sketches and returns an estimate of the Jaccard index, *P* value, and Mash distance (estimates rate of sequence mutation under a simple evolutionary model). Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* Mash publication: DOI: 10.1186/s13059-016-0997-x
* Source: http://mash.readthedocs.org
* Completed on: 11Jan2021
* Platform: Ceres, conda
* Commands ran:
1. Load environment and install mash
```
salloc
module load miniconda
conda create --name mashenv
source activate mashenv
conda install -c bioconda mash
```

2. Make a new directory mash_all that has all polished genomes and reference genomes together. Then sketch all the genomes together.
```
mash sketch -p 1 -o /project/fsepru/kmou/FS19C/polished_genomes_100X/mash_all/ *.fasta
```
Reference genomes include:
* E. coli MG1655 aka https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3
  * Saved as Ecoli_K-12_MG1655.fasta
* E. coli HS aka https://www.ncbi.nlm.nih.gov/nuccore/NC_009800.1
  * Saved as Ecoli_HS.fasta
* E. coli Nissle 1917 aka https://www.ncbi.nlm.nih.gov/nuccore/CP007799.1
  * Saved as Ecoli_Nissle1917.fasta
* E. coli O157:H7 str. NADC 6564 aka https://www.ncbi.nlm.nih.gov/nuccore/CP017251
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6408774/#CR28
  * Saved as Ecoli_NADC6564.fasta
* E. coli O157:H7 EDL933 aka https://www.ncbi.nlm.nih.gov/nuccore/CP008957.1
  * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6408774/
  * Saved as Ecoli_O157H7_EDL933.fasta

Output generated: .msh

3. Run distance estimates comparing sketch genomes with itself. This runs really fast (less than a second)
```
mash dist -p 1 .msh .msh > distances.tab
```

4. distances.tab columns: reference_id, query_id, mash_distance, pvalue, matching_hashes

* Files generated:
  * distances.tab
  * FS19Cmashdistances.xlsx


#### Visualize ANI pairwise genome-genome similarity calculations with MDS, heatmap
* Summary: Made distance matrix from mash and fastani output to create heatmap and MDS to visualize clustering and identify any outliers. The MDS was a bit hard to decipher what was an outlier, so I ran a heatmap to see how fastANI and mash compared and whether the pairwise comparisons were similar between the two, including heatmap of pearson correlation coefficients.
* Completed on: 15Jan2021
* Platform: R
* Commands ran:
1. See scripts/qc_mds.R for details
```
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
```

2. Analyzed fastANI and mash MDS plots. The plots look very different and clustering was variable between the two. Will need to ask Jules what counts as outliers in each plot, and does it matter which MDS plot I choose for determining what counts as outliers before doing pangenome analysis as long as I make note of what I chose?

3. Before asking Jules, I googled how to determine outliers in MDS plot and came across this paper: https://doi.org/10.1093/bioinformatics/btz964. They use heatmaps of distance matrices to type bacteria, so I decided to try generating a heatmap to see how it compares between fastANI and mash.

4. Followed this site for generating heatmap from distance matrix and correlation: https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
* This paper does something similar that ran mash and tested linear correlation with pearson's correlation coefficient test. Purpose was to look at whole genome similarity to identify bacterial meningitis causing species: https://doi.org/10.1186/s12879-018-3324-1

5. Added heatmap and Pearson's correlation coefficient scripts to scripts/qc_mds.R.

6. I compared FS19C_fastani_heatmap.tiff and FS19C_mash_heatmap.tiff and they're very similar (order of samples listed on bottom and side of heatmap aren't in exact order, but they're in the same color blocks). The scale of key is different because I notice that the ANI values used in mash are 1.0 and below while for fastani it's 100 and below (difference by a factor of 100).

7. When I compared FS19C_fastani_correlation_heatmap.tiff with FS19C_mash_correlation_heatmap.tiff, the figures were quite different.

8. I had short video meeting with Jules on 15Jan2021 from 1-1:30pm to go over my MDS plots.
  1. I made two main errors:
    1. Doing distance calculation on mash file when I didn't need to (it already did the calculation for me, it outputs the distances, I just need to make and run distance matrix).
    2. For fastANI file, the values are reported as similarity values in percent (ANI % 1-100). I needed to convert the percent to a proportion, then do 1-X from that decimal to get distance matrix. Then run distance calculation. I originally made a similarity matrix.
  2. Other things to note:
    1. Jules said that changing order of sample names in reference_id shouldn't change ordination. After correcting my errors (see below), I ran mash_tab with and without change to reference_id order and I got the same figures either way. So, I don't need to worry about order of samples listed in distance matrix.
    2. How to determine outliers in MDS? Put another organism (like Salmonella) in the mix that you know is an absolute outlier and re-run ordination to see where the true outliers lie.

9. I ran qc_mds.R script with the corrections and saw that there weren't any obvious outliers, like what Jules had hinted when he ran the code on his end.

10. Next step: run prokka. I can use UnitProt pangenome E. coli to run annotation. I asked if I chould do this or run a single reference genome for annotation and Jules said I should try the pangenome method first. Will need to find the E. coli pangenome annotation (protein).

#### Genome annotation
* Summary: Identify annotation reference from UniProt (E. coli pangenome) and use prokka to annotate all 95 samples.
* Began on: 20Jan2021
* Completed on:
* Platform:
1. Find E. coli pangenome (pan proteome) on UniProt:
  * [Escherichia coli (strain K12) (Strain: K12 / MG1655 / ATCC 47076)](https://www.uniprot.org/proteomes/UP000000625)
  * [Escherichia coli O157:H7 (Strain: O157:H7 / Sakai / RIMD 0509952 / EHEC)](https://www.uniprot.org/proteomes/UP000000558)
2. Read Prokka paper
  * DOI: 10.1093/bioinformatics/btu153
3. Find prokka github page: https://github.com/tseemann/prokka
4. install prokka in fastanienv conda environment on Ceres in FS19C/polished_genomes_100X
```
salloc
module load miniconda
source activate fastanienv
conda install -c conda-forge -c bioconda -c defaults prokka
prokka #this test didn't work...
```
  * How else to install and run prokka on conda??
  * See https://github.com/tseemann/prokka/issues/448
  * Try this: https://github.com/tseemann/prokka/issues/508
5. (21Jan2021) Tried the following commands based on [this issues page](https://github.com/tseemann/prokka/issues/508) on Ceres to remove prokka_env, make new prokka_env with different way of installing prokka. It works!
```
conda env remove --name prokka_env #remove prokka_env and try it again with the next command
conda create -n prokka_env -c conda-forge -c bioconda prokka
source activate prokka_env
prokka
prokka --version 1.14.6
prokka --listdb
#Looking for databases in: /home/kathy.mou/.conda/envs/prokka_env/db
#Kingdoms: Archaea Bacteria Mitochondria Viruses
#Genera: Enterococcus Escherichia Staphylococcus
#HMMs: HAMAP
#CMs: Archaea Bacteria Viruses
```
6. Made a new directory polishedgenomesprokka/, copied fasta files from polished_genomes_100X to this directory, renamed .fasta to .fna
```
rename .fasta .fna *.fasta
```
7. Download [E. coli strain K12 proteome (fasta)](https://www.uniprot.org/proteomes/UP000000625) - Try this annotation first.
8. Need to ask Jules what kind of fasta file do I need to use for annotation first. The fasta file I have doesn't have the ~~~ symbols that prokka says is needed for annotation tag formats: https://github.com/tseemann/prokka/blob/master/README.md#fasta-database-format
9. (22Jan2021) At lab meeting, Jules recommended not spending too much time finding the "perfect" annotation (custom or publicly-available ones) because there will always be something missing. Instead, find annotations that have the genes you want. So I looked up what genes the project plan had mentioned (see 01_Background) and came up with this list:
  * EA utilization, eut operon for utilization ethanolamine as a nitrogen source (see project plan)
  * genes involved in utilizing the following carbon/sugar substrates: glucose, sucrose, galactose, arabinose, lactose, fucose, maltose, hexuronate, mannose, ribose, N-acetylglucosamine, N-acetylgalactosamine, N-acetylneuraminate, sialic acid and D-gluconate.
	* bacteriocins
  * Look up papers from project plan, record in 01_Background.md, read and find exact gene or operon names
10. (27Jan2021)  Look through UniProt for metabolic pathway genomes
  * Complete list: ethanolamine, glucose, sucrose, galactose, arabinose, lactose, fucose, maltose, hexuronate, mannose, ribose, N-acetylglucosamine, N-acetylgalactosamine, N-acetylneuraminate, sialic acid and D-gluconate
  * List of substrate-related genes not in [Escherichia coli K12 substr. MG1655 EMBL file](https://www.ebi.ac.uk/ena/browser/view/U00096): missing N-acetylgalactosamine
      * File: U00096.3.txt
  * List of substrate-related genes not in [Escherichia coli O157:H7 str. Sakai EMBL file](https://www.ebi.ac.uk/ena/browser/view/BA000007): it has genes related to all substrates (not sure if the complete metabolic pathways are present for each substrate, but there are proteins associated with all these substrates)
      * File: BA000007.3.txt
11. Looked through each EMBL file and print out all proteins associated with each of the substrates:
  * Escherichia coli O157:H7 str. Sakai
    ```
    grep -f carbohydratelist.txt BA000007.3.txt > completecarblistO157H7.txt
    ```
  * Escherichia_coli_str_K12_MG1655
    ```
    grep -f carbohydratelist.txt U00096.3.txt > completecarblistK12.txt
    ```
12. Compared the two files and printed out genes unique to each strains
  ```
  comm -1 -3 completecarblistO157H7.txt completecarblistK12.txt > genesnotpresentinO157H7.txt #not present in O157H7
  comm -1 -3 completecarblistK12.txt completecarblistO157H7.txt > genesnotpresentinK12.txt
  ```
13. So, at least one of the strains' proteomes have genes associated with each of the metabolites listed in OSQR plan. I will also go through review papers referenced in OSQR plan to see what other genes to include.
 * Maltby et al.: EMBL files missing glucuronate, galacuronate; however, based on a paper referenced by Maltby et al., glucoronate, galacuronate, and hexuronate are one sugar
 *    
14. (28Jan2021) Decide to run prokka on Ceres using the for loop to run prokka on 95 genomes, taken from: https://doi.org/10.3389/fvets.2020.582297
```
salloc
module load miniconda
source activate prokka_env
for file in *.fna; do tag=$file%.fna; prokka –prefix “$tag” –locustag “$tag” –genus Escherichia –strain “$tag” –outdir “$tag”_prokka –force –addgenes “$file”; done
```
Get error message:
```
[09:45:50] This is prokka 1.14.6
[09:45:50] Written by Torsten Seemann <torsten.seemann@gmail.com>
[09:45:50] Homepage is https://github.com/tseemann/prokka
[09:45:50] Local time is Thu Jan 28 09:45:50 2021
[09:45:50] You are kathy.mou
[09:45:50] Operating system is linux
[09:45:50] You have BioPerl 1.007002
[09:45:50] System has 72 cores.
[09:45:50] Will use maximum of 8 cores.
[09:45:50] Annotating as >>> Bacteria <<<
[09:45:50] '–prefix' is not a readable non-empty FASTA file
```
15. Looked up issues page on prokka: https://github.com/tseemann/prokka/issues/86 and checked that my files are fasta format, they are readable, and are more than 0 bytes. Emailed Jules how to fix error. He says it looks like my prokka command has some weird formatting applied and some special characters were inserted where they shouldn't be. This is a danger when copying things from the internet. There's a difference between long and short dash. Shell scripts don't take too kindly to things like long dashes and it confuses them, along with other hidden special characters. So I retyped the command on bbedit, and prokka ran past the part I got stuck on in the last step. However, it didn't like how long the length of contig IDs in my fasta are, so it suggested renaming contigs or try --centre X --compliant to generate clean contig names. Need to shorten them to less than or equal to 37 characters long.
```
[10:20:41] Contig ID must <= 37 chars long: NODE_1_length_378381_cov_16.779345_pilon_pilon_pilon
[10:20:41] Please rename your contigs OR try '--centre X --compliant' to generate clean contig names.
```
Added --centre X --compliant command line options
```
for file in *.fasta; do tag=$file%.fasta; prokka -prefix "$tag" -locustag "$tag" -genus Escherichia -strain "$tag" -outdir "$tag"_prokka -force -addgenes "$file" -centre X -compliant; done
```
It is working! Started around 10:30AM, finished at 11pm. All samples have a new directory with .err, .faa, .fnn, .fna, .fsa, .gbk, .gff, .log, .sqn, .tbl, .tsv, .txt files

16. (29Jan2021) Read Roary paper and look up how to run roary


## WGS submission to SRA
* Must complete Biosample entry before you can complete and Bioproject entry
* [SRA Quick Start Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submit/) and a more detailed [SRA Submission Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submitbio/)
