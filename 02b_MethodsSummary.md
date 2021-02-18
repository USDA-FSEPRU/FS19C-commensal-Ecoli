# Methods Summary

Summary of sequence analyses methods performed on FS19C samples 1-96 in this repo (has some information repeated in 02_Methods.md). This summary skips the errors I encountered in 02_Methods.md and only shows the final corrected code + procedures for each of the tools used.

## Sample Collection
### OneNote FS19C lab notebook entries on E. coli isolates 1-96 from 24Apr2020 to 11June2020: streaking plates, DNA extraction, re-extraction of DNA, DNA cleanup, nanodrop, Qubit, gel
* **20May2020** DNA extraction of isolates #1-24 using DNeasy Blood and Tissue Kit
* **14May2020** DNA extraction of isolates #25-48 using DNeasy Blood and Tissue Kit
* **21May2020** DNA extraction of isolates #49-72 using DNeasy Blood and Tissue Kit
* **28May2020** DNA extraction of isolates #73-96 using DNeasy Blood and Tissue Kit
* Lists of samples had additional DNA re-extractions found on entries: 28May2020, 29May2020, 3June2020, 11June2020
* **11June2020** Submission of FS19C 96-well plate of E. coli isolates 1-96 DNA to David Alt for NovaSeq sequencing.

### Additional files of importance
* FS19C Samples 1-96 Final Data.xlsx
* Hannah Sorbitol-positive isolates - MALDI, list for sequencing.xlsx
* Sorbitol-negative isolates - agglutination, MALDI, list for sequencing.xlsx
* KathyMou_NovaSeq_Submission_Form_8June2020.xlsx
* FS19C_metadata.xlsx
* FS19C 96 S-S+ E. coli gDNA gels.pdf

## Conda environments
### fastanienv
* FastANI, multiqc, prokka

### mashenv
* Mash
  * Directions: http://mash.readthedocs.org

### prokka_env
* prokka, PPanGGOLiN

## Sequence assembly
* Summary: QC and assemble sequences using BBMap and spades
* Platform: Ceres

1. Make text file replace.sh
```
#usr/bin/bash
while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
```

2. Run replace.sh to make slurm script for each of the 96 isolates.
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

4. Run the slurm script (SRAassemblyPipeline.FS19C.SLURM_TEMPLATE) to assemble sequences for each isolate (e.g. 10-434FEN3.slurm, 65-440REC.slurm, etc.) except for isolate #20 (assembly failed).
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

##### Files/directories generated (for each isolate if indicated with a '*')
* polishedfasta.txt
* *_1.fastq.gz
* *_2.fastq.gz
* *_pol.fasta
* *.slurm
* *.names
* *_covstats.txt
* *.fasta
* *_spades_out/

## QC with FastQC
* Summary: Ran fastqc on FS19C sequence data to assess sequence quality of individual reads for each sample, and to use output for multiqc
* Platform: Ceres, slurm

1. Run the following fastqc.batch slurm script (aka fastqc.slurm) on Ceres:
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

2. I looked at the fastqc reports for samples 95 and 96 individually (the re-sequenced isolates 95 and 96).
* **95** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9
  * 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
  * 1-H11-95-440FED_S1_L002_R1_001_fastqc.html
* **96** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9
  * 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
  * 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html

##### Files generated (for each isolate):
  * *fastqc.zip
  * *fastqc.html

## QC with MultiQC
* Tutorial: https://www.youtube.com/watch?v=qPbIlO_KWN0
* Summary: Ran multiqc with fastqc output of FS19C sequence data to assess quality of all sequences for samples 1-94. I did not include samples 95 and 96 in the multiqc run as these samples were ran on a second sequencing run, so their coverage is different than the first sequencing run that had all samples (Samples 1-94 are from the first run). In addition, examining fastqc output for samples 95 and 96 was fine enough.
* Platform: Ceres, fastanienv conda environment

1. Command ran on Ceres to run fastanienv conda environment:
  ```
  salloc
  module load miniconda
  source activate fastanienv
  conda install -c bioconda multiqc
  multiqc *.fastqc.zip
  conda deactivate
  ```

2. Ran multiqc on isolates #1-94 to generate FS19_1-94_multiqc_report.html and FS19_1-94_data directory.

##### Files generated:
 * FS19all_multiqc_report.html
 * FS19all_multiqc_data directory
 * FS19_1-94_multiqc_report.html
 * FS19_1-94_multiqc_data directory


## QC with FastANI
* Summary: ran fastANI on FS19C sequence data by running in conda environment to estimate Average Nucleotide Identity (ANI) using alignment-free approximate sequence mapping. It calculates distance between 2 sequences. Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* FastANI publication: DOI: 10.1038/s41467-018-07641-9
* Platform: Ceres, fastanienv conda environment

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

1. List of reference genomes to include in fastANI analysis:
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

2. Combine all file names of 95 isolates and 5 reference genomes in a querylist2.txt and run fastani with this file listed as query and as reference lists.
```
ls -dv "$PWD"/* > quertylist2.txt
fastANI --ql querylist2.txt --rl querylist2.txt -o fs19cfastanioutput2.out
```

##### Files generated:
  * FS19CfastANIoutput.xlsx
  * FS19CfastANIoutput2.xlsx
  * fs19cfastanioutput.out
  * fs19cfastanioutput2.out

## QC with Mash
* Summary: ran Mash on FS19C sequence data by running in conda environment to compare results with fastANI. The *sketch* function converts a sequence or collection of sequences into a MinHash sketch. The *dist* function compares two sketches and returns an estimate of the Jaccard index, *P* value, and Mash distance (estimates rate of sequence mutation under a simple evolutionary model). Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* Mash publication: DOI: 10.1186/s13059-016-0997-x
* Source: http://mash.readthedocs.org
* Platform: Ceres, mashenv conda environment

1. Install and load mash in mashenv conda environment
```
salloc
module load miniconda
conda create --name mashenv
source activate mashenv
conda install -c bioconda mash
```

2. Make a new directory `mash_all` that has all polished genomes and reference genomes together. Then sketch all the genomes together.
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
* * (Enterobacteriaceae) Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000006945.2
  * Saved as Styphimurium_LT2.fna.gz
* (Non-Enterobacteriaceae, same phylum) Campylobacter jejuni subsp. jejuni NCTC 11168 aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000009085.1
  * Saved as Cjejuni_11168.fna.gz
* (Non-Proteobacteria, a Firmicutes) Clostridium saccharoperbutylacetonicum N1-4(HMT) aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000340885.1
  * Saved as Clostridium_N1-4.fna.gz

3. Run distance estimates comparing sketch genomes with itself. This runs really fast (less than a second)
```
mash dist -p 1 .msh .msh > distances.tab
```
* distances_thirdrun.tab columns: reference_id, query_id, mash_distance, pvalue, matching_hashes

##### Files generated:
  * distances_thirdrun.tab
  * FS19Cmashdistances_thirdrun.xlsx

## QC: Visualize ANI pairwise genome-genome similarity calculations with MDS, heatmap
* Summary: Made distance matrix from mash and fastani output to create heatmap and MDS to visualize clustering and identify any outliers. The MDS was a bit hard to decipher what was an outlier, so I ran a heatmap to see how fastANI and mash compared and whether the pairwise comparisons were similar between the two, including heatmap of pearson correlation coefficients.
* Platform: R Studio

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

    ######## Make distance matrix ########
  #mash_dist with mash_tab
  mash_dist <- mash_tab %>%
  arrange(reference_id) %>%
  pivot_wider(names_from = reference_id, values_from=ani_dist) %>%
  column_to_rownames(var='query_id')

    ######## Heatmap of mash distance matrix as another way of looking at clustering ########
  #using mash_dist
  mash_heatmap <- pheatmap(mash_dist, scale = "row", main = "mash distance matrix heatmap")
  mash_heatmap
  ggsave("FS19C_mash_heatmap.tiff", plot=mash_heatmap, width = 13, height = 14, dpi = 500, units =c("in"))

    ######## Pairwise correlation between samples (columns) ########
  #using mash_dist
  cols.cor.mash <- cor(mash_dist, use = "pairwise.complete.obs", method='pearson')
  #pairwise correlation between samples (rows)
  rows.cor.mash <- cor(t(mash_dist), use = 'pairwise.complete.obs', method='pearson')

    ######## Plot heatmap of pairwise correlations ########
  #using mash_dist
  mash_corr_heatmap <- pheatmap(mash_dist, scale = 'row',
                              clustering_distance_cols = as.dist(1 - cols.cor.mash),
                              clustering_distance_rows = as.dist(1 - rows.cor.mash),
                              main = "mash pairwise correlation heatmap")
  mash_corr_heatmap
  ggsave("FS19C_mash_correlation_heatmap.tiff", plot=mash_corr_heatmap, width = 11, height = 14, dpi = 500, units =c("in"))

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
  plot_mash_mdsA <-
  mash_mdsA %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=reference_id), max.overlaps = 50) +
  labs(x='MDS1', y='MDS2') +
  ggtitle("Mash MDS")
  plot_mash_mdsA
  ggsave("FS19C_mashMDS.tiff", plot=plot_mash_mdsA, width = 9, height = 8, dpi = 500, units =c("in"))
  ```

2. I googled how to determine outliers in MDS plot and came across this paper: https://doi.org/10.1093/bioinformatics/btz964. They use heatmaps of distance matrices to type bacteria, so I decided to try generating a heatmap to see how it compares between fastANI and mash. Followed this site for generating heatmap from distance matrix and correlation: https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
* This paper does something similar that ran mash and tested linear correlation with pearson's correlation coefficient test. Purpose was to look at whole genome similarity to identify bacterial meningitis causing species: https://doi.org/10.1186/s12879-018-3324-1

##### Files generated:
* fastANImashMDSheatmaps.pptx
* FS19C_fastaniMDS.tiff
* FS19C_mashMDS.tiff
* qc_mds.R

## Genome annotation with prokka, pan-genome analysis with roary, genomic island identification with gifrop
* Summary: The pan_pipe slurm script, provided by Jules, runs prokka, roary, and gifrop (developed by Julian Trachsel. Gifrop2 = gifrop version 2) via slurm on Ceres. It will annotate all fasta genomes (`.fna` format) with prokka in parallel (will do 24 genomes at a time, each with 1 thread), run roary and generate a core genome alignment, and with gifrop, it will extract, classify, and cluster genomic islands
* Github: https://github.com/Jtrachsel/gifrop
* Platform: Ceres


#### Files generated (for each isolate):
* Prokka
  * *_pol.fasta%.fasta_prokka/
    * *_pol.fasta%.fasta.err
    * *_pol.fasta%.fasta.faa
    * *_pol.fasta%.fasta.ffn
    * *_pol.fasta%.fasta.fna
    * *_pol.fasta%.fasta.fsa
    * *_pol.fasta%.fasta.gbk
    * *_pol.fasta%.fasta.gff
    * *_pol.fasta%.fasta.log
    * *_pol.fasta%.fasta.sqn
    * *_pol.fasta%.fasta.tbl
    * *_pol.fasta%.fasta.tsv
    * *_pol.fasta%.fasta.txt
* Roary
* gifrop
