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
* FastANI, multiqc

### mashenv
* Mash


## Sequence Analysis
### Sequence assembly
* Summary:
* Completed on:
* Platform:
* Commands:
1. Make text file replace.batch
```
#usr/bin/bash
while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
```

2. Run replace.batch
```
sh replace.batch
```

3. Rename *.fastq.gz files to fit the '*_1.fastq.gz' or '*_2.fastq.gz' format using the following bash script to create a directory "linked", take the fastq.gz files in current directory, shorten the names from something like 1-H12-96-441FEC_S103_L001_R2_001.fastq.gz to 1-H12-96-441FEC_2.fastq.gz, and copy over the short-named fastq.gz files to "linked" directory while also linking the short-named fastq.gz files to the original long-named fastq.gz files.
```
mkdir linked
for myfile in ./*.fastq.gz; do
        target=$(echo $myfile|sed -r 's/1\-[A-H][0-9]+\-([0-9]+\-[0-9]+[A-Z]+[0-9]?[A-Z]?\_)[A-Z]+[0-9]+\_[A-Z]+[0-9]+\_R([0-9])\_[0-9]+/\1\2/g')
        echo "$myfile is being renamed to $target"
        ln -sr $myfile linked/$target
done
```


4. Run the follow slurm script  
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
* fastqc.batch slurm script:

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
* Began on: 11Jan2021
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

* Files generated:
  * distances.tab
  * FS19Cmashdistances.xlsx


#### MDS
1. Make distance matrix from distances.tab and fs19cfastanioutput2.out
