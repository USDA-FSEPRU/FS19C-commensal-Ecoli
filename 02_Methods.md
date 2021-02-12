# Methods
Details of sequence analyses methods performed on FS19C samples 1-96 in this repo. Includes lab notes of how methods performed

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

## Conda environment (updated on 11Feb2021 - condensed to one conda environment)
### **Make sure when calling environments, to use this path: /project/fsepru/kmou/dot_files/.conda/envs/**
See: https://scinet.usda.gov/guide/conda/#user-installed-software-on-ceres-with-conda

### /project/fsepru/kmou/prokka_env
* fastani, multiqc, mash prokka, PPanGGOLiN
* created .condarc file in home directory to re-direct downloading packages in package cache from home directory to project directory
```
pkgs_dirs:
  - /project/fsepru/kmou/my_pkg_cache
channels:
  - r
  - conda-forge
  - bioconda
  - defaults
```

## shortcuts in .bashrc
alias debug='salloc -N 1 -p debug -t 01:00:00'
alias myjobs='squeue | grep kathy.mo'

## Sequence assembly
* Summary: QC and assemble sequences using BBMap and spades
* Began on: 16Sept2020
* Completed on: 18Dec2020
* Platform: Ceres

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

### Ran BBmap + spades assembly slurm script with the following genomes:
**2Nov2020**
1. Read ProkEvo

2. Made bash shell script to rename *fastq.gz filenames (shorten to wellID-SampleID_1/2.fastq.gz) and saved regex bash shell script as "renamefiles.batch"

3. Copy over samples.txt to Ceres project folder

4. Ran on SLURM the following genomes:
* 2-437RN1A.slurm: Submitted batch job 5262838
* 3-437REN7A.slurm: Submitted batch job 5262849
* 4-437REN7B.slurm: Submitted batch job 5262852
* 5-437REN3A.slurm: Submitted batch job 5262865
* 6-437REN3B.slurm: Submitted batch job 5262868
* 7-426FEN5A.slurm: Submitted batch job 5262869
* 8-426FEN5B.slurm: Submitted batch job 5262872
* 9-429FN2A.slurm: Submitted batch job 5262877
* 10-434FEN3.slurm: Submitted batch job 5262878

5. Checked the assemblies on 3Nov2020 and filled out `FS19C slurm progress.xlsx`

**10Nov2020**
1. Ran the next 9 (#11-19).
* 11-434FEN3.slurm: Submitted batch job 5275036
* 12-435FEN3.slurm: Submitted batch job 5275046
* 13-435FEN3.slurm: Submitted batch job 5275047
* 14-437FEN5.slurm: Submitted batch job 5275050
* 15-437FEN5.slurm: Submitted batch job 5275052
* 16-427REC.slurm: Submitted batch job 5275053
* 17-427RED.slurm: Submitted batch job 5275054
* 18-429REC.slurm: Submitted batch job 5275056
* 19-429RED.slurm: Submitted batch job 5275057

**12Nov2020**
1. Submitted the following jobs to slurm
* 21-429FEC.slurm: Submitted batch job 5281670
* 22-429FED.slurm: Submitted batch job 5281671
* 23-430FEC.slurm: Submitted batch job 5281672
* 24-430FED.slurm: Submitted batch job 5281679
* 25-427FED.slurm: Submitted batch job 5281681
* 26-433FEC.slurm: Submitted batch job 5281686
* 27-433FED.slurm: Submitted batch job 5281687
* 28-436FEC.slurm: Submitted batch job 5281688
* 29-436FED.slurm: Submitted batch job 5281689
* 30-440FEC.slurm: Submitted batch job 5281690

**18Nov2020**
1. Submitted the following jobs to slurm
* 31-440FED.slurm: Submitted batch job 5291498
* 32-436REC.slurm: Submitted batch job 5291535
* 33-436RED.slurm: Submitted batch job 5291539
* 34-440REC.slurm: Submitted batch job 5291543
* 35-440RED.slurm: Submitted batch job 5291559
* 36-441REC.slurm: Submitted batch job 5291560
* 37-441RED.slurm: Submitted batch job 5291562
* 38-426REC.slurm: Submitted batch job 5291563
* 39-426RED.slurm: Submitted batch job 5291565
* 40-427REC.slurm: Submitted batch job 5291566

**19Nov2020**
1. Previous jobs ran successfully. Submitted the following jobs to slurm
* 41-427RED.slurm: Submitted batch job 5292113
* 42-428REC.slurm: Submitted batch job 5292114
* 43-428RED.slurm: Submitted batch job 5292122
* 44-429REC.slurm: Submitted batch job 5292127
* 45-429RED.slurm: Submitted batch job 5292128
* 46-430REC.slurm: Submitted batch job 5292129
* 47-430RED.slurm: Submitted batch job 5292130
* 48-431RED.slurm: Submitted batch job 5292131
* 49-432REC.slurm: Submitted batch job 5292132
* 50-432RED.slurm: Submitted batch job 5292133

**24Nov2020**
1. Previous jobs ran successfully. Submitted the following jobs to slurm:
* 51-433REC.slurm: Submitted batch job 5297378
* 52-433RED.slurm: Submitted batch job 5297379
* 53-434REC.slurm: Submitted batch job 5297380
* 54-434RED.slurm: Submitted batch job 5297382
* 55-435REC.slurm: Submitted batch job 5297384
* 56-435RED.slurm: Submitted batch job 5297385
* 57-436REC.slurm: Submitted batch job 5297386
* 58-436RED.slurm: Submitted batch job 5297387
* 59-437REC.slurm: Submitted batch job 5297388
* 60-437RED.slurm: Submitted batch job 5297389

**30Nov2020**
1. Previous jobs ran successfully. Submitted the following jobs to slurm:
* 61-438REC.slurm: Submitted batch job 5300684
* 62-438RED.slurm: Submitted batch job 5300685
* 63-439REC.slurm: Submitted batch job 5300686
* 64-439RED.slurm: Submitted batch job 5300687
* 65-440REC.slurm: Submitted batch job 5300688
* 66-440RED.slurm: Submitted batch job 5300689
* 67-441REC.slurm: Submitted batch job 5300690
* 68-441RED.slurm: Submitted batch job 5300691
* 69-426FEC.slurm: Submitted batch job 5300692
* 70-426FED.slurm: Submitted batch job 5300693

**4Dec2020**
1. Previous jobs ran successfully. Submitted the following jobs to slurm:
* 71-427FED.slurm: Submitted batch job 5307961
* 72-428FEC.slurm: Submitted batch job 5307963
* 73-428FED.slurm: Submitted batch job 5307964
* 74-429FEC.slurm: Submitted batch job 5307965
* 75-429FED.slurm: Submitted batch job 5307966
* 76-430FEC.slurm: Submitted batch job 5307978
* 77-430FED.slurm: Submitted batch job 5307979
* 78-431FEC.slurm: Submitted batch job 5307980
* 79-431FED.slurm: Submitted batch job 5307981
* 80-432FEC.slurm: Submitted batch job 5307982
* 81-432FED.slurm: Submitted batch job 5307984
* 82-433FEC.slurm: Submitted batch job 5308053
* 83-433FED.slurm: Submitted batch job 5308054
* 84-434FEC.slurm: Submitted batch job 5308056
* 85-434FED.slurm: Submitted batch job 5308057
* 86-435FEC.slurm: Submitted batch job 5308058
* 87-435FED.slurm: Submitted batch job 5308059
* 88-436FEC.slurm: Submitted batch job 5308060
* 89-436FED.slurm: Submitted batch job 5308061
* 90-437FED.slurm: Submitted batch job 5308062
* 91-438FEC.slurm: Submitted batch job 5308063
* 92-438FED.slurm: Submitted batch job 5308064
* 93-439FEC.slurm: Submitted batch job 5308065
* 95-440FED.slurm: Submitted batch job 5308067

2. 95 did not assemble as expected because Mike Baker said for samples 95 and 96, they could not get any sequences.

3. In total, 20, 95, and 96 did not assemble. Ignore #20? Wait for 95 and 96.

4. Finished assembling! Will need to transfer polished fasta sequences (#11-93) to `polished_genomes_100X directory`. Is it possible to move files from (11-93)_pol.fasta to directory? Learned to use grep.

**18Nov2020** Received email from Darrell of re-sequenced isolates 95-96 fasta files.
1. Darrell email:
```
Your NovaSeq data has been downloaded from the ISU Sequencing Center and has been archived here at the Center.  According to David, there were a couple of the original samples, 95 and 96 (H11 and H12) in the original plate, which did not provide data when completed.  This new data is the result after Dr. Baker has repeated the library preparation and rerun those samples.  
 
I’ve placed the data on the Q: drive at: Q:\_TempTransfer\DBayles\mou.  Additionally, I put a copy on Ceres where you can access it at: /90daydata/shared/mou_201216/."
```

2. I initially didn't see the two files and asked Darrell. He said the following:
```
There is one file, but when unpacked contains multiple sample FASTQ files for both the repeats.  That’s just how the Sequencing Center bundles them.  I know it can be confusing when ISU uses the name of a single sample for their archive name and then includes multiple other samples inside that archive.
 
When you unpack the archive, you will find the following files contain within:
1-H11-95-440FED_S1_L002_R1_001.fastq.gz
1-H11-95-440FED_S1_L002_R2_001.fastq.gz
1-H12-96-441FEC_S2_L002_R1_001.fastq.gz
1-H12-96-441FEC_S2_L002_R2_001.fastq.gz
 
You should have all the replacement data."
```

3. Yep! I un-tarred the files, transferred them to `/project/fsepru/kmou/FS19C/lane1/`
```
1_33298_01_1-A01-1-428RN3A_HVHJT_1839.tar
1-H11-95-440FED_S1_L002_R1_001.fastq.gz
1-H11-95-440FED_S1_L002_R2_001.fastq.gz
1-H12-96-441FEC_S2_L002_R1_001.fastq.gz
1-H12-96-441FEC_S2_L002_R2_001.fastq.gz
```

4. Rename old 96 fastq.gz files with "_old" so I don't confuse the re-sequenced files with the old sequence files.
```
mv 96-441FEC_1.fastq.gz 96-441FEC_1_old.fastq.gz
mv 96-441FEC_2.fastq.gz 96-441FEC_2_old.fastq.gz
```

5. Copy 95 and 96 fastq.gz files to `linked` directory. Rename them to
```
95-440FED_1.fastq.gz
95-440FED_2.fastq.gz
96-441FEC_1.fastq.gz
96-441FEC_2.fastq.gz
```

6. Made new slurm script for #96 because it didn't exist in linked/ directory.

7. Submitted the following jobs to slurm:
* 95-440FED.slurm: Submitted batch job 5361148
* 96-441FEC.slurm: Submitted batch job 5361152

8. Jobs completed with no issues. Copied polished fasta genomes to `polished_genomes` directory from `linked` directory with `mv *pol.fasta` command

9. Ran this to get complete list of polished genomes
```
ls -v *pol.fasta > polishedfasta.txt
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
* Summary: Ran fastqc on FS19C sequence data to assess sequence quality of individual reads for each sample, and to use output for multiqc (forgot to do this prior to sequence assembly)
* Began on: 6Jan2021
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

##### Files generated (for each isolate):
  * *fastqc.zip
  * *fastqc.html

## QC with MultiQC
* Tutorial: https://www.youtube.com/watch?v=qPbIlO_KWN0
* Summary: Ran multiqc with fastqc output of FS19C sequence data to assess quality of all sequences for samples 1-94 (forgot to do this prior to sequence assembly). I did not include samples 95 and 96 in the multiqc run as I realized these samples were ran on a second sequencing run, so their coverage is different than the first sequencing run that had all samples (but I would only consider  samples 1-94 from that first run). In addition, examining fastqc output for samples 95 and 96 was fine enough.
* Began and Completed on: 7Jan2021
* Platform: Ceres, fastanienv conda environment

1. Command ran:
  ```
  salloc
  module load miniconda
  source activate fastanienv
  conda install -c bioconda multiqc
  multiqc *.fastqc.zip
  conda deactivate
  ```

2. Why are the plots flat plot (not interactive)?
  * From: https://multiqc.info/docs/
  "Flat plots
  Reports with large numbers of samples may contain flat plots. These are rendered when the MultiQC report is generated using MatPlotLib and are non-interactive (flat) images within the report. The reason for generating these is that large sample numbers can make MultiQC reports very data-intensive and unresponsive (crashing people's browsers in extreme cases). Plotting data in flat images is scalable to any number of samples, however.
  Flat plots in MultiQC have been designed to look as similar to their interactive versions as possible. They are also copied to multiqc_data/multiqc_plots"

3. I noticed the FS19all_multiqc_report.html report showed samples 95 and 96 having huge number of reads, per sequence GC content had several large peaks, and samples 95 and 96 making up the majority of the overrepresented sequences.

4. To try to eliminate the discrepancies due to samples 95 and 96 (these two samples had so many reads as a result of re-sequencing them) and are therefore skewing the report stats, I moved samples 95 and 96 fastqc.zip to Fastqc_Sample95_96 directory, and ran multiqc to generate FS19_1-94_multiqc_report.html and FS19_1-94_data directory.

5. multiqc report of samples 1-94 look a lot better with the sequence count ranges being a lot closer among all samples, sequence quality histograms all in green, per sequence quality scores in green, less than 1% of reads making up overrepresented sequences, and a single peak for per sequence GC content

6. I also looked at the fastqc reports for samples 95 and 96 individually.
  * **95** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9
  * **96** (both reads): quality scores are green through entire position, some sequence duplication levels starting at 9, peak at >10 and ends at >500; per base sequence content is iffy from positions 1-9

7. Note: Jules says if you're able to get assemblies to work, you can be sure that the sequences' qualities are good

##### Files generated:
 * FS19all_multiqc_report.html
 * FS19all_multiqc_data directory
 * FS19_1-94_multiqc_report.html
 * FS19_1-94_multiqc_data directory
 * 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
 * 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R1_001_fastqc.html

## QC with FastANI
* Summary: ran fastANI on FS19C sequence data by running in conda environment to estimate Average Nucleotide Identity (ANI) using alignment-free approximate sequence mapping. It calculates distance between 2 sequences. Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* FastANI publication: DOI: 10.1038/s41467-018-07641-9
* Began on: 28Dec2020
* Completed on: 11Jan2021
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

2. (28Dec2020) Testing fastani on sample genomes provided by https://github.com/ParBLiSS/FastANI
```
fastANI -q Shigella_flexneri_2a_01.fna -r Escherichia_coli_str_K12_MG1655.fna -o testfastani.out
```
Viewed output testfastani.out and it looked the same as output from https://github.com/ParBLiSS/FastANI. Awesome!

3. (28Dec2020) Made query list in polished_genomes/ with this command and checked to make sure query list only contains sample names (tells the path to where reference genome files are):
```
ls -dv "$PWD"/* > quertylist.txt
```

4. (28Dec2020) Made reference list with this (tells the path to where reference genome files are):
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

5. (28Dec2020) Moved querylist.txt and referencelist.txt to same directory and ran fastani:
```
fastANI --ql querylist.txt --rl referencelist.txt -o fs19cfastanioutput.out
```

6. (28Dec2020) Downloaded fs19cfastanioutput.out and made FS19CfastANIoutput.xlsx

7. (29Dec2020, 11Jan2021) After discussion with Jules about how I'd run mash, I realized I ran fastani incorrectly. I should be measuring distances of all genomes (references and samples) between each other. So combine all file names in a new querylist2.txt and run fastani with this file listed as query and as reference lists.
```
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
* Began and Completed on: 11Jan2021
* Platform: Ceres, mashenv conda environment

1. Load environment and install mash
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

Output generated: .msh

3. Run distance estimates comparing sketch genomes with itself. This runs really fast (less than a second)
```
mash dist -p 1 .msh .msh > distances.tab
```

4. distances.tab columns: reference_id, query_id, mash_distance, pvalue, matching_hashes

5. (10Feb2021) Run Mash with all 95 isolates, and 5 reference genomes including TW14588 on slurm with `mash.slurm` script
```
#!/bin/bash
#SBATCH --job-name=mash                           # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
set -e
module load miniconda
source activate /project/fsepru/kmou/dot_files/.conda/envs/mashenv
mash sketch -p 1 -o /project/fsepru/kmou/FS19C/polished_genomes_100X/mash_all/ *.fasta
mash dist -p 1 .msh .msh > distances_secondrun.tab
```
```
Submitted batch job 5533997
```

6. Downloaded `distances_secondrun.tab` to local computer, moved file to `Files/` and import to `mash_mds.R` script and created MDS. The plot is pretty much the same as previous mash results with TW14588 buried within the large cluster. Not sure why roary is giving so many errors.



##### Files generated:
  * distances.tab
  * FS19Cmashdistances.xlsx

## QC: Visualize ANI pairwise genome-genome similarity calculations with MDS, heatmap
* Summary: Made distance matrix from mash and fastani output to create heatmap and MDS to visualize clustering and identify any outliers. The MDS was a bit hard to decipher what was an outlier, so I ran a heatmap to see how fastANI and mash compared and whether the pairwise comparisons were similar between the two, including heatmap of pearson correlation coefficients.
* Began on: 14Jan2021
* Completed on: 15Jan2021
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

10. (22Jan2021) At microbe meeting, Crystal pointed out that it's good we're seeing some differences of the commensal isolates (like isolates collected from EDL933 group cluster separately from EDL933 isolate) from the STEC isolates in this plot. Regardless of whether we do find any differences in metabolic genes or not, there are other differences we can explore too.

11. Next step: run prokka. I can use UnitProt pangenome E. coli to run annotation. I asked if I chould do this or run a single reference genome for annotation and Jules said I should try the pangenome method first. Will need to find the E. coli pangenome annotation (protein).

##### Files generated:
* fastANImashMDSheatmaps.pptx
* FS19C_fastaniMDS.tiff
* FS19C_mashMDS.tiff
* qc_mds.R

## Genome annotation with prokka
* Summary: Identify annotation reference from UniProt (E. coli pangenome) and use prokka to annotate all 95 samples.
* Prokka publication: DOI: 10.1093/bioinformatics/btu153
* Began on: 20Jan2021
* Completed on: 28Jan2021
* Platform: Ceres, prokka_env conda environment

1. (20Jan2021) Find E. coli pangenome (pan proteome) on UniProt:
  * [Escherichia coli (strain K12) (Strain: K12 / MG1655 / ATCC 47076)](https://www.uniprot.org/proteomes/UP000000625)
  * [Escherichia coli O157:H7 (Strain: O157:H7 / Sakai / RIMD 0509952 / EHEC)](https://www.uniprot.org/proteomes/UP000000558)

2. (20Jan2021) Read Prokka paper
  * DOI: 10.1093/bioinformatics/btu153

3. (20Jan2021) Find prokka github page: https://github.com/tseemann/prokka

4. (20Jan2021) install prokka in fastanienv conda environment on Ceres in FS19C/polished_genomes_100X
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

6. (21Jan2021) Made a new directory polishedgenomesprokka/, copied fasta files from polished_genomes_100X to this directory, renamed .fasta to .fna
```
rename .fasta .fna *.fasta
```

7. (21Jan2021) Download [E. coli strain K12 proteome (fasta)](https://www.uniprot.org/proteomes/UP000000625) - Try this annotation first.

8. (21Jan2021) Need to ask Jules what kind of fasta file do I need to use for annotation first. The fasta file I have doesn't have the ~~~ symbols that prokka says is needed for annotation tag formats: https://github.com/tseemann/prokka/blob/master/README.md#fasta-database-format

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

11. (27Jan2021) Looked through each EMBL file and print out all proteins associated with each of the substrates:
  * Escherichia coli O157:H7 str. Sakai
    ```
    grep -f carbohydratelist.txt BA000007.3.txt > completecarblistO157H7.txt
    ```
  * Escherichia_coli_str_K12_MG1655
    ```
    grep -f carbohydratelist.txt U00096.3.txt > completecarblistK12.txt
    ```

12. (27Jan2021) Compared the two files and printed out genes unique to each strains
  ```
  comm -1 -3 completecarblistO157H7.txt completecarblistK12.txt > genesnotpresentinO157H7.txt #not present in O157H7
  comm -1 -3 completecarblistK12.txt completecarblistO157H7.txt > genesnotpresentinK12.txt
  ```

13. (27Jan2021) So, at least one of the strains' proteomes have genes associated with each of the metabolites listed in OSQR plan. I will also go through review papers referenced in OSQR plan to see what other genes to include.
 * Maltby et al.: EMBL files missing glucuronate, galacuronate; however, based on a paper referenced by Maltby et al., glucoronate, galacuronate, and hexuronate are one sugar

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

15. (28Jan2021) Looked up issues page on prokka: https://github.com/tseemann/prokka/issues/86 and checked that my files are fasta format, they are readable, and are more than 0 bytes. Emailed Jules how to fix error. He says it looks like my prokka command has some weird formatting applied and some special characters were inserted where they shouldn't be. This is a danger when copying things from the internet. There's a difference between long and short dash. Shell scripts don't take too kindly to things like long dashes and it confuses them, along with other hidden special characters. So I retyped the command on bbedit, and prokka ran past the part I got stuck on. However, it didn't like the length of contig IDs in my fasta files, so it suggested renaming contigs or try '--centre X --compliant' to generate clean contig names. Need to shorten them to less than or equal to 37 characters long.
```
[10:20:41] Contig ID must <= 37 chars long: NODE_1_length_378381_cov_16.779345_pilon_pilon_pilon
[10:20:41] Please rename your contigs OR try '--centre X --compliant' to generate clean contig names.
```
I added --centre X --compliant command line options
```
for file in *.fasta; do tag=$file%.fasta; prokka -prefix "$tag" -locustag "$tag" -genus Escherichia -strain "$tag" -outdir "$tag"_prokka -force -addgenes "$file" -centre X -compliant; done
```

16. (28Jan2021) It is working! Started around 10:30AM, finished at 11pm. All samples have a new directory with .err, .faa, .fnn, .fna, .fsa, .gbk, .gff, .log, .sqn, .tbl, .tsv, .txt files

17. (8Feb2021) Need to re-run prokka to get same annotation for 6 E. coli reference fasta files and 95 isolates. Results from this prokka run will be used for ppangolin (use gbk files). Waiting to hear back from Jules on how to run prokka to get correct annotation and output for roary. Ran prokka.slurm (ask slurm to load miniconda and prokka_env).
```
#!/bin/bash
#SBATCH --job-name=prokka                            # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
set -e #this causes the script to exit on error: https://marcc-hpc.github.io/tutorials/shortcourse_slurm.html
module load miniconda
source activate ~/dot_files/.conda/envs/prokka_env
for file in *.fasta; do tag=$file%.fasta; prokka -prefix "$tag" -locustag "$tag" -genus Escherichia -strain "$tag" -outdir "$tag"_prokka -force -addgenes "$file" -centre X -compliant; done
```
```
Submitted batch job 5531180
```

18. (9Feb2021) Prokka didn't like how long contig IDs were (over 37) for reference genome MG1655 and didn't continue on to do Ecoli_NADC6564, Ecoli_Nissle1917, Ecoli_O157H7_EDL933, and Ecoli_TW14588. So I'm shortening both file names and the header line in fasta files to see if that makes a difference. Will run `debug` to test if prokka will accept those files before running with all other 95 isolates.
```
EDL933.fasta
>CP008957.1 EDL933, complete genome # header of EDL933.fasta
TW14588.fasta
>NZ_CM000662.1 TW14588, WGS shotgun # header of TW14588.fasta
MG1655.fasta
>NC_000913.3 MG1655, complete genome # header of MG1655.fasta
NADC6564.fasta
>CP017251.1 NADC6564, complete sequence # header of NADC6564.fasta
Nissle1917.fasta
>CP007799.1 Nissle 1917, complete genome # header of Nissle1917.fasta
```

19. (9Feb2021) Moved fasta files to test directory. Run debug mode and test following code. Seemed to work fine, although I can't run it long enough to see if the same contig ID error message will show up. But each of the reference strains got a prokka_fasta folder made.
```
module load miniconda
source activate ~/dot_files/.conda/envs/prokka_env
for file in *.fasta; do tag=$file%.fasta; prokka -prefix "$tag" -locustag "$tag" -genus Escherichia -strain "$tag" -outdir "$tag"_prokka -force -addgenes "$file" -centre X -compliant; done
```

20. (9Feb2021) Submitted slurm job using `prokka.slurm`. Job #5532526, started at 11:09am and finished at 3:25pm. Stderr and stdout logs say job was completed with no errors.

#### Files generated (for each isolate):
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

## Pangenome analysis with roary
* Summary: ran Roary on FS19C gff data by running in Ceres to generate pangenome analysis of *E. coli* isolates.
* Roary publication DOI: 10.1093/bioinformatics/btv421
* Began on: 29Jan2021
* Completed on:
* Platform: Ceres

1. (29Jan2021) Notes from Roary publication (including supplemental info)
  * pass in the flag '-e' to get multi-fasta file to use with RAxML or FastTree to generate phylogenetic tree based on SNPs in core genes. The file you want for those applications is called *core_gene_alignment.aln*. This flag also generates *pan_genome_reference.fa* file
  * Access this tutorial for step-by-step: https://github.com/microgenomics/tutorials/blob/master/pangenome.md

2. (29Jan2021) Roary is on Ceres. I will copy all gff from each folder and make a new folder of only gff files. Then I can roary with them with some kind of loop on slurm.
```
find . -name *.gff -exec cp '{}' "./prokka_gff/" ";"
```
Adapted from this [forum](https://unix.stackexchange.com/questions/67503/move-all-files-with-a-certain-extension-from-multiple-subdirectories-into-one-di).

3. (29Jan2021) Ran md5sum checksum to compare isolate 96's gff file from prokka_gff/ with 96's gff file from 96-441FEC_pol.fasta%.fasta_prokka/
```
md5sum 96-441FEC_pol.fasta%.fasta.gff ../prokka_gff/96-441FEC_pol.fasta%.fasta.gff > cksum96.txt
md5sum -c cksum96.txt
# 96-441FEC_pol.fasta%.fasta.gff: OK
# ../prokka_gff/96-441FEC_pol.fasta%.fasta.gff: OK
```

4. (29Jan2021) The find command worked, so I will generate slurm script (roary.slurm) to run roary on the 95 isolates.
```
#!/bin/bash
#SBATCH --job-name=roary                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 12:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
module load roary
roary -f ./roary_output -e -n -v *.gff
#End of file
```

5. (29Jan2021) Submitted job on SLURM
```
sbatch roary.slurm
Submitted batch job 5495686
```

6. (1Feb2021) Job didn't finish because my job reached the time limit of the partition: short partition for 48 hours. Job got cancelled. I will need to re-run roary using a partition with longer simulation time (medium, 7 days). I emailed Jules to ask about how to know how much of each resource from Ceres to use for a job.

7. (1Feb2021) Jules asked what command I used with roary, so I sent him my slurm script. After looking through, he says my initial resource request should be more than enough. The problem is I'm not telling roary to use more than one cpu/core/thread. Check out roary documentation to see how to tell roary to use more than one processor/thread. I looked up and found an option to include:
```
roary -f ./roary_output -e -n -v -p 16 *.gff
```

8. (1Feb2021) Modified roary.slurm script to the following:
```
#!/bin/bash
#SBATCH --job-name=roary                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                   # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                     # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
module load roary
roary -f ./roary_output -e -n -v -p 16 *.gff
#End of file
```

9. (1Feb2021) Submitted job on SLURM
```
sbatch roary.slurm
Submitted batch job 5507099
```

10. (1Feb2021) Job ran from 10:24am to 1:10pm and completed! Copied `roary_output` directory to local desktop
```
scp -r ceres:~/fsepru_kmou/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gff/roary_output ./
```

11. (2Feb2021) See gene differences between groups of isolates using `query_pan_genome` command. Copy all files from `roary_output` directory to same directory as gff files because `query_pan_genome` uses some of those output files (as I found out the first time I ran the slurm script below. The job cancelled because it couldn't find `clustered_proteins` file).
```
#!/bin/bash
#SBATCH --job-name=roary                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
module load roary
#roary -f ./roary_output -e -n -v -p 16 *.gff
query_pan_genome  -o pan_genome_results_union -v -a union *.gff
query_pan_genome  -o pan_genome_results_core -v -a intersection *.gff
query_pan_genome  -o pan_genome_results_accessory -v -a complement *.gff
#End of file
```
```
Submitted batch job 5511086
```

12. (2Feb2021) Job completed, downloaded `pan_genome_results_core`, `pan_genome_results_accessory`, and `pan_genome_results_union` to local `roary_output/querypangenome_output` directory.

13. (3Feb2021) Looked at the following roary output files:
* `gene_presence_absence.csv`
  * Order within fragment, combined with the genome fragment this gives an indication of the order of genes within the graph. In Excel, sort on Column G and H.
  * Accessory order with fragment, combined with the Accessory fragment this gives an indication of the order of genes within the accessory graph. In Excel, sort on columns I and J.
  * QC column: many that have "investigate" or "hypothetical"
  * How to group gene names by pathways?
* `gene_presence_absence.Rtab`
  * simple tab delimited binary matrix with presence and absence of each gene in each sample.
  * The first row is the header containing the name of each sample, and the first column contains the gene name. A 1 indicates the gene is present in the sample, a 0 indicates it is absent.
  * read Rtab file in `roary.R`

14. (3Feb2021) Run `create_pan_genome_plots.R` from [Roary github](https://github.com/sanger-pathogens/Roary/blob/master/bin/). Saved as `roary.R`. Generated many plots, below are comments of each plot made.
* `number_of_new_genes.Rtab`: 95 genomes show very few new genes
* `number_of_conserved_genes.Rtab`: 95 genomes suggest a little over 3K conserved genes
* `number_of_genes_in_pan_genome.Rtab`: 95 genomes show approaching 14000 genes in pan-genome
* `number_of_unique_genes.Rtab`: 95 genomes show between 2000-3000 unique genes
* `blast_identity_frequency.Rtab`: 95 genomes show below 5000 blastp results with 98%, 99% blastp identity; between 5000-10000 blastp results with 95-97% blastp identity, and more than 15000 blastp results with 100% blastp identity.
* `conserved_vs_total_genes.png` graph showing number of conserved genes to number of total genes in pan-genome (no surprise, more total genes than conserved genes)
* `unique_vs_new_genes.png` graph showing number of new genes vs unique genes (more unique than new genes)

15. (3Feb2021) Looking through FS19C plan, I realized I forgot to add other reference genomes gff files to roary analysis as landmarks for comparing with FS19C isolates. I downloaded GFF files from the following isolates via accessing their Genbank file (see hyperlink above in fastANI section), click on "Assembly" under "Related Information" on the right-side panel, then click on the only organism listed on results page, then click on "Download Assembly" on upper right side of page, select RefSeq (Source database) and Genomic GFF (File type), and download:
* E. coli MG1655 https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2
  * Saved as Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gff.gz
* E. coli HS https://www.ncbi.nlm.nih.gov/assembly/GCF_000017765.1
  * Saved as Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gff.gz
* E. coli Nissle 1917 https://www.ncbi.nlm.nih.gov/assembly/GCF_000714595.1
  * Saved as Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gff.gz
* E. coli O157:H7 str. NADC 6564 https://www.ncbi.nlm.nih.gov/assembly/GCF_001806285.1
  * Saved as Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gff.gz
* E. coli O157:H7 EDL933 https://www.ncbi.nlm.nih.gov/assembly/GCF_000732965.1
  * Saved as EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gff.gz
* TW14588 https://www.ncbi.nlm.nih.gov/assembly/GCF_000155125.1/
  * Saved as Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gff.gz
* Notes:
  * What is the difference between GenBank and RefSeq genome assembly? https://www.ncbi.nlm.nih.gov/datasets/docs/gca-and-gcf-explained/
  * [RefSeq FAQS](https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.what_is_a_reference_sequence_r)

16. (3Feb2021) No annotated assemblies available for these two E. coli isolates used in FS19 studies:
* RM6067 - no assembly publicly available, also no annotation found (only found RNAseq-related info about this strain)
* FRIK1989 - no assembly or annotation publicly available

17. (3Feb2021) Uploaded 6 E. coli reference genomes' gff files to `project/fsepru/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gff/EcoliReferenceGenomes` directory on Ceres. Unzipped files with `gzip -d *.gz` and copied to `prokka_gff` directory. Ran slurm job of roary on gff files of 95 isolates + 6 reference E. coli using `roary.slurm` script. Job completed in 3h49m.
```
module load roary
roary -f ./roary_95isolates_6referencestrains_output -e -n -v -p 16 *.gff
Submitted batch job 5518586
```

18. (4Feb2021) Copy files from `roary_95isolates_6referencestrains_output` to `prokka_gff` and run `query_pan_genome` with `roary.slurm` script.
```
module load roary
#roary -f ./roary_95isolates_6referencestrains_output -e -n -v -p 16 *.gff
query_pan_genome  -o pan_genome_results_union -v -a union *.gff
query_pan_genome  -o pan_genome_results_core -v -a intersection *.gff
query_pan_genome  -o pan_genome_results_accessory -v -a complement *.gff
```
```
Submitted batch job 5521850
#stderr.5521850.ceres14-compute-33.roary
2021/02/04 12:26:43 Could not extract any protein sequences from Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:26:45 Could not extract any protein sequences from Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:26:47 Could not extract any protein sequences from Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:26:48 Could not extract any protein sequences from EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:26:50 Could not extract any protein sequences from Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:26:52 Could not extract any protein sequences from Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:28 Could not extract any protein sequences from Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:29 Could not extract any protein sequences from Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:31 Could not extract any protein sequences from Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:33 Could not extract any protein sequences from EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:34 Could not extract any protein sequences from Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:33:36 Could not extract any protein sequences from Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:11 Could not extract any protein sequences from Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:12 Could not extract any protein sequences from Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:14 Could not extract any protein sequences from Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:15 Could not extract any protein sequences from EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:17 Could not extract any protein sequences from Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gff. Does the file contain the assembly as well as the annotation?
2021/02/04 12:40:19 Could not extract any protein sequences from Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gff. Does the file contain the assembly as well as the annotation?
```

19. (4Feb2021) Looked up error message and found this issue post: https://github.com/sanger-pathogens/Roary/issues/417. Tested with the gff and fna files of Mycoplasma pneumoniae FH that was mentioned in post: https://www.ncbi.nlm.nih.gov/assembly/GCF_001272835.1 (see script below). Seemed to work? I will download the fna files and cat the gff and fna files together to make gff3 and run them with roary.
```
cat GCF_001272835.1_ASM127283v1_genomic.gff GCF_001272835.1_ASM127283v1_genomic.fna > GCF_001272835.1_ASM127283v1_genomic.gff3
```

20. (5Feb2021) Need to annotate reference strain genomes with isolate genomes - get the same annotation. Asked Jules for specifics.
(9Feb2021) Response from Jules
```
I mentioned using genbank files during the prokka annotation step.  This is optional and only necessary if you want the genes to be annotated in a similar way to an established reference genome.

For example, when I annotate Salmonella genomes I will often call prokka using something like “--proteins LT2.gbk”.  This will make prokka try and annotate the experimental genomes using the annotations from the LT2 genome as a first priority.  I don’t believe this will change what genes are identified, just what they are called.  Again, this is totally optional.
```
* **I will ask E. coli group if they want me to annotate with a specific reference genome. For now, I will keep going and analyze results as is.**

21. (5Feb2021) Moved dot files from home directory on ceres to /project/fsepru/kmou/dot_files and made symbolic link to home
```
ln -s /project/fsepru/kmou/dot_files/software/ . # this is for software folder
ln -s /project/fsepru/kmou/dot_files/ . # this is for conda
```

22. (8Feb2021) Remove prokka_gbk/ and prokka_gff/, polishedgenomesprokka/.
Check FNA files (same as fasta files on Ceres?). Also uploaded Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.fna.gz fasta file to Ceres and unzipped with `gzip -d Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.fna.gz`

* FNA files:
* E. coli MG1655: already uploaded to Ceres (Ecoli_K-12_MG1655.fasta)
  * Ran md5sum checksum to compare fna file (downloaded from assembly page: RefSeq) with fasta file. They're both the same.
```
md5sum 96-441FEC_pol.fasta%.fasta.gff ../prokka_gff/96-441FEC_pol.fasta%.fasta.gff > cksum96.txt
md5sum -c cksum96.txt
# Ecoli_K-12_MG1655.fasta: OK
#Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.fna: OK
```
* E. coli HS: already uploaded to Ceres (Ecoli_HS.fasta)
* E. coli Nissle 1917: already uploaded to Ceres (Ecoli_Nissle1917.fasta)
* E. coli O157:H7 str. NADC 6564: already uploaded to Ceres (Ecoli_NADC6564.fasta)
* E. coli O157:H7 EDL933: already uploaded to Ceres (Ecoli_O157H7_EDL933.fasta)
* TW14588 https://www.ncbi.nlm.nih.gov/assembly/GCF_000155125.1/
  * Saved as Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.fna.gz
  * Renamed this file to Ecoli_TW14588.fasta on Ceres

23. (9Feb2021) Copied gff files to `prokka_gff` and place in `project/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gff` directory
```
find ../polishedgenomesprokka_95isolates5refgenomes/ -name *.gff -exec cp '{}' "./" ";"
```

24. (9Feb2021) Ran roary.slurm script to the following:
```
#!/bin/bash
#SBATCH --job-name=roary                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                   # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                     # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
module load roary
roary -f ./roary_95isolates_5referencestrains_output -e -n -v -p 16 *.gff
#End of file
```
```
Submitted batch job 5533863
```

25. (10Feb2021) Job 5533863 had more errors! I forgot to put `set -e` in slurm script, so it kept going. Some example error messages listed below. I looked at a few gff files and they have the annotation and the sequence data. I looked up the first error message (sequence without letters) and found this forum post: https://www.biostars.org/p/365927/. Also this issue page on roary: https://github.com/sanger-pathogens/Roary/issues/229. Could be that I have some reference genomes that are unrelated to my 95 isolates. Maybe it's TW14588? I'll try running mash to see what the distances are when I add TW14588 to the mix. In the meantime, I will also run `roary.slurm` without TW14588 (deleted `TW14588.fasta%.fasta.gff` from directory) and see how that goes.
```
#line 5:
2021/02/09 23:34:18 Extracting proteins from GFF files
Warning: unable to close filehandle $bed_fh properly: Disk quota exceeded at /usr/share/perl5/Bio/Roary/BedFromGFFRole.pm line 41.
2021/02/10 11:43:18 Could not extract any protein sequences from /project/fsepru/kmou/FS19C/polished_genomes_100X/prokka_gff/10-434FEN3_pol.fasta%.fasta.gff.

2021/02/09 23:34:21 Could not extract any protein sequences from /project/fsepru/kmou/FS19C/polished_genomes_100X/prokka_gff/10-434FEN3_pol.fasta%.fasta.gff. Does the file contain the assembly as well as the annotation?
--------------------- WARNING ---------------------
MSG: Got a sequence without letters. Could not guess alphabet
---------------------------------------------------
```
```
Submitted batch job 5533972
```

26. (10Feb2021) Same error messages popped up for job 5533972. Will test run roary on just the 95 isolates. I was able to get it to work the first time I ran roary with just the 95 isolates. Moved reference strain fasta files to `refgenomes/`. Run `roary.slurm`
```
set -e
set -u
module load roary
roary -f ./roary_95isolates_testrun_output -e -n -v -p 16 *.gff
```
```
Submitted batch job 5534101
```

27. (10Feb2021) Job 5534101 did not complete. Same error messages popped up. Something related to disk quota exceeding? I tried to find `/usr/share/perl5/Bio/Roary/BedFromGFFRole.pm` but could not find it on my end. Emailed vsrc support for help.

28. (11Feb2021) Re-installed conda environment `prokka_env` in correct directory (`/project/fsepru/kmou/`) following SciNet best practices: https://scinet.usda.gov/guide/conda/#user-installed-software-on-ceres-with-conda. This may help with disk quota in home directory? Will try running `roary.slurm` again and see if that changes.
```
#!/bin/bash
#SBATCH --job-name=roary                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#SBATCH --account fsepru
#Enter commands here:
set -e
set -u
module load roary
roary -f ./roary_95isolates5genomes_testrun_output -e -n -v -p 16 *.gff
```
```
Submitted batch job 5545690
```
Looked at stderr and stdout and no error messages were found. It is strange that stdout is so short, but the output files are all there. `summary_statistics.txt` looks good too:
```
Core genes	(99% <= strains <= 100%)	3135
Soft core genes	(95% <= strains < 99%)	329
Shell genes	(15% <= strains < 95%)	3413
Cloud genes	(0% <= strains < 15%)	9383
Total genes	(0% <= strains <= 100%)	16260
```

29. Renamed `roary_95isolates5genomes_testrun_output` to `roary_95isolates5refgenomes_final`. I also deleted output directories from previous failed roary runs (`roary_95isolates_5referencestrains_output/`, `roary_95isolates_5referencestrains_output_BAD/`)

30. Run the the `query_pan_genome` commands on slurm. Made sure to make soft link for `clustered_proteins` in same location as *.gff files.
```
module load roary
query_pan_genome  -o pan_genome_results_union -v -a union *.gff
query_pan_genome  -o pan_genome_results_core -v -a intersection *.gff
query_pan_genome  -o pan_genome_results_accessory -v -a complement *.gff
```
```
Submitted batch job 5550887
```

31. Download `roary_95isolates5refgenomes_final/` files and analyze.
Looked at the following roary output files:
* `gene_presence_absence.csv`
  * Order within fragment, combined with the genome fragment this gives an indication of the order of genes within the graph. In Excel, sort on Column G and H.
  * Accessory order with fragment, combined with the Accessory fragment this gives an indication of the order of genes within the accessory graph. In Excel, sort on columns I and J.
  * QC column: "investigate", "hypothetical" or blank
  * How to group gene names by pathways?
* `gene_presence_absence.Rtab`
  * simple tab delimited binary matrix with presence and absence of each gene in each sample.
  * The first row is the header containing the name of each sample, and the first column contains the gene name. A 1 indicates the gene is present in the sample, a 0 indicates it is absent.
  * read Rtab file in `roary.R`
* `pan_genome_results_union`: Union of genes found in isolates
* `pan_genome_results_core`: Intersection of genes found in isolates (core genes)
* `pan_genome_results_accessory`: Complement of genes found in isolates (accessory genes)

14. (3Feb2021) Run `create_pan_genome_plots.R` from [Roary github](https://github.com/sanger-pathogens/Roary/blob/master/bin/). Saved as `roary.R`. Generated many plots, below are comments of each plot made.
* `number_of_new_genes.Rtab`: 95 genomes show very few new genes
* `number_of_conserved_genes.Rtab`: 95 genomes suggest a little over 3K conserved genes
* `number_of_genes_in_pan_genome.Rtab`: 95 genomes show approaching 14000 genes in pan-genome
* `number_of_unique_genes.Rtab`: 95 genomes show between 3000-4000 unique genes
* `blast_identity_frequency.Rtab`: The 95 genomes show that between 5000-10,000 blastp results have 95-97% blast identity, less than 5000 blastp hits with 98 and 99% blastp identity; over 20,000 blastp results with 100% blastp identity.
* `conserved_vs_total_genes.png` graph showing number of conserved genes to number of total genes in pan-genome (no surprise, more total genes than conserved genes)
* `unique_vs_new_genes.png` graph showing number of new genes vs unique genes (more unique than new genes)

#### Files generated:
* accessory.header.embl			
* core_alignment_header.embl
* accessory.tab				
* core_gene_alignment.aln
* accessory_binary_genes.fa		
* gene_presence_absence.Rtab
* accessory_binary_genes.fa.newick
* gene_presence_absence.csv
* accessory_graph.dot			
* number_of_conserved_genes.Rtab
* blast_identity_frequency.Rtab		
* number_of_genes_in_pan_genome.Rtab
* clustered_proteins			
* number_of_new_genes.Rtab
* core_accessory.header.embl		
* number_of_unique_genes.Rtab
* core_accessory.tab			
* pan_genome_reference.fa
* core_accessory_graph.dot		
* summary_statistics.txt
* pan_genome_results_accessory
* pan_genome_results_core
* pan_genome_results_union


## Pangenome analysis with PPanGGOLiN
* Summary: ran PPanGGOLiN on FS19C gff data by running in prokka_env conda environment on Ceres to generate pangenome analysis of *E. coli* isolates.
* PPanGGOLiN publication DOI: https://doi.org/10.1371/journal.pcbi.1007732
* Github: https://github.com/labgem/PPanGGOLiN/wiki/Basic-usage-and-practical-information, https://github.com/labgem/PPanGGOLiN
* Began on: 3Feb2021
* Completed on:
* Platform: Ceres, prokka_env conda

1. (3Feb2021) Copy gbk files from 95 isolates on Ceres and place in `project/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gbk` directory
```
find . -name *.gbk -exec cp '{}' "./prokka_gbk/" ";"
```

2. (3Feb2021) Download gbk files from GenBank of the 6 E. coli reference strains via RefSeq (source database) and Genomic GenBank format .gbff (file type). PPanGGOLiN says can use .gbff/.gbk files. Uploaded 6 E. coli reference genomes' gbff files to `project/fsepru/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gbk/Ecolireferencegenomes` directory on Ceres. Unzipped files with `gzip -d *.gz` and copied to `prokka_gbk` directory.
* E. coli MG1655 https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2
  * Saved as Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gbff.gz
* E. coli HS https://www.ncbi.nlm.nih.gov/assembly/GCF_000017765.1
  * Saved as Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gbff.gz
* E. coli Nissle 1917 https://www.ncbi.nlm.nih.gov/assembly/GCF_000714595.1
  * Saved as Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gbff.gz
* E. coli O157:H7 str. NADC 6564 https://www.ncbi.nlm.nih.gov/assembly/GCF_001806285.1
  * Saved as Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gbff.gz
* E. coli O157:H7 EDL933 https://www.ncbi.nlm.nih.gov/assembly/GCF_000732965.1
  * Saved as EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gbff.gz
* TW14588 https://www.ncbi.nlm.nih.gov/assembly/GCF_000155125.1/
  * Saved as Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gbff.gz

3. (3Feb2021) Created text file listing all filenames in `prokka_gbk` directory and downloaded text file `Ecoligbk.txt`. Removed unnecessary filenames from list and opened in R. In R, made organism list with 1st column containing unique organism name (use the gbk file name) and second column as path to location of gbk file. Save as txt file. See this [page](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) for details. R code used:
```
library(tidyverse)
tsv <- read_tsv('Ecoligbk.txt', col_names = FALSE)
tsv$X3 <- paste(tsv$X2, tsv$X1, sep= "/")
head(tsv$X3)
tsv2 <- tsv[c("X1", "X3")]
head(tsv$X3)
write_tsv(tsv2, "Ecoligbkpath.txt")
```

4. (3Feb2021) Upload `Ecoligbkpath.txt` to Ceres.

5. (4Feb2021) Install PPanGGOLiN in prokka_env environment on Ceres
```
salloc
module load miniconda
source activate prokka_env
conda install -c bioconda ppanggolin
```

6. (4Feb2021) Tested ppanggolin.slurm script on Ceres in the prokka_gbk directory using test.slurm script via allocate a debug node to see if ppanggolin.slurm script will run on salloc.
```
ppanggolin workflow --anno ORGANISMS_ANNOTATION_LIST
```

7. (4Feb2021) Ran ppanggolin.slurm on Ceres as a slurm job.
```
#!/bin/bash
#SBATCH --job-name=ppanggolin                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
module load miniconda
source activate prokka_env
ppanggolin workflow --anno Ecoligbkpath.txt
```
```
Submitted batch job 5524245
```

8. (5Feb2021) Will need to re-run ppanggolin by re-running prokka on 95 isolates with E. coli reference genomes to get same annotation. Then can take the gbk files and run through ppanggolin.

9. (9Feb2021) Copy gbk files from `polishedgenomesprokka_95isolates5refgenomes/` and place in `project/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gbk` directory
```
find ../polishedgenomesprokka_95isolates5refgenomes/ -name *.gbk -exec cp '{}' "./" ";"
```

10. (9Feb2021) Created text file listing all filenames in `prokka_gbk` directory and downloaded text file `Ecoligbk.txt`.

11. (10Feb2021) Opened `Ecoligbk.txt` in R, duplicated 1st column to make second column. Goal is to have 1st column containing unique organism name (use the gbk file name) and second column as path to location of gbk file. Save as txt file. See this [page](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) for details. R code used:
```
library(tidyverse)
tsv <- read_tsv('Ecoligbk.txt', col_names = FALSE)
tsv$X2 <- tsv$X1 #duplicated X1 column and name as X2
write_tsv(tsv, "Ecoligbkpath.txt")
```
In Excel, deleted the path name `/lustre/project/fsepru/kmou/FS19C/polished_genomes_100X/prokka_gbk/` in first column, and deleted `/lustre/` in second column. Saved tsv. Uploaded to Ceres. Ran `ppanggolin.slurm` and found out it was looking for `prokka_env` in home, which didn't exist. I accidentally deleted `.conda` from project directory so I had to reinstall prokka, PPanGGOLiN in `prokka_env` in project directory following the correct guidelines on SciNet: https://scinet.usda.gov/guide/conda/#user-installed-software-on-ceres-with-conda

12. (11Feb2021) Ran ppanggolin.slurm on Ceres as a slurm job.
```
#!/bin/bash
#SBATCH --job-name=ppanggolin                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#SBATCH --account fsepru
#Enter commands here:
set -e
set -u
module load miniconda
source activate /project/fsepru/kmou/prokka_env
ppanggolin workflow --anno Ecoligbkpath.txt
```
```
Submitted batch job 5545691
```
Job cancelled. Looked at stderr file and saw the only line:
`/project/fsepru/kmou/prokka_env/etc/conda/activate.d/java_home.sh: line 1: JAVA_HOME: unbound variable`
How to address this??

#### Files generated:



## Create phylogenetic tree with RAxML-NG
* Summary: ran ... on FS19C ... data by running in ... to generate phylogenetic tree of the core genome of isolates. See how related isolates are.
* RAxML-NG publication DOI: 10.1093/bioinformatics/btz305
* Github: https://github.com/amkozlov/raxml-ng
* Began on: 3Feb2021
* Completed on:
* Platform: Ceres ...

1. It is on bioconda: https://bioconda.github.io/recipes/raxml-ng/README.html

#### Files generated:



## Metabolic pathways in pangenome with gapseq

## Screen for AMR or virulence genes with abricate


## WGS submission to SRA
* Must complete Biosample entry (which will generate biosample entry in tandem)
* [SRA Quick Start Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submit/) and a more detailed [SRA Submission Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submitbio/)
