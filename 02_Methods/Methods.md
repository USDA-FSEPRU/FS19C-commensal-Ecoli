# Methods
Details of sequence analyses methods performed on FS19C samples 1-96 in this repo. Includes lab notes of how methods performed

## Notes
GFF file description: http://gmod.org/wiki/GFF3

## Conda environment (updated on 11Feb2021 - condensed to one conda environment)
* **Make sure when calling environments, to use this path: /project/fsepru/kmou/conda_envs/**
See: https://scinet.usda.gov/guide/conda/#user-installed-software-on-ceres-with-conda
* Also can use conda environments loaded on `/project/fsepru/conda_envs/`

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

### Ceres
* Learned that the more resource parameters you specify on Ceres, it will try to make those work rather than doing what works best. So limit your parameters to let Ceres decide how to allocate resources to run job.
* Logical core (https://scinet.usda.gov/guide/ceres/#partitions-or-queues) includes hyperthreading
* Make sure when you run a job from a project directory, the owner of the directory is `proj-fsepru` and not your user name (will exceed home directory disk quota)!

## shortcuts in .bashrc
alias debug='salloc -N 1 -p debug -t 01:00:00'
alias myjobs='squeue | grep kathy.mo'

## (1) Sample Collection
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

## (2) Sequence assembly
* Summary: QC and assemble sequences using BBMap and spades
* Began on: 16Sept2020
* Completed on: 18Dec2020
* Platform: Ceres
  * /project/fsepru/kmou/FS19C/**

1. Make text file replace.sh

<details><summary>`replace.sh` script</summary>
```
#usr/bin/bash
while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
```
</details>

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

4. Run SRAassemblyPipeline.FS19C.SLURM_TEMPLATE script.

<details><summary>SRAassemblyPipeline.FS19C.SLURM_TEMPLATE script</summary>

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
</details>

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

## (3) QC with FastQC
* Summary: Ran fastqc on FS19C sequence data to assess sequence quality of individual reads for each sample, and to use output for multiqc (forgot to do this prior to sequence assembly)
* Began on: 6Jan2021
* Completed on: 7Jan2021
* Platform: Ceres, slurm
  * /project/fsepru/kmou/FS19C/**

1. Run `fastqc.slurm` script

<details><summary>fastqc.slurm script</summary>

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
</details>

##### Files generated (for each isolate):
  * *fastqc.zip
  * *fastqc.html

## (4) QC with MultiQC
* Tutorial: https://www.youtube.com/watch?v=qPbIlO_KWN0
* Summary: Ran multiqc with fastqc output of FS19C sequence data to assess quality of all sequences for samples 1-94 (forgot to do this prior to sequence assembly). I did not include samples 95 and 96 in the multiqc run as I realized these samples were ran on a second sequencing run, so their coverage is different than the first sequencing run that had all samples (but I would only consider  samples 1-94 from that first run). In addition, examining fastqc output for samples 95 and 96 was fine enough.
* Began and Completed on: 7Jan2021
* Platform: Ceres, fastanienv conda environment
  * /project/fsepru/kmou/FS19C/**

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

## (5) QC with FastANI
* Summary: ran fastANI on FS19C sequence data by running in conda environment to estimate Average Nucleotide Identity (ANI) using alignment-free approximate sequence mapping. It calculates distance between 2 sequences. Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* FastANI publication: DOI: 10.1038/s41467-018-07641-9
* Began on: 28Dec2020
* Completed on: 11Jan2021
* Platform: Ceres, fastanienv conda environment
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/**

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

## (6) QC with Mash
* Summary: ran Mash on FS19C sequence data by running in conda environment to compare results with fastANI. The *sketch* function converts a sequence or collection of sequences into a MinHash sketch. The *dist* function compares two sketches and returns an estimate of the Jaccard index, *P* value, and Mash distance (estimates rate of sequence mutation under a simple evolutionary model). Also need to include reference genomes to see how all sequences cluster relative to one another and if there are any outliers. Jules had mentioned fastANI is more accurate than Mash, but Mash is faster.
* Mash publication: DOI: 10.1186/s13059-016-0997-x
* Source: http://mash.readthedocs.org
* Began and Completed on: 11Jan2021
* Platform: Ceres, mashenv conda environment
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/**

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
* Output generated: .msh

3. Run distance estimates comparing sketch genomes with itself. This runs really fast (less than a second)
```
mash dist -p 1 .msh .msh > distances.tab
```

4. distances.tab columns: reference_id, query_id, mash_distance, pvalue, matching_hashes

5. (10Feb2021) Run Mash with all 95 isolates, and 5 reference genomes including TW14588 on slurm with `mash.slurm` script.

<details><summary>mash.slurm script</summary>

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
</details>

6. Downloaded `distances_secondrun.tab` to local computer, moved file to `Files/` and import to `mash_mds.R` script and created MDS. The plot is pretty much the same as previous mash results with TW14588 buried within the large cluster. Not sure why roary is giving so many errors.

7. (17Feb2021) Re-run mash to also include additional reference genomes to contrast related vs unrelated strains with the 95 E. coli isolates. Downloaded RefSeq fna files of the following 3 organisms and uploaded to `/project/fsepru/kmou/FS19C/polished_genomes_100X/referencegenomes/`
* (Enterobacteriaceae) Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000006945.2
  * Saved as Styphimurium_LT2.fna.gz
* (Non-Enterobacteriaceae, same phylum) Campylobacter jejuni subsp. jejuni NCTC 11168 aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000009085.1
  * Saved as Cjejuni_11168.fna.gz
* (Non-Proteobacteria, a Firmicutes) Clostridium saccharoperbutylacetonicum N1-4(HMT) aka https://www.ncbi.nlm.nih.gov/assembly/GCF_000340885.1
  * Saved as Clostridium_N1-4.fna.gz

8. (17Feb2021) Ran `rename .fna .fasta *.fna` to change *.fna extension to *.fasta. Created softlink of the three genomes from `/project/fsepru/kmou/FS19C/polished_genomes_100X/referencegenomes/` to `/project/fsepru/kmou/FS19C/polished_genomes_100X/mash_all`. Run `mash.slurm`

<details><summary>mash.slurm script</summary>

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
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
module load miniconda
source activate /project/fsepru/kmou/prokka_env
mash sketch -p 1 -o /project/fsepru/kmou/FS19C/polished_genomes_100X/mash_all/ *.fasta
mash dist -p 1 .msh .msh > distances_thirdrun.tab
```
```
Submitted batch job 5572662
```
</details>

##### Files generated:
  * distances.tab
  * FS19Cmashdistances.xlsx

## (7) QC: Visualize ANI pairwise genome-genome similarity calculations with MDS, heatmap
* Summary: Made distance matrix from mash and fastani output to create heatmap and MDS to visualize clustering and identify any outliers. The MDS was a bit hard to decipher what was an outlier, so I ran a heatmap to see how fastANI and mash compared and whether the pairwise comparisons were similar between the two, including heatmap of pearson correlation coefficients.
* Began on: 14Jan2021
* Completed on: 15Jan2021
* Platform: R Studio on local computer

1. See scripts/qc_mds.R for details

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

## (8) Genome annotation with prokka
* Summary: Identify annotation reference from UniProt (E. coli pangenome) and use prokka to annotate all 95 samples.
* Prokka publication: DOI: 10.1093/bioinformatics/btu153
* Began on: 20Jan2021
* Completed on: 28Jan2021
* Platform: Ceres, prokka_env conda environment
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/**

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

<details><summary>prokka.slurm script</summary>

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
</details>

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

## (9) Pangenome analysis with roary
* Summary: ran Roary on FS19C gff data by running in Ceres to generate pangenome analysis of *E. coli* isolates.
* Roary publication DOI: 10.1093/bioinformatics/btv421
* Began on: 29Jan2021
* Completed on: 24Feb2021 (completed with gifrop)
* Platform: Ceres
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/**

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

<details><summary>roary.slurm script</summary>

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
</details>

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

8. (1Feb2021) Modified roary.slurm script.

<details><summary>roary.slurm script</summary>

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
</details>

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

<details><summary>roary.slurm script</summary>

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
</details>

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

<details><summary>roary.slurm script</summary>

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
</details>

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
* **I will ask E. coli group if they want me to annotate with a specific reference genome. For now, I will keep going and analyze results as is.** Found out it was K-12, then EDL933.

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

24. (9Feb2021) Ran roary.slurm script

<details><summary>roary.slurm script</summary>

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
</details>

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

<details><summary>roary.slurm script</summary>

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
</details>

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

<details><summary>roary.slurm script</summary>

```
module load roary
query_pan_genome  -o pan_genome_results_union -v -a union *.gff
query_pan_genome  -o pan_genome_results_core -v -a intersection *.gff
query_pan_genome  -o pan_genome_results_accessory -v -a complement *.gff
```
```
Submitted batch job 5550887
```
</details>

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

32. (3Feb2021) Run `create_pan_genome_plots.R` from [Roary github](https://github.com/sanger-pathogens/Roary/blob/master/bin/). Saved as `roary.R`. Generated many plots, below are comments of each plot made.
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


## (10) Pangenome analysis with PPanGGOLiN
* Summary: ran PPanGGOLiN on FS19C gff data by running in prokka_env conda environment on Ceres to generate pangenome analysis of *E. coli* isolates.
* PPanGGOLiN publication DOI: https://doi.org/10.1371/journal.pcbi.1007732
* Github: https://github.com/labgem/PPanGGOLiN/wiki/Basic-usage-and-practical-information, https://github.com/labgem/PPanGGOLiN
* Began on: 3Feb2021
* Completed on:
* Platform: Ceres, prokka_env conda
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka/renamed_contigs/prokka_gbk/**

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

7. (4Feb2021) Ran `ppanggolin.slurm` on Ceres as a slurm job.

<details><summary>ppanggolin.slurm script</summary>

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
</details>

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

12. (11Feb2021) Ran `ppanggolin.slurm` on Ceres as a slurm job.

<details><summary>ppanggolin.slurm script</summary>

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
</details>

* Job cancelled. Looked at stderr file and saw the only line:
`/project/fsepru/kmou/prokka_env/etc/conda/activate.d/java_home.sh: line 1: JAVA_HOME: unbound variable`
How to address this??

13. (12Feb2021) Modified `ppanngolin.slurm` slurm script. I've encountered the unbound variable issue before and it was fixed when I put in `set +eu` in slurm script. Also modified conda package path. Jules made one for ppanggolin in `/project/fsepru/conda_envs` directory so I will use that.
```
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/conda_envs/ppanggolin
ppanggolin workflow --anno Ecoligbkpath.txt
```
```
Submitted batch job 5566707
```
Job completed successfully.

14. (24Feb2021) Re-run PPanGGOLiN with new gbk files generated by prokka. Moved old copy of `Ecoli_K12_MG1655.gbk` to `/project/kmou/FS19C/polished_genomes_100X/referencegenomes/`. Copied gbk files from 95 isolates + 6 reference genomes *_pol/ to `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/prokka_gbk` directory.
```
find . -name *.gbk -exec cp '{}' "./prokka_gbk/" ";"
```
  * (13Apr2021) This `Ecoli_K12_MG1655.gbk` gbk file is the actual gbk file from NCBI (not the one generated from prokka). Try running gifrop prokka without calling a specific priority strain??

15. (24Feb2021) Created text files listing all filenames and their paths in `prokka_gbk` directory and downloaded text files `Ecoligbklist.txt` and `Ecoligbkpath.txt`.
```
ls -d -1 "$PWD/"*.gbk > Ecoligbkpath.txt
ls *.gbk > Ecoligbklist.txt
```

16. (24Feb2021) Opened `Ecoligbklist.txt` and `Ecoligbkpath.txt` and names of files in 1st column (containing unique organism name) and full directory name of files in second column (as path to location of gbk file). Save as `Ecoligbk.txt`. See this [page](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) for details. Uploaded to ceres.

17. (24Feb2021) Ran `ppanggolin.slurm` with this slurm script:
```
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/conda_envs/ppanggolin
ppanggolin workflow --anno Ecoligbk.txt
```
```
Submitted batch job 5589844
```

17. (24Feb2021) Job ran successfully. Downloaded files and examined `tile_plot.html`. It looked slightly different from previously ran ppangolin, but not too different.

#### Files generated:
* gene_presence_absence.Rtab       
* organisms_statistics.tsv  
* pangenomeGraph_light.gexf  
* projection/
* matrix.csv                       
* pangenomeGraph.gexf       
* pangenome.h5               
* tile_plot.html
* mean_persistent_duplication.tsv  
* pangenomeGraph.json       
* partitions/
* Ushaped_plot.html


## (11) Create phylogenetic tree with RAxML
* Summary: ran core_gene_alignment.aln (generated from roary) in raxml to generate phylogenetic tree of the core genome of isolates. Goal is to identify units of horizontal gene transfer within very closely related strains using a pangenome framework.
* RAxML-NG publication DOI: 10.1093/bioinformatics/btu033
* Github: https://github.com/stamatak/standard-RAxML/blob/master/README
* Began on: 12Feb2021
* Completed on: 24Feb2021
* Platform: Ceres
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/**

1. (12Feb2021) Ran raxml on Ceres with `raxml.slurm` script adapted from Jules, available here: `/project/fsepru/shared_resources/SLURMS`. Job completed successfully.

<details><summary>raxml.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=RAXML                              # name of the job submitted
#SBATCH -p short                                       # name of the queue you are submitting to
#SBATCH -N 1                                         # number of nodes in this job
#SBATCH -n 40                                       # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=120G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
# ENTER COMMANDS HERE:
module load raxml
# -T is threads
# -x takes a random seed and turns on rapid-bootstrapping
# -N is how to do bootstrapping 'autoMRE' is an algorithm that will automatically determine when enough bootstraps have been done
# -m is the model to use
# -f select algorithm: "-f a": rapid Bootstrap analysis and search for best-scoring ML tree in one program run
# -n is the name to use for the outputs
# -p is a random seed to use for parsimony searches
# -o is the name of the genome you want to use as an outgroup, must match exactly, not required
raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f a -n core_genome_tree_1 -s core_gene_alignment.aln -T 40 -x 7 -N autoMRE -p 7
#End of file
```
```
Submitted batch job 5566413
```
</details>

2. (18Feb2021) Re-run RAxML on Ceres with slurm script `raxml.slurm` in `project/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates5refgenomes/pan` where the latest `core_gene_alignment.aln` file is located
```
Submitted batch job 5573447
```

3. (19Feb2021) Visualized 3 trees in FigTree with `RAxML_bestTree.core_genome_tree_1`, `RAxML_bipartitions.core_genome_tree_1`, `RAxML_bootstrap.core_genome_tree_1`.
What's the difference between ML and bootstrap-created trees? https://www.biostars.org/p/226407/
```
ML search and bootstrap search are different things in RAxML. The ML search is performed to find the best-scoring ML tree among 20 ML trees calculated from different stepwise-addition parsimony starting trees.This is based on the original alignment. On the other hand, bootstrapping resamples the alignment positions randomly with replacement to the length of the original alignment (100 times in this case) and generates a single tree on each round to infer statistical support of the branches. Afterwards, the bootstrap values are placed on the corresponding branches of the best-scoring ML tree.
```
Defines output files: https://cme.h-its.org/exelixis/resource/download/NewManual.pdf

4. (24Feb2021) Re-run RAxML with new `core_gene_alignment.aln` that was generated from 19Feb2021 (latest `gifrop.slurm` run):
```
Submitted batch job 5589843
```

5. (24Feb2021) Downloaded trees and open with FigTree. Options selected on left panel of FigTree:
  * Trees > Order nodes
  * Select Tip Labels, Scale Bar
I notice the `RAxML_bestTree.core_genome_tree_1` and `RAxML_bipartitions.core_genome_tree_1` produce very similar trees (order nodes) compared to `RAxML_bootstrap.core_genome_tree_1`.

6. Next: Highlight isolates with no stx genes, virulence genes. Save png for tree from `RAxML_bestTree.core_genome_tree_1` as `bestTree_Xisolates.png`.

#### Files generated:
* RAxML_bestTree.core_genome_tree_1
* RAxML_bipartitionsBranchLabels.core_genome_tree_1
* RAxML_bipartitions.core_genome_tree_1
* RAxML_bootstrap.core_genome_tree_1
* RAxML_info.core_genome_tree_1

## (12) Extract genomic islands with gifrop
* Summary: ran gifrop2 (developed by Julian Trachsel. Gifrop2 = version 2 of gifrop) via slurm on Ceres to identify ‘genomic islands’ from roary pangenomes. See how related isolates are.
* Github: https://github.com/Jtrachsel/gifrop
* Began on: 12Feb2021
* Completed on: 24Feb2021
* Platform: Ceres
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/**

1. (12Feb2021) Made `gifrop2.slurm` script and ran in `/project/fsepru/kmou/FS19C/polished_genomes_100X/prokka_gff` where `*.gff` files and soft link for `gene_presence_absence.csv` is.

<details><summary>gifrop2.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=gifrop2                             # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/conda_envs/gifrop2
gifrop --get_islands
```
```
Submitted batch job 5566653
```
</details>

Job ran successfully. Downloaded `gifrop_out/` to local computer. I looked at `clustered_island_info.csv` and noticed all entries under columns `island_type,	RESISTANCE,	res_type,	vir_type,	plasmid_type,	viro_type,	megares_type` were `NA`. Jules said he'll check it out.

2. (17Feb2021) Jules found my error, woo hoo!
```
It looks like you made a little mistake on your bash substitutions when you called prokka, that’s why all your files have the weird .fasta%.fasta.gff
Suffix.  I’m not certain yet but I think that might be throwing a wrench into things.
Your error is in the prokka slurm script when you are assigning the bash variable ‘tag’
you need to use quotes and curly brackets for bash substitutions
you had:
tag=$file%.fasta
you need:
tag=“${file%.fasta}”
Once I renamed the fastas I was able to run gifrop just fine and got lots of output.
If you want to re-annotate these genomes once you have renamed them, Maybe you could try out the ‘pan_pipe’ script I include in gifrop.  
I have made a SLURM script for you with an example of how I would do things and I placed it at:
/project/fsepru/kmou/FS19C/polished_genomes_100X/GIFROP.slurm
```

#### Files generated:
* clustered_island_info.csv  
* gifrop.log                 
* my_islands/
  * abricate/
    * All_islands.megares2  
    * All_islands.ncbi  
    * All_islands.plasmidfinder  
    * All_islands.vfdb
    * All_islands.viroseqs
  * All_islands.fasta  
  * island_info.csv
* pan_with_island_info.csv
* figures/
  * island_length_histogram.png          
  * islands_per_isolate.png   
  * Number_of_occurances_secondary.png
  * islands_per_isolate_no_unknowns.png  
  * Number_of_occurances.png
* islands_pangenome_gff.csv  
* pan_only_islands.csv  
* sequence_data/
  * *_pol.fasta%.fasta_short.gff
  * _pol.fasta%.fasta.fna


## (13a) Run pan_pipe from gifrop to run prokka, roary, and gifrop altogether with Ecoli_K12_MG1655 for prokka annotation argument.
  * Summary: Re-run prokka because I forgot to set E. coli MG1655 genbank file as priority annotation when I first ran prokka. This slurm script was provided by Jules, which runs prokka, roary, and gifrop (developed by Julian Trachsel. Gifrop2 = gifrop version 2) via slurm on Ceres. It will annotate all with prokka in parallel (will do 24 genomes at a time, each with 1 thread), run roary and generate a core genome alignment, and with gifrop, it will extract, classify, and cluster genomic islands
  * Github: https://github.com/Jtrachsel/gifrop
  * Began on: 17Feb2021
  * Completed on: 24Feb2021
  * Platform: Ceres
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/**

1. (17Feb2021) Ran `gifrop.slurm`

<details><summary>gifrop.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=prokka                            # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 24                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#Enter commands here:
set -e
module load miniconda
source activate /project/fsepru/conda_envs/gifrop2
# you should run this from a folder containing all your assemblies in fasta format with the suffix ".fna"
# change "REFERENCE.GBK" to whatever the name of your desired reference genome is.
# will annotate all with prokka in parallel (will do 24 genomes at a time, each with 1 thread)
# run roary and generate a core genome alignment
# run gifrop and extract classify and cluster genomic islands
pan_pipe --prokka_args "--proteins Ecoli_K12_MG1655.gbk --cpus 1" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
```
```
Submitted batch job 5570779
```
</details>

Job failed because had message `Please rename your contigs OR try '--centre X --compliant' to generate clean contig names.` Need to add that argument to `pan_pipe.slurm`

2. (18Feb2021) Revised `GIFROP.slurm` with the following, removed all previous prokka output files including `panpipe_logs` and ran slurm job:
```
pan_pipe --prokka_args "--proteins Ecoli_K12_MG1655.gbk --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
```
```
Submitted batch job 5573240
```
Job completed successfully! Downloaded `gifrop_out/` and roary output.

3. (19Feb2021) Gifrop output still have NA for `clustered_island_info`. In lab meeting, Jules says he has a bash script for manually changing the contig names.

4. (19Feb2021) See `BashScriptLesson.md` for details about `rename_contigs` and for-loop bash script `~/scripts/renamecontigs.sh` to run `rename_contigs` in desired fasta file directory. I saved `rename_contigs` in `~/.bashrc`, copied *.fna files from `/project/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/` to a new directory `/project/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/` and ran `~/scripts/renamecontigs.sh` in this new directory.

5. (19Feb2021) Moved `GIFROP.slurm` to new directory and modified script via removing `--centre X --compliant`. Run job on slurm.

<details><summary>gifrop.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=panpipe                            # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 24                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=32G   # memory
  #SBATCH --mail-user=kathy.mou@usda.gov
  #Enter commands here:
  set -e
  module load miniconda
  source activate /project/fsepru/conda_envs/gifrop2

  pan_pipe --prokka_args "--proteins Ecoli_K12_MG1655.gbk --cpus 1" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
  ```
```
Submitted batch job 5576687
```
</details>

6. (24Feb2021) Downloaded `gifrop_out/` and roary output. Gifrop worked (yay!), `clustered_island_info.csv` showed virulence genes, etc. Will need to go through `clustered_island_info.csv` to search for specific gene groups (sugar utilization, virulence genes, etc)

<details><summary>gifrop summary</summary>

1. Need to rename suffix of `.fasta` files to `.fna`. Also need to modify contig IDs in fasta files for gifrop to work properly (able to call virulence genes, etc. from various databases). Jules showed me his `rename_contigs` script:
```
rename_contigs() {
        FILE=$1
        BASE=$2
        awk -v basev="$BASE" '/^>/{print ">"basev"_"++i;next}{print}' "$FILE" > "$BASE"_rename.fasta
        rm $FILE
        mv "$BASE"_rename.fasta $FILE
}
export -f rename_contigs
```

2. I added this script in my `.bashrc` profile in home directory on Ceres. Do `source ~/.bashrc` from wherever you are to update bashrc profile.

3. I created a for-loop to run `rename_contigs` on all fasta files of desired directory. See `BashScriptLesson.md` for-loop bash script `~/scripts/renamecontigs.sh` and other details to run `rename_contigs` in desired fasta file directory.

4. I copied *.fna files from `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/` to a new directory `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/` and ran `~/scripts/renamecontigs.sh` in this new directory.

5. (19Feb2021) Moved `GIFROP.slurm` to `/project/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/` and ran the following script on slurm. Used Ecoli str. K-12 substr. MG1655 as priority annotation for prokka.

  ```
  #!/bin/bash
  #SBATCH --job-name=panpipe                            # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 24                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=32G   # memory
  #SBATCH --mail-user=kathy.mou@usda.gov
  #Enter commands here:
  set -e
  module load miniconda
  source activate /project/fsepru/conda_envs/gifrop2

  pan_pipe --prokka_args "--proteins Ecoli_K12_MG1655.gbk --cpus 1" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
  ```


6. Download `gifrop_out/` and roary output. `clustered_island_info.csv` shows virulence genes, etc.
</details>

#### Files generated:
* **_pol/ or Ecoli_*/
  * *_pol.err
  * *_pol.faa
  * *_pol.ffn
  * *_pol.fna
  * *_pol.fsa
  * *_pol.gbk
  * *_pol.gff
  * *_pol.log
  * *_pol.sqn
  * *_pol.tbl
  * *_pol.tsv
  * *_pol.txt
  * proteins.faa
  * proteins.pdb
  * proteins.pot
  * proteins.ptf
  * proteins.pto
* pan/
  * *.gff
  * accessory_binary_genes.fa
  * accessory_binary_genes.fa.newick
  * _accessory_clusters
  * _accessory_clusters.clstr
  * accessory_graph.dot
  * accessory.header.embl
  * accessory.tab
  * blast_identity_frequency.Rtab
  * _blast_results
  * _clustered
  * _clustered.clstr
  * clustered_proteins
  * _combined_files
  * _combined_files.groups
  * core_accessory_graph.dot
  * core_accessory.header.embl
  * core_accessory.tab
  * core_alignment_header.embl
  * core_gene_alignment.aln
  * core_gene_alignment.aln.reduced
  * gene_presence_absence.csv
  * gene_presence_absence.Rtab
  * gifrop_out/
    * clustered_island_info.csv
    * figures/
      * island_length_histogram.png
      * islands_per_isolate_no_unknowns.png
      * islands_per_isolate.png
      * Number_of_occurances.png
      * Number_of_occurances_secondary.png
    * gifrop.log
    * islands_pangenome_gff.csv
    * my_islands/
      * abricate/
        * All_islands.megares2
        * All_islands.ncbi
        * All_islands.plasmidfinder
        * All_islands.vfdb
        * All_islands.viroseqs
      * island_info.csv
      * All_islands.fasta
    * pan_only_islands.csv
    * pan_with_island_info.csv
    * sequence_data/
      * *.fna
      * *_short.gff
  * _inflated_mcl_groups
  * _inflated_unsplit_mcl_groups
  * _labeled_mcl_groups
  * M7lUUryBzC/
    * *.gff.proteome.faa
  * number_of_conserved_genes.Rtab
  * number_of_genes_in_pan_genome.Rtab
  * number_of_new_genes.Rtab
  * number_of_unique_genes.Rtab
  * pan_genome_reference.fa
  * pan_genome_sequences/
  * summary_statistics.txt
  * _uninflated_mcl_groups

## (13b) Analyze gifrop output to narrow down list of commensal E. coli isolates that don't possess any virulence factors (LEE, stx, hemolysin).

<details><summary>Details of attempts to exclude commensal E. coli isolates with virulence genes</summary>

1. (19Mar2021) Went through gifrop csv files:
* Clustered_island_info.csv: virulence and ARG genes, clusters of genes wrapped together (no annotation)
* Pan_with_island_info.csv: annotation of genes identified in each isolate (wide-format), similar info as Islands_pangenome_gff.csv but in different orientation. Missing some columns like source, type, start, end, score, strand, phase, attributes, ID, product, num_locus_tag, loc_tag_order, only_island, seqid_len, flanking_genes
* Islands_pangenome_gff.csv: annotation of genes identified in each isolate (long-format), similar info as pan_with_island_info.csv but in different orientation. Missing some columns like Pcluster, all_Sclusters, all_Tclusters, all_Qclusters
* Other less important files:
  * Pan_only_islands.csv: same as pan_with_island_info.csv (16345 rows), but missing pan info (11781 rows)
  * Island_info.csv: only info about islands (no annotation, 12361 rows)
* Jules mentioned in lab meeting on 19Mar2021 that gifrop is not designed to find all the metabolic pathways. It looks for genomic islands (not including core genome that all strains share) using the following databases:
  * Megares (v2.0, 2021-Jan-20) (contains ~8000 AMR genes and annotation): [annotations with notes](https://megares.meglab.org/download/megares_v2.00/megares_full_annotations_with_notes_v2.00.csv)
  * PlasmidFinder (v2.1?? 2020-Apr-19): plasmid typing -- identifies plasmids in total or partial sequences isolates of bacteria
    ```
    PlasmidFinder and pMLST: in silico detection and typing of plasmids.
    Carattoli A, Zankari E, Garcia-Fernandez A, Voldby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
    Antimicrob. Agents Chemother. 2014. April 28th.
    ```
  * vfdb (2020-Apr-19): virulence factors of pathogenic bacteria
  * ProphET aka viroseqs (2021-Jan-20): prophage sequence prediction tool that Jules customized to convert nucleotide tmp files into abricate formatted database. He uses viroseqs to detect phage genes in bacterial assemblies using abricate
  * NCBI (2020-Apr-19)
* So I will need to use other tools to find the relevant metabolic genes.

2. (19Mar2021) Screen for stx genes
* Which isolates are stx-negative? Looked at `gene_presence_absence.csv`
  * When I presented to CRIS group, I mentioned I finding 14 stx- isolates. Now the list has grown bigger.
  * Found 36 isolates that were stx-negative: 1-15, 26, 27, 30, 31, 55-62, 69, 72, 73, 79, 84, 86, 87, 95, 96

3. (8Apr2021 and 12Apr2021) Screen for hemolysin and LEE operon genes. Looked at `gene_presence_absence.csv` and `clustered_island_info.csv`

| Operon | virulence genes in ORF | Present in which of the 95 isolates? |
| -- | -- | -- |
| LEE1 | ler (LEE-encoded regulator)| not detected in `gene_presence_absence.csv` or `clustered_island_info.csv` |
| LEE1 | escRSTU (T3SS) | not detected in `gene_presence_absence.csv` |
| LEE2 | sepZ (T3SS)| not detected in `gene_presence_absence.csv` |
| LEE3 | NA | NA |
| Tir | tir promoter controls intimin expression via polycistronic operon containing these genes: tir | yes (only the Ecoli reference strains NADC6564, O157:H7, TW14588) |
| Tir | eae (attaching and effacing aka intimin) | yes (only the Ecoli reference strains NADC6564, O157:H7, TW14588) |
| Tir | escD | not detected in `gene_presence_absence.csv` |
| Tir | cesT ((T3SS LEE chaperone)) | yes (only the Ecoli reference strains NADC6564, O157:H7, TW14588) |
| LEE4 | espADB | not detected in `gene_presence_absence.csv` |
| LEE4 | espF | yes (only the Ecoli reference strains NADC6564, O157:H7, TW14588) |
| EHEC and EPEC plasmids | hlyCABD (hemolysin) | yes, 23 have it including: 18, 23, 24, 28, 29, 34, 35, 44, 45, 47, 51, 52, 74-77, 79, 80-83, 91, 92, TW14588. None of these isolates are in the stx-negative list |
| EHEC and EPEC plasmids | tagA | not detected in `gene_presence_absence.csv` or `clustered_island_info.csv` |
| EHEC and EPEC plasmids | espC | not detected in `gene_presence_absence.csv` or `clustered_island_info.csv`|
| EHEC and EPEC plasmids | bfp (bundle forming pili) | not detected in `gene_presence_absence.csv` or `clustered_island_info.csv` |
| NA | hlyE and hlyE_2 (hemolysin E) | yes, only these don't have hlyE: 19, 21-24, 38, 39, 42, 43, 46, 49, 50, 70, 76, 80, 81, 91, 92, Ecoli_Nissle1917. Only 79 is stx-negative; not detected in `clustered_island_info.csv`|

* No isolates have none of these virulence genes.

4. (8Apr2021) Created script `Leehly_virulence_gene_search.R` to automate finding isolates that don't possess LEE and hemolysin genes. Create heatmap to quickly visualize results.

5. (8Apr2021) Why was ler not detected in roary? Not present in `gene_presence_absence.csv`. Vijay says I can use EDL933 `ler` gene for blast search. How to do a parallel blast search? Ask Jules.

6. (12Apr2021) Search through vfdb for complete list of virulence genes in E. coli. Downloaded `Escherichia_VFs_comparison.xls`.

7. (12Apr2021) Sent an email to Vijay with this list of genes, asking which of these to focus on to narrow search. The list he sent me are the following (see [01_Background](https://github.com/k39ajdM2/Notebook/blob/main/01_Background.md) for complete email correspondence):
  ```
  The modified list would be: ler, escRSTU, sepZ, escD, espADB, tir, eae, cesT, espF, stcE, espC, hlyCABD, hlyE, pchABC, espP, efa1/lifA, toxB, stx1, and stx2 for you to use in Blast search...

  Yes, EDL933 will be a good reference strain to use in your Blast search.
  ```
  * 29 genes total

8. (12Apr2021) Created a bash script `search.sh` to look up and see if the genes are found in Ecoli_O157H7_EDL933 using results from `clustered_island_info.csv`.

```
#!/bin/bash
echo "Beginning of file" > isolates.txt
for line in `cat Ecolivirulencegene.txt`;
do
	echo $line >> output.txt
	grep -n "$line" O157h7.txt >> isolates.txt
	echo "END_OF_GENE" >> isolates.txt
done
echo "End of file" >> isolates.txt
```

9. (12Apr2021) Opened `output.txt` in BBEdit and did the following:
  * Find `\n` and replace all with ` , `. This causes everything to be on one line separated by ` , `
  * Find `, END_OF_GENE ,` with `\n`. This removes all instances of `, END_OF_GENE ,` and makes a new line. This way it is a lot easier to read the output.
  * Examined output and saw the following genes were not identified in Ecoli_O157H7_EDL933 based on `clustered_island_info.csv`:
  ```
  stcE
 espC
 hlyC
 hlyA
 hlyB
 hlyD
 hlyE
 pchA
 pchB
 pchC
 espP
 efa1
 lifA
 toxB
  ```

10. (12Apr2021) Looked for these genes in Ecoli_O157H7_EDL93 `gene_presence_absence.csv`:

| Question | Genes |
| -- | -- |
| Did not find these listed in `Gene` column| ler, espC, pchA, pchB, pchC, efa1, lifA, toxB |
| Not present in Ecoli_O157H7_EDL93 | stcE, hlyC, hlyA, hlyB, hlyD,  espP |
| Present in Ecoli_O157H7_EDL93 | hlyE |

11. (12Apr2021) Looked through Escherichia_VFs_comparison.xls and saw that EDL933 in their database:
* Should not have: espC, pchA, pchB, pchC
  * **Need to ask Vijay about this**
  * Perhaps espC, pchABC found in other E. coli O157:H7 strains but not in the ones I looked at?
* Should have: ler, efa1, lifA, toxB (plasmid only)
  * Why does my Ecoli_O157H7_EDL933 not have these genes?

12. (12Apr2021) Changed approach to search all the virulence genes Vijay requested with `search.sh` in all isolates in  `clustered_island_info.csv` file. The genes I did not find in `clustered_island_info.csv` for all isolates were:
```
ler
espC
hlyE
pchA
pchB
pchC
efa1
lifA
```

The genes I did not find in `gene_presence_absence.csv` for all isolates were:
```
ler
escR
escS
escT
escU
sepZ
escD
espA
espD
espB
espC
pchA
pchB
pchC
efa1
lifA
toxB
```

13. (12Apr2021) Saw on uniprot that `ler` can also be `hns` (https://www.uniprot.org/uniprot/Q7DB51). Found hns gene in csv but Ecoli_Nissle1917 has it listed along with other genes under the "genes" column.
```
hns_2|group_2882|group_2879|group_3880|group_3881|group_3882|group_5366|group_3953|group_3954|group_3955|group_5559|group_5558|group_9873|group_9874|group_5725|group_13768|group_13769|group_13770|group_13771|group_13772|group_13773|group_13774|group_13775|group_13776|group_13777|group_2880
```

14. (13Apr2021) Spoke with Jules in the morning (9-9:30AM) and found out the following:
  * gifrop and roary may not call all genes. With gifrop, only genes that pass some threshold would show up. With roary, the genes may be close to being called a known gene name, but instead given an obsure name like `group_1805` because can't definitively call it the common gene name (like ler).
  * can run `roary -e n z` to get multialignment fasta of core genes
  * explained column names of `clustered_island_info.csv`: percent_island, only_island, genes, res_type, megares_type, vir_type, viro_type
  * why vfdb didn't call some of the expected virulence genes in EDL933: again, the sequences may not match close enough to call it that gene name. May call it group_xxxx.
  * Try running blast using custom blast database of virulence genes to your isolates to make sure they are or aren't there.
  * **File that has pan-genome fasta core and accessory genes of all isolates: pan_genome_reference.fa**

15. (13Apr2021) E. coli CRIS meeting had some suggestions:
  * Does prokka annotation file have gene names? Check out EDL933 for ler
  * Try running prokka with EDL933 as priority strain to see if that fixes the miscalled gene issue. Also try running prokka without setting priority strain (just specify Escherichia coli) and see what that does for results.
    * Do we need to run two different annotation pipelines - one with K-12 and the other with EDL933?
  * Try BioCyc to characterize E. coli metabolic pathways - subscription required… Iowa State has one, but not sure how we can access.
  * DRAM - email paper to everyone, DONE.
  * Ask Vijay or Indira if they can recommend any papers that did Nissle & K-12 comparative genomics/metabolic pathways to check out.

16. (13Apr2021) Looked for ler in prokka annotation file of EDL933 in `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/pan`. I noticed `Ecoli_O157H7_EDL933.gff` has `similar to AA sequence:Ecoli_K12_MG1655.gbk` for each gene annotation line. Will need to re-run gifrop with EDL933 as priority strain.

17. (13Apr2021) Found EDL933 gbff file (gbk) locally: `EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gbff.gz`
  * It has the genes `ler, toxB`
  * It does not have the genes: `pchABC, espC, efa1, lifA`

18. (13Apr2021) Checked a new gbff copy from here and did not find the `pchABC, espC, efa1, lifA` genes. ASM73296v1 (UC San Diego): https://www.ncbi.nlm.nih.gov/assembly/GCF_000732965.1
  * It does have:
    ```
    Ler
    escD
    escRSTU
    espF
    espADB
    Tir
    Eae
    cesT
    hlyD
    hlyE
    toxB
    stxA1
    stcE
    espP
    ```
  * It doesn't have:
    ```
    hlyA
    hlyB
    pchABC
    Efa1/lifA
    stxA2, stxB1, stxB2
    sepZ
    espC
    ```

19. (14Apr2021) Look for other versions of EDL933 genome assembly:
  * ASM666v1 (Wisconsin): https://www.ncbi.nlm.nih.gov/assembly/GCF_000006665.1
    * no `pchABC, espC, efa1, lifA`
  * ASM302871v1 (Rehman Medical Institute): https://www.ncbi.nlm.nih.gov/assembly/GCF_003028715.1
    * no `ler, pchABC, espC, efa1, lifA`
  * ASM976618v1 (Ohio State): https://www.ncbi.nlm.nih.gov/assembly/GCF_009766185.1
    * no `pchABC, espC, efa1, lifA`
  * ASM94844v1 (Canadian Food Inspection Agency): https://www.ncbi.nlm.nih.gov/assembly/GCF_000948445.1
    * no `pchABC, espC, efa1, lifA`
  * If none of these have the virulence genes, maybe we can create a custom blast database with genes from EDL933 that it does have, and the other genes from other E. coli strains. Or look through vfdb `Escherichia_VFs_comparison.xls` for potential strains that could have all the genes.

20. (13Apr2021) Looked up E. coli O157:H7 str. Sakai https://www.ncbi.nlm.nih.gov/assembly/GCF_000008865.2#/def:
  * It does have:
  ```
  Ler
  escRSTU
  escD
  espADB
  Tir
  Eae
  cesT
  espF
  hlyCABD
  hlyE
  pchABC
  Efa1/lifA
  toxB
  stx1
  stx2
  ```
  * It doesn't have:
  ```
  sepZ
  stcE
  espC
  espP
  ```

21. (13Apr2021) Looked up espC in https://www.ncbi.nlm.nih.gov/gene/?term=espC+AND+escherichia+coli and search results are all discontinued except for this: https://www.ncbi.nlm.nih.gov/gene/12682595
  * **Ask Vijay for help**

22. Try running prokka without setting priority strain (just specify Escherichia coli) and see what that does for results? Caveat naming scheme may not be the same (https://github.com/tseemann/prokka)
  ```
  Option: --proteins

  The --proteins option is recommended when you have good quality reference genomes and want to ensure gene naming is consistent. Some species use specific terminology which will be often lost if you rely on the default Swiss-Prot database included with Prokka.

  If you have Genbank or Protein FASTA file(s) that you want to annotate genes from as the first priority, use the --proteins myfile.gbk. Please make sure it has a recognisable file extension like .gb or .gbk or auto-detect will fail. The use of Genbank is recommended over FASTA, because it will provide /gene and /EC_number annotations that a typical .faa file will not provide, unless you have specially formatted it for Prokka.
  ```
Also run search.sh on all EDL933 to see which ones have more than the other.
23. (14Apr2021) Also check out `Escherichia_VFs_comparison.xls` to find strain
  * Found 25 genes except for `pchABC, stcE`

24. (14Apr2021) Met with Bradd Haley, Jo Ann Kessel, Chris, and Crystal. I discussed the issue of selectively choosing a strain to prioritize annotation (Ecoli_K12_MG1655), but that causes some genes to be missed (virulence genes) and not called by the common gene name. Another issue is even if I choose EDL933, I still have the same issue that some genes are missing or they're not called the same among other EDL933 strains. Do I make my own custom annotation database or ...? Bradd suggested running prokka against standard database, or I can specify Escherichia coli, and then compare the two pan-genomes (STEC vs non-STEC) and run fisher exact test to see which genes are enriched in one group or another (less biased approach).
  * new approach:
      * run roary of commensals and STEC separate but make sure they have the same annotation approach?
      * How to do Fisher exact test?
      * DRAM
  * current approach:
      * run prokka with Escherichia coli of commensals + STECs via gifrop
      * parallel blast to make sure virulence genes exist in the strains
        * custom database of virulence genes: track which source (EDL933, Sakai, etc.)
        * some genes...
</details>

## (13c) Run pan_pipe from gifrop to run prokka, roary, and gifrop altogether using Escherichia coli prokka argument.
  * Summary: Re-run prokka setting Escherichia coli as priority annotation. This slurm script was provided by Jules, which runs prokka, roary, and gifrop (developed by Julian Trachsel. Gifrop2 = gifrop version 2) via slurm on Ceres. It will annotate all with prokka in parallel (will do 24 genomes at a time, each with 1 thread), run roary and generate a core genome alignment, and with gifrop, it will extract, classify, and cluster genomic islands
  * Github: https://github.com/Jtrachsel/gifrop
  * Began on: 14Apr2021
  * Completed on:
  * Platform: Ceres
  * /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/gifropWithEcoliAnnotation/

1. (14Apr2021) Renamed `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/` to `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/` to make it clear that the `*.fna` files in this directory have not been annotated by prokka. They use to be `*.fasta` but were changed to `*.fna` to match with prokka commands from github page (see Section 8, #6 (21Jan2021) entry for details).

2. (14Apr2021) Make a new directory called `gifropWithEcoliAnnotation` and link `*.fna` files from `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs` to `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/gifropWithEcoliAnnotation/`
```
 ln -s /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/*.fna .
```

3. (14Apr2021) Run gifrop with Escherichia coli specification for prokka argument. See `gifrop.slurm` for script details. Job # 5751498.

<details><summary>gifrop.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=Ecolipanpipe                           # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 24                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=32G   # memory
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
module load miniconda

source activate /project/fsepru/conda_envs/gifrop

pan_pipe --prokka_args "--genus Escherichia --species coli --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
```
</details>

4. (14Apr2021) Make `blastparallel.slurm` and modify script. Need to determine how to make database.

<details><summary>blastparallel.slurm script</summary>

  ```
  #!/bin/bash

  #SBATCH --job-name=virgene_BLAST                          # name of the job submitted
  #SBATCH -p short                                        # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 72                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name, %x adds job name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=64G   #
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  # ENTER COMMANDS HERE:

  module load blast+
  module load parallel

  mkdir virulencegenes

  find /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/gifropWithEcoliAnnotation/ -name "*.fna" | parallel 'blastn -db ~/blastdb/METAL_ISLAND.fasta -query {} -out ./virulencegenes/{/.}.virulence -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore"'

  cd virulencegenes

  # collects and combines all blast results into one file
  find . -name "*virulence" | xargs -n 100 cat >> ALL_virulence.results

  # removes all files ending with 'metal'
  find . -name "*virulence" | xargs rm

  #End of file
  ```
</details>

5. (16Apr2021) Looked at `gene_presence_absence.csv` and `clustered_island_info.csv` for the respective missing virulence genes from the first run. Didn't find the virulence genes. Could make custom blast database of desired virulence genes, but will follow Bradd Haley's approach instead.

6. (19Apr2021) Downloaded `EDL933.fasta` from `/project/fsepru/data_transfer/O157_challenge_strains/` to `/project/fsepru/kmou/FS19C/EDL933prokkatest/` to see when I run prokka on this via `prokka.slurm`, will `ler` and other virulence genes show up. Submitted job 5758134.

<details><summary>prokka.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=prokka                            # name of the job submitted
#SBATCH -p short                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:

module load prokka
prokka --genus Escherichia --species coli --cpus 1 --centre X --compliant --outdir prokka_out EDL933.fasta
```
</details>


7. (20Apr2021) The job failed at
`Could not run command: tbl2asn -V b -a r10k -l paired-ends -M n -N 1 -y 'Annotated using prokka 1.14.5 from https://github.com/tseemann/prokka' -Z prokka_out\/PROKKA_04192021\.err -i prokka_out\/PROKKA_04192021\.fsa 2> /dev/null`. I looked through github issue page: https://github.com/tseemann/prokka/issues/139 and saw it was the version of tbl2asn used. In my stderr file, I saw a couple lines saying the version I had was outdated and needed current version. I emailed VSRC for help on how to proceed.

8. (21Apr2021) Yasasvy fixed the issue. I'm re-running `prokka.slurm` job 5759514. Looked through `PROKKA_04212021.gbk` and could not find ler, espC, pch, sepZ, etc. So this did not work... will have to result to blast.

9. (26Apr2021) Save Bradd Haley's R script as `BraddHaleyOriginalScript.R` and `FisherExactTest.R` in `scripts` folder.

10. (26Apr2021) Created `EcoliAnnotation_gene_presence_absence.Rtab.csv` from `gene_presence_absence.Rtab` generated in `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/gifropWithEcoliAnnotation/pan`. Added a row called `type` that designates each sample as `commensal` or `STEC`. Column name for genes is labeled as `ID`.

11. (26Apr2021) Run `EcoliAnnotation_gene_presence_absence.Rtab.csv` through `FisherExactTest.R`.

12. (27Apri2021) Emailed Bradd about what the `PanGenomeEcoli.csv` file is from the line `genomes <- read.csv("PanGenomeEcoli.csv",header=FALSE,skip=2)` in `FisherExactTest.R`.

13. (27Apr2021) Look for STEC genomes on NCBI (RefSeq)
* Read https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/README.txt and followed directions on how to grab assembled sequence data of STEC strains from NCBI.
* Looked for latest assembly versions: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/README.txt
```
all_assembly_versions & latest_assembly_versions are not listed for any species with more than 1,000 assemblies.
Use the assembly_summary.txt file in the species directory to find assemblies of interest
and then access the data using the ftp_path from column #20 of the file.
```
* Download `assembly_summary.txt` via `wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/assembly_summary.txt`
* Saved as `assembly_summary.xlsx`. Sheets include:
  * `assembly_summary`: all assembly data
  * `Complete Genome`: filtered from `assembly_summary` sheet, only includes complete genomes based on column L `assembly_level`.
  * `O157H7`: filtered from `Complete Genome` sheet, only includes  `Escherichia coli O157:H7` based on column H `organism_name`.
  * total of 133 STEC genomes

14. (27Apr2021) Fetch 133 STEC genomes from NCBI
* Tested wget with the first two and the files don't exist. Tried running slurm script to see what can be fetched - nothing was picked up. Asking Jules if he's encountered this or should I contact NCBI help desk.

<details><summary>wget and stec.ftp.slurm script details</summary>
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/307/215/GCF_001307215.1_ASM130721v1
--2021-04-27 15:48:41--  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/307/215/GCF_001307215.1_ASM130721v1
           => ‘GCF_001307215.1_ASM130721v1’
Resolving ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)... 165.112.9.229, 130.14.250.7, 2607:f220:41e:250::11, ...
Connecting to ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)|165.112.9.229|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /genomes/all/GCF/001/307/215 ... done.
==> SIZE GCF_001307215.1_ASM130721v1 ... done.
==> PASV ... done.    ==> RETR GCF_001307215.1_ASM130721v1 ...
No such file ‘GCF_001307215.1_ASM130721v1’.
```
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/995/GCF_001558995.2_ASM155899v2
--2021-04-27 15:49:03--  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/995/GCF_001558995.2_ASM155899v2
           => ‘GCF_001558995.2_ASM155899v2’
Resolving ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)... 165.112.9.229, 130.14.250.7, 2607:f220:41e:250::11, ...
Connecting to ftp.ncbi.nlm.nih.gov (ftp.ncbi.nlm.nih.gov)|165.112.9.229|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /genomes/all/GCF/001/558/995 ... done.
==> SIZE GCF_001558995.2_ASM155899v2 ... done.
==> PASV ... done.    ==> RETR GCF_001558995.2_ASM155899v2 ...
No such file ‘GCF_001558995.2_ASM155899v2’.
```
  ```
  #!/bin/bash
  #SBATCH --job-name=stec                             # name of the job submitted
  #SBATCH -p short                                 # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 2                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  #Enter commands here:
  wget -i stecftp.txt
  ```
</details>

15. (28Apr2021) Jules figured out I was fetching the directory and not any particular file. `wget` expects to fetch a file. He also pointed out that I can access the directory depending on what browser I used. He used chrome but couldn't access it. So he switched to Microsoft Edge and it worked. I used safari and got "can't open page". I assumed the link was broke, but when I tried chrome, I could access the directory. Found out I want the file that ends with `_genomic.fna.gz`. Edited `assembly_summary.xlsx` via creating columns `desired_file_from_directory` and `new_ftp_path` and combining columns to create new path. Copied new paths to `stecftp.txt` and uploaded to Ceres.

16. (28Apr2021) Ran `stec.slurm` to copy `*genomic.fna.gz` from RefSeq to `/project/fsepru/kmou/FS19C/STECgenomes`. Unzip files via `gzip -d *genomic.fna.gz`.

17. (28Apr2021) Accidentally deleted all contents of directory `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/`. I thought I wanted to sym link the `*pol.fna` files within `renamed_contigs` in `/project/fsepru/kmou/FS19C/STECgenomes` but then realized I wanted another copy of the `*pol.fna` in this directory. So I deleted the sym link with `rm -r 96commensalEcoli/`. Whoops. I should've done `rm 96commensalEcoli` or use `unlink 96commensalEcoli`. I did manage to save the `*_pol.fna` files from `rename_contigs/` (where the contig names were renamed to work with prokka in gifrop). I copied them over from `/project/fsepru/kmou/FS19C/STECgenomes` to `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/`. However, everything else is gone. The only things that I'm sad I lost were gifrop results (setting `--prokka_args "--genus Escherichia --species coli`) and the original `*_pol.fna` before I changed the contig names.

18. Run `gifrop.slurm` on `*pol.fna` files in `/project/fsepru/kmou/FS19C/STECgenomes`. Job 5801953.

<details><summary>gifrop.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=Ecolipanpipe                           # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 24                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=32G   # memory
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov
  #Enter commands here:
  set -e
  module load miniconda
  source activate /project/fsepru/conda_envs/gifrop
  pan_pipe --prokka_args "--genus Escherichia --species coli --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
  ```
</details>

19. (29Apr2021) Job completed. For the files I don't want transferred to local computer, I placed them in a directory `notransfer/`. For the rest of the files, I did `rsync -av --exclude notransfer -e ssh ceres:/project/fsepru/kmou/FS19C/STECgenomes/pan ./`

## 14. Screen for bacteriocins, microcins
* What are the genes for bacteriocins, microcins?
* BACTIBASE or Bagel4 for bacteriocin ID


## (15) Metabolic pathways in pangenome with gapseq
* Github: https://github.com/jotech/gapseq
* Began on: 15Mar2021
* Completed on:
* Platform: personal MacOS via HomeBrew installation, then on Ceres
* gapseq publication DOI: https://doi.org/10.1186/s13059-021-02295-1
* https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02295-1
* Initially didn't think to try gapseq on Ceres, so I tried it on my local computer. Ran into many issues with gapseq readLink and exonerate

<details><summary>gapseq on local computer</summary>

1. (15Mar2021) Followed installation instructions: https://gapseq.readthedocs.io/en/latest/install.html. Hit a couple snags and fixed with these commands:
```
brew update     # <= always update brew before installing anything new
brew upgrade
brew install coreutils binutils git glpk blast bedtools r brewsci/bio/barrnap grep bc
vim ~/.Rprofile
  options(repos=structure(c(CRAN="https://mirror.las.iastate.edu/CRAN/")))
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "CHNOSZ", "jsonlite"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

2. (15Mar2021) This might be the closest tutorial? https://gapseq.readthedocs.io/en/latest/tutorials/yogurt.html. Need to look up how to run on all 95 isolates at once. I found in issues page for gapseq you can run parallel: https://github.com/jotech/gapseq/issues/52. I installed parallel like so:
```
brew install parallel
# Next command taken from source
parallel ./gapseq find -e 5.4.99.17 {} ::: toy/*.fna.gz
```

3. (15Mar2021) Prepare fasta files from Ceres as a tar file. On Ceres, created new directory `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/polishedgenomes/` and do the following:
```
cp /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/*.fna polishedgenomes/
tar -czvf fs19cpolishedgenomes.tar.gz polishedgenomes/
du -sh polishedgenomes/   # <= see size of polishedgenomes/
```
Tar file size is 156M.
`polishedgenomes/` is 513M on Ceres.

4. (15Mar2021) Download tar file to local Mac, untar via `tar -xvf fs19cpolishedgenomes.tar.gz`.

5. (17Mar2021) Run gapseq with parallel, search for metabolic pathways:
```
parallel ./gapseq find -p all {} ::: ./gapseq/polishedgenomes/*.fna
```
Gapseq options:
* test = testing dependencies and basic functionality of gapseq
```
gapseq test
```
* find = pathway analysis, try to find enzymes based on homology
```
gapseq find (-p pathway | -e enzymes [-b bitscore] (genome)
gapseq find -p all toy/myb71.fna.gz
```

6. Encountered an error:
```
Invalid file: <path to file/*.fna>
readlink: illegal option --f
usage: readlink [-n] [file ...]
```
Looked up Issues page: https://github.com/jotech/gapseq/issues/28
Tried the following:
```
ln -s /usr/local/bin/greadlink /usr/local/bin/readlink
ln -s /usr/local/bin/gstat /usr/local/bin/stat
ln -s /usr/local/bin/ggrep /usr/local/bin/grep
```
7. Tried running gapseq on one file
```
gapseq find -p all <path to file/91-438FEC_pol.fna>
```
Encountered another error that gapseq could not find fastafetch or fastaindex. Found out I need to install exonerate: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate. I didn't think that I could do this on gapseq, but after speaking with some colleagues, I learned it is ok to try on Ceres as long as I don't need administrative permissions.
</details>

### On Ceres
1. (17Mar2021) Go to fsepru project directory, make installation directory (`programs/`), cd to `programs`, then `git clone https://github.com/jotech/gapseq` in `programs`.

2. (17Mar2021) Edit bashrc profile to include this line (install Exonerate and gapseq to specific path instead of to usr/local/bin):
```
export PATH="/project/fsepru/kmou/programs/bin/:/project/fsepru/kmou/programs/gapseq:$PATH"
source ~/.bashrc
```

3. (17Mar2021) Install Exonerate in `programs` and run the following commands to compile C code:
```
wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0.tar.gz
tar -xvf exonerate-2.2.0.tar.gz
./configure --prefix=<installationpath>
make
make check
make install
make clean    # <= removes temp files during install
```
You can also check out the `INSTALL` file for directions

4. (17Mar2021) Tested that gapseq and exonerate programs work (got an intro page).

5. (17Mar2021) Go to `project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes` and made `template.sh` template slurm script and `gapseq.sh` for loop with `gapseq` and `sbatch` commands to generate a slurm script with `gapseq` and `sbatch` commands for each of the 95 isolates + 6 reference genomes.

<details><summary>template.slurm and gapseq.sh scripts</summary>

**template.slurm**

  ```
  #!/bin/bash
  #SBATCH --job-name=gapseq                            # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 2                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-user=kathy.mou@usda.gov
  #Enter commands here:
  ```

**gapseq.sh**

  ```
  #!/bin/bash

  for file in ./*.fna
    do
      #cp template.sh $file.gapseq.sh
      #echo "gapseq find -p all $file" >> $file.gapseq.sh
      #echo "#End of file" >> $file.gapseq.sh
       sbatch "$file.gapseq.sh"
    done
  ```
</details>

6. (17Mar2021) Ran gapseq.sh, comment out `sbatch "$file.gapseq.sh"`. Then run gapseq.sh again, commenting out the following lines and only run `sbatch "$file.gapseq.sh"`
```
#cp template.sh $file.gapseq.sh
#echo "gapseq find -p all $file" >> $file.gapseq.sh
#echo "#End of file" >> $file.gapseq.sh
```

7. (17Mar2021) Submitted all 101 jobs on slurm at once (jobs 5616514-5616614).

8. (18Mar2021) Jobs completed, but browsed through Ecoli_TW14588, Ecoli_K12_MG1655 Pathways.tbl files and a few of the 95 isolates and did not detect any "true" predictions, there were no pathways detected... need to do more digging. See how gapseq examples were run. What were their fasta files like? Could it be because of how my fasta file is organized? Test again with Ecoli_TW14588, use different gapseq parameters? Do I need complete assembled genomes?
* `gapseq find-transport *.fna`

9. Try running gapseq on Streptococcus_thermophilus_strain_ATCC_19258.fna
`$gapseq find -p all -b 200 -m Bacteria Streptococcus_thermophilus_strain_ATCC_19258.fna`

## (16) DRAM
* Summary:
* Began on: 19Apr2021
* Completed on:
* DRAM doi: https://academic.oup.com/nar/article/48/16/8883/5884738
* https://github.com/shafferm/DRAM
* Platform: Ceres, miniconda

1. (19Apr2021) Before installing DRAM, I rearranged directories so that there's a `programs` directory for all non-conda things and a `conda_envs` for all things conda.

2. (19Apr2021) Install DRAM via miniconda: https://github.com/shafferm/DRAM/wiki/2.-How-to-Install-and-Set-Up-DRAM in .
```
cd /project/fsepru/kmou/conda_envs
wget https://raw.githubusercontent.com/shafferm/DRAM/master/environment.yaml
```

3. (19Apr2021) Edited `environment.yaml` to the following:
```
pkgs_dirs:
  - /project/fsepru/kmou/my_pkg_cache
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.*
  - pandas
  - pytest
  - scikit-bio
  - prodigal
  - mmseqs2!=10.6d92c
  - hmmer!=3.3.1
  - trnascan-se >=2
  - sqlalchemy
  - barrnap
  - altair >=4
  - openpyxl
  - networkx
  - ruby
  - parallel
  - dram
environment.yaml
```

4. (19Apr2021) Create DRAM conda environment and activate to test.
```
conda env create -f environment.yaml --prefix /project/fsepru/kmou/conda_envs/DRAM # <= -f creates an environment from environment.yaml
conda activate /project/fsepru/kmou/conda_envs/DRAM
```

5. (19Apr2021 and 20Apr2021) Made `dram.slurm` script to set up databases. Asked Chris for input and ran job (5758517).

<details><summary>dram.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH --mem=550gb
#SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs/DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data
```
</details>

<details><summary>Advice from Chris and Jules about DRAM and running slurm jobs</summary>

### Chris:
```
DRAM will make the output directory wherever you run the script from (don’t run it in your HOME directory though as there wont be enough storage). When you submit a job, you probably want to tell it how much total memory to allocate to the job with the –mem parameter. It says to setup DRAM you need at least ~512GB of RAM. The short partition doesn’t allow that much RAM it looks like, so you need to use the mem partition of the server (-p parameter). Usually when I request threads, I use --ntasks-per-core. I think you will typically want to request a number of threads instead of a number of cores (which is -n). Most bioinformatics software can only be run on a single core is my understanding but can use multiple threads on that core. If you want to change the default number of threads that DRAM uses, you need to add the --threads 16 to the command.
```
```
I would recommend setting up with the uniref database, as their paper shows that gave them the best results. It will take awhile to download though and needs a lot of RAM when you start the job – “Setting up DRAM can take a long time (up to 5 hours) and uses a large about of memory (512 gb) by default.”

To install uniref:
DRAM-setup.py prepare_databases --output_dir DRAM_data

When we ran DRAM, we just wanted something quicker, so we skipped uniref and setup the databases using this command:
DRAM-setup.py prepare_databases --output_dir DRAM_data --skip_uniref
```

### Jules:
```
If you are unsure you have specified your SLURM resource requests correctly, you can always check the number of processors available with `nproc` and the amount of memory available with `free -h`
I don’t often use the “--ntasks-per-core” flag, I don’t really understand what it does.  If I want 16 processors on a single node I will use `-N 1` and `-n 16`
If you think others in the unit would like to use this software (I know I would) you should consider installing the conda environment in the /project/fsepru/conda_envs/ directory, that way everyone will know it’s available there.
If you are having trouble getting your job to start because the ‘mem’ partition is busy, I recommend trying the “scavenger” partition, it allows you to use the nodes that others have paid for priority access to but aren’t using.
```

### Lit reading:
* Learned that the more resource parameters you specify on Ceres, it will try to make those work rather than doing what works best. So limit your parameters to let Ceres decide how to allocate resources to run job.
* logical core includes hyperthreading
</details>

6. (20Apr2021) slurm job failed twice, with same error message:
```
Traceback (most recent call last):
  File "/project/fsepru/kmou/conda_envs/DRAM/bin/DRAM-setup.py", line 146, in <module>
    args.func(**args_dict)
  File "/project/fsepru/kmou/conda_envs/DRAM/lib/python3.9/site-packages/mag_annotator/database_processing.py", line 457, in prepare_databases
    output_dbs['uniref_db_loc'] = download_and_process_uniref(uniref_loc, temporary, uniref_version=uniref_version,
  File "/project/fsepru/kmou/conda_envs/DRAM/lib/python3.9/site-packages/mag_annotator/database_processing.py", line 110, in download_and_process_uniref
    download_file(uniref_url, uniref_fasta_zipped, verbose=verbose)
  File "/project/fsepru/kmou/conda_envs/DRAM/lib/python3.9/site-packages/mag_annotator/utils.py", line 27, in download_file
    run_process(['wget', '-O', output_file, url], verbose=verbose)
  File "/project/fsepru/kmou/conda_envs/DRAM/lib/python3.9/site-packages/mag_annotator/utils.py", line 38, in run_process
    return subprocess.run(command, check=check, shell=shell, stdout=subprocess.PIPE,
  File "/project/fsepru/kmou/conda_envs/DRAM/lib/python3.9/subprocess.py", line 528, in run
    raise CalledProcessError(retcode, process.args,
subprocess.CalledProcessError: Command '['wget', '-O', 'DRAM_data/database_files/uniref90.fasta.gz', 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz']' returned non-zero exit status 3.
```
* Potential fix: https://github.com/shafferm/DRAM/issues/22
* Try running command separately ...? issue with mmseqs2? I already restarted database processing step. Is UniRef temporarily offline?  Download uniref and give that to `DRAM-setup.py`

7. (20Apr2021) Worked with Chris and found out it was a disk quota exceeded issue. Jules pointed out that the folder I was working out of `/project/fsepru/kmou/conda_envs` was under the ownership of my username instead of `proj-fsepru` and I was hitting my disk quota instead of project quota. Changed ownership of `/project/fsepru/kmou/conda_envs` via: `chown -R --reference=FS9/ conda_envs/`.

8. (20Apr2021) Installed conda environment `DRAM` to `/project/fsepru/conda_envs` with
`conda env create -f environment.yaml --prefix /project/fsepru/conda_envs/DRAM` from my `/project/fsepru/kmou/conda_envs` directory because I already have the `environment.yaml` file there (don't need to do `wget https://raw.githubusercontent.com/shafferm/DRAM/master/environment.yaml`).

9. (20Apr2021) Ran `DRAM.slurm` again, job 5758613. Seems to work. How to know when it stops?

10. (21Apr2021) Database setup still going after 19 hrs. Ran `DRAM-setup.py print_config` and got the following:
```
KEGG db: None
KOfam db: None
KOfam KO list: None
UniRef db: None
Pfam db: None
Pfam hmm dat: None
dbCAN db: None
dbCAN family activities: None
RefSeq Viral db: None
MEROPS peptidase db: None
VOGDB db: None
VOG annotations: None
Description db: None
Genome summary form: /Users/shafferm/lab/DRAM/data/genome_summary_form.tsv
Module step form: /Users/shafferm/lab/DRAM/data/module_step_form.tsv
ETC module database: /Users/shafferm/lab/DRAM/data/etc_module_database.tsv
Function heatmap form: /Users/shafferm/lab/DRAM/data/function_heatmap_form.tsv
AMG database: /Users/shafferm/lab/DRAM/data/amg_database.tsv
```

11. (21Apr2021) I reached out to SCINet staff Jennifer if the reason my job is taking so long is because I didn't specify number of threads. My logic was that I thought if I didn't limit what resources slurm should use, it would be smart to know it has freedom to choose whatever resources are available to run the job faster. Turns out not true. She recommended I make a new folder and submit a new slurm with a larger number of threads. So this was what I did: made `/project/fsepru/kmou/conda_envs/dram2/` and copied over `/project/fsepru/kmou/conda_envs/DRAM.slurm`. Renamed slurm script to `dram2.slurm`. Added 16 threads to the script/. Ran job 5759510.

<details><summary>dram2.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram2                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 16
#SBATCH --mem=550gb
#SBATCH -t 24:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs/DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data2
```
</details>



12. (22Apr2021) Checked progress of job 5759510 and at 22 hours and still not done. Did `salloc` and activated dram conda environment and ran `DRAM-setup.py print_config`. Got the following for dram2 and DRAM_data:
```
KEGG db: None
KOfam db: None
KOfam KO list: None
UniRef db: None
Pfam db: None
Pfam hmm dat: None
dbCAN db: None
dbCAN family activities: None
RefSeq Viral db: None
MEROPS peptidase db: None
VOGDB db: None
VOG annotations: None
Description db: None
Genome summary form: /Users/shafferm/lab/DRAM/data/genome_summary_form.tsv
Module step form: /Users/shafferm/lab/DRAM/data/module_step_form.tsv
ETC module database: /Users/shafferm/lab/DRAM/data/etc_module_database.tsv
Function heatmap form: /Users/shafferm/lab/DRAM/data/function_heatmap_form.tsv
AMG database: /Users/shafferm/lab/DRAM/data/amg_database.tsv
```

13. (22Apr2021) Let the jobs continue running. dram and dram2 job failed because reached time limit. Re-run `dram2.slurm` with the following modified slurm script.

<details><summary>dram2.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram2                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=550gb
#SBATCH -t 96:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs/DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data2
```
I increased the number of cores to 32 (32 cores * 16GB mem per core = 512 GB memory total... should be ok?). Job 5765159.
</details>

14. (23Apr2021) Ran `DRAM-setup.py print_config` in `/project/fsepru/kmou/conda_envs/dram2/` and got same message as yesterday. I asked Chris and showed him my slurm script of the latest job. It turns out I didn't request threads in `DRAM-setup.py` and need to request threads for slurm itself. Chris gave me the corrected slurm script (called it `dram3.slurm`. Ran job 5770291.

<details><summary>dram3.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram3                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=16
#SBATCH --mem=550gb
#SBATCH -t 96:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs/DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data3 --threads 16
```
</details>

15. (26Apr2021) Both dram2 and dram3 jobs finished. Ran `DRAM-setup.py print_config` in `/project/fsepru/kmou/conda_envs/dram2/DRAM_data2` and `/project/fsepru/kmou/conda_envs/DRAM_data3`

<details><summary>dram2 and dram3 stdout details</summary>

* DRAM_data2: Weird that it only lists dram3 files...

```
KEGG db: None
KOfam db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/kofam_profiles.hmm
KOfam KO list: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/kofam_ko_list.tsv
UniRef db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/uniref90.20210423.mmsdb
Pfam db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/pfam.mmspro
Pfam hmm dat: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/Pfam-A.hmm.dat.gz
dbCAN db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/dbCAN-HMMdb-V9.txt
dbCAN family activities: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/CAZyDB.07302020.fam-activities.txt
RefSeq Viral db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/refseq_viral.20210423.mmsdb
MEROPS peptidase db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/peptidases.20210423.mmsdb
VOGDB db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/vog_latest_hmms.txt
VOG annotations: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/vog_annotations_latest.tsv.gz
Description db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/description_db.sqlite
Genome summary form: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/genome_summary_form.20210423.tsv
Module step form: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/module_step_form.20210423.tsv
ETC module database: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/etc_mdoule_database.20210423.tsv
```
  * DRAM_data2 `stdout.5765159.ceres19-mem-4.dram2 `:
  ```
  2021-04-22 15:12:52.471425: Database preparation started
  3:24:12.036707: UniRef database processed
  4:21:15.971722: PFAM database processed
  4:21:20.145036: dbCAN database processed
  4:21:34.134308: RefSeq viral database processed
  4:22:16.861004: MEROPS database processed
  4:25:10.715304: VOGdb database processed
  4:30:19.642823: KOfam database processed
  4:30:23.212796: KOfam ko list processed
  4:30:26.074825: PFAM hmm dat processed
  4:30:26.332474: dbCAN fam activities processed
  4:30:27.158386: VOGdb annotations processed
  4:30:29.653302: DRAM databases and forms downloaded
  4:30:29.715845: Files moved to final destination
  4:30:29.716273: Setting database paths
  4:30:29.719822: Database locations added to CONFIG
  4:30:31.164532: Database connection established
  4:30:31.164567: KEGG descriptions added to description database
  2 days, 9:43:30.210830: UniRef descriptions added to description database
  2 days, 9:43:31.869414: PFAM descriptions added to description database
  2 days, 9:43:32.021867: dbCAN descriptions added to description database
  2 days, 9:43:38.696772: RefSeq viral descriptions added to description database
  2 days, 9:44:01.249890: MEROPS descriptions added to description database
  2 days, 9:44:04.666929: VOGdb descriptions added to description database
  2 days, 9:44:04.666987: Description database populated
  2 days, 9:44:04.669741: Database descriptions updated
  2 days, 9:44:04.729725: Database locations set
  2 days, 9:44:27.444800: Database preparation completed
  ```

* DRAM_data3

```
KEGG db: None
KOfam db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/kofam_profiles.hmm
KOfam KO list: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/kofam_ko_list.tsv
UniRef db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/uniref90.20210423.mmsdb
Pfam db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/pfam.mmspro
Pfam hmm dat: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/Pfam-A.hmm.dat.gz
dbCAN db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/dbCAN-HMMdb-V9.txt
dbCAN family activities: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/CAZyDB.07302020.fam-activities.txt
RefSeq Viral db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/refseq_viral.20210423.mmsdb
MEROPS peptidase db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/peptidases.20210423.mmsdb
VOGDB db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/vog_latest_hmms.txt
VOG annotations: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/vog_annotations_latest.tsv.gz
Description db: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/description_db.sqlite
Genome summary form: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/genome_summary_form.20210423.tsv
Module step form: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/module_step_form.20210423.tsv
ETC module database: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/etc_mdoule_database.20210423.tsv
Function heatmap form: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/function_heatmap_form.20210423.tsv
AMG database: /lustre/project/fsepru/kmou/conda_envs/DRAM_data3/amg_database.20210423.tsv
```
  * DRAM_data3 `stdout.5770291.ceres19-mem-4.dram3` :
  ```
  2021-04-23 13:21:45.932741: Database preparation started
  6:20:31.395180: UniRef database processed
  9:10:33.936828: PFAM database processed
  9:10:38.359894: dbCAN database processed
  9:11:13.738898: RefSeq viral database processed
  9:13:06.458784: MEROPS database processed
  9:16:46.102607: VOGdb database processed
  9:22:20.828967: KOfam database processed
  9:22:24.315550: KOfam ko list processed
  9:22:26.972572: PFAM hmm dat processed
  9:22:27.242745: dbCAN fam activities processed
  9:22:28.067740: VOGdb annotations processed
  9:22:30.818197: DRAM databases and forms downloaded
  9:22:30.885684: Files moved to final destination
  9:22:30.886492: Setting database paths
  9:22:30.890521: Database locations added to CONFIG
  9:22:31.551252: Database connection established
  9:22:31.551289: KEGG descriptions added to description database
  2 days, 10:32:38.774039: UniRef descriptions added to description database
  2 days, 10:32:40.127383: PFAM descriptions added to description database
  2 days, 10:32:40.245665: dbCAN descriptions added to description database
  2 days, 10:32:46.807567: RefSeq viral descriptions added to description database
  2 days, 10:33:08.729547: MEROPS descriptions added to description database
  2 days, 10:33:12.582132: VOGdb descriptions added to description database
  2 days, 10:33:12.582173: Description database populated
  2 days, 10:33:12.584098: Database descriptions updated
  2 days, 10:33:12.651200: Database locations set
  2 days, 10:33:27.046711: Database preparation completed
  ```
I'll stick with DRAM_data3 to run annotation.
</details>

16. (26Apr2021) Ran annotation in `/project/fsepru/kmou/conda_envs/dram3.slurm`, job 5772721.

<details><summary>dram3.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram3                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 1                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH --ntasks-per-core=20
#SBATCH --mem=550gb
#SBATCH -t 96:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs/DRAM
DRAM.py annotate -i '/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/*.fna' -o annotation --threads 20
```
</details>

17. (28Apr2021) Job cancelled (because I accidentally deleted `*pol.fna` files). Rerun `dram3.slurm` on `*pol.fna` files in `/project/fsepru/kmou/FS19C/STECgenomes`, with some changes to `dram3.slurm` script. Job 5801952.

<details><summary>details of why *pol.fna files were deleted, from Section 13c, entry 28Apr2021</summary>

Accidentally deleted all contents of directory `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/`. I thought I wanted to sym link the `*pol.fna` files within `renamed_contigs` in `/project/fsepru/kmou/FS19C/STECgenomes` but then realized I wanted another copy of the `*pol.fna` in this directory. So I deleted the sym link with `rm -r 96commensalEcoli/`. Whoops. I should've done `rm 96commensalEcoli` or use `unlink 96commensalEcoli`. I did manage to save the `*_pol.fna` files from `rename_contigs/` (where the contig names were renamed to work with prokka in gifrop). I copied them over from `/project/fsepru/kmou/FS19C/STECgenomes` to `/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/`. However, everything else is gone. The only things that I'm sad I lost were gifrop results (setting `--prokka_args "--genus Escherichia --species coli`) and the original `*_pol.fna` before I changed the contig names.
</details>

<details><summary>dram3.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=dram3                            # name of the job submitted
  #SBATCH -p mem                                    # name of the queue you are submitting to
  #SBATCH -N 1
  #SBATCH -n 2
  #SBATCH --ntasks-per-core=32
  #SBATCH --mem=550gb
  #SBATCH -t 96:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  #Enter commands here:
  set -e
  set -u
  set +eu

  module load miniconda
  source activate /project/fsepru/kmou/conda_envs/DRAM
  #DRAM-setup.py prepare_databases --output_dir DRAM_data3 --threads 16
  #DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 16 --uniref_loc /path/to/uniref90.fasta.gz
  DRAM.py annotate -i '/project/fsepru/kmou/FS19C/STECgenomes/*.fna' -o annotation_v2 --threads 32
  ```
</details>

18. (3May2021) Job 5801952 failed because it exceeded time limit (96h). I will exclude the Ecoli STEC reference genomes from next job via move `Ecoli_TW14588.fna`, `Ecoli_NADC6564.fna`, `Ecoli_O157H7_EDL933.fna` to `project/fsepru/kmou/FS19C/STECgenomes/repeatstecgenomes/`. Submitted job 5817486.
* I have 133 STEC genomes + 95 commensal E. coli + 3 commensal reference strains = 231 genomes
* DRAM used 20 processors for 80 MAGS and completed in ~17h. I have 231 genomes. If I use 32 cores, with 16 GB per core = 512 GB memory. Set memory at 550GB.

<details><summary>dram3.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=dram3                            # name of the job submitted
  #SBATCH -p mem                                    # name of the queue you are submitting to
  #SBATCH -N 1
  #SBATCH -n 32
  #SBATCH --mem=550gb
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  #Enter commands here:
  set -e
  set -u
  set +eu

  module load miniconda
  source activate /project/fsepru/kmou/conda_envs/DRAM
  DRAM.py annotate -i '/project/fsepru/kmou/FS19C/STECgenomes/*.fna' -o annotation_v3 --threads 32
  ```
</details>

19. (10May2021) Job cancelled because reached time limit. Annotated 213/231 genomes. Need to create a file with all genome file names, and then search from the list which ones are present and not present in directory (total of 18: 5 commensal E. coli and 13 STEC).

20. (11May2021) Manually looked for the ones that are missing. Print lines from `/project/fsepru/kmou/conda_envs/annotation_v3/working_dir` and `/project/fsepru/kmou/FS19C/STEC_genomes` directories and saved as `genomes.txt` and `steccommensalgenomes.txt`, respectively. Compared and printed out missing lines using `comm -13 steccommensalgenomes.txt genomes.txt > needtorun2.txt`, but this didn't work because it also printed lines that were present in both. I did end up finding commensal E. coli files that are missing in `genomes.txt`. As for the STEC genomes, I did `sed '/.fna/d' steccommensalgenomes.txt > steccommensalgenomesnofna.txt` to remove `*.fna` from list. Then opened Excel, copied lines from `steccommensalgenomesnofna.txt` and `genomes.txt`, ran `=IF($A=$D, "Match", "")` to find which rows match and which didn't. Ended up with this list, saved as `secondlist.txt`:
  ```
  25-427FED_pol.fna
  45-429RED_pol.fna
  63-439REC_pol.fna
  75-429FED_pol
  77-430FED_pol
  GCF_001651945.2_ASM165194v2_genomic
  GCF_002208865.2_ASM220886v2_genomic
  GCF_003722195.1_ASM372219v1_genomic
  GCF_013167135.1_ASM1316713v1_genomic
  GCF_013167335.1_ASM1316733v1_genomic
  GCF_013167535.1_ASM1316753v1_genomic
  GCF_013167615.1_ASM1316761v1_genomic
  GCF_013167675.1_ASM1316767v1_genomic
  GCF_017164775.1_ASM1716477v1_genomic
  GCF_017164835.1_ASM1716483v1_genomic
  GCF_017165135.1_ASM1716513v1_genomic
  GCF_017165375.1_ASM1716537v1_genomic
  GCF_017165455.1_ASM1716545v1_genomic
  ```

21. (11May2021) Upload `secondlist.txt` to Ceres. Run `rsync -a /project/fsepru/kmou/FS19C/STECgenomes --files-from=/project/fsepru/kmou/FS19C/STECgenomes/secondlist.txt /project/fsepru/kmou/FS19C/STECgenomes/seconddramannotation` to copy missing `*.fna` files to `seconddramannotation`.

22. (11May2021) Run `dram3.slurm`. Job 5857976

  <details><summary>dram3.slurm script</summary>

    ```
    #!/bin/bash
    #SBATCH --job-name=dram3                            # name of the job submitted
    #SBATCH -p mem                                    # name of the queue you are submitting to
    #SBATCH -N 1
    #SBATCH -n 32
    #SBATCH --mem=550gb
    #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
    #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
    #SBATCH --account fsepru
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=kathy.mou@usda.gov

    #Enter commands here:
    set -e
    set -u
    set +eu

    module load miniconda
    source activate /project/fsepru/kmou/conda_envs/DRAM
    DRAM.py annotate -i '/project/fsepru/kmou/FS19C/STECgenomes/seconddramannotation/*.fna' -o annotation_v4 --threads 32
    ```
  </details>

24. (11May2021) Figured out how to rename `/project/fsepru/kmou/conda_envs/annotation_v3/working_dir/30-440FEC_pol/genes.annotated.gff` to `<directoryName>.gff`. Source: https://unix.stackexchange.com/questions/183021/rename-file-by-its-folder-name. Practiced on example `test` directory with several subdirectories, each containing a `genes.annotated.gff` file (`touch genes.annotated.gff`).
  ```
  #!/bin/bash
  for pathname in */genes.annotated.gff; do
    cp "$pathname" "$( basename "$( dirname "$pathname" )" ).gff"
  done
  ```

25. (11May2021) On Ceres, make `/project/fsepru/kmou/FS19C/STECgenomes/gfftest`. Make above shell script in `gfftest` folder, save as `renamegff.sh` and test run script. It works!

26. () On Ceres, copy over `/project/fsepru/kmou/conda_envs/annotation_v3/working_dir` and `/project/fsepru/kmou/conda_envs/annotation_v4/working_dir` to `/project/fsepru/kmou/FS19C/STECgenomes/gfftest`.

## WGS submission to SRA
* Must complete Biosample entry (which will generate biosample entry in tandem)
* [SRA Quick Start Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submit/) and a more detailed [SRA Submission Guide](https://www.ncbi.nlm.nih.gov/sra/docs/submitbio/)
