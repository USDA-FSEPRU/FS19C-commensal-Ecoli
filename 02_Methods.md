# Methods

Written summary of methods performed in this repo. Includes lab notes of how methods performed

## Sample Collection
See OneNote FS19C lab notebook entries xxx
* FS19C Samples 1-96 Final Data.xlsx
* Hannah Sorbitol-positive isolates - MALDI, list for sequencing.xlsx
* Sorbitol-negative isolates - agglutination, MALDI, list for sequencing.xlsx
* KathyMou_NovaSeq_Submission_Form_8June2020.xlsx
* FS19C_metadata.xlsx

## Sequence Analysis
### QC
#### FastQC
* Ran fastqc on FS19C sequence data by generating a slurm script (fastqc.bash) and running on slurm
* Completed on: 1/7/2021
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
* Ran on FS19C sequence data in conda environment (1/7/2021)
* Completed on: 1/7/2021
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
* Ran fastANI on FS19C sequence data by running in conda environment
* Completed on: 12/29/2020
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

5. Moved querylist.txt and referencelist.txt to same directory and ran fastani:
```
fastANI --ql querylist.txt --rl referencelist.txt -o fs19cfastanioutput.out
```

6. Downloaded fs19cfastanioutput.out and made FS19CfastANIoutput.xlsx

* Files generated:
  * FS19CfastANIoutput.xlsx

#### Mash


* Files generated:
  * xxx
