# Methods

Written summary of methods performed in this repo. Includes lab notes of how methods performed

## Sample Collection
Need to do

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
* Why are the plots flat plot (not interactive)?
  * From: https://multiqc.info/docs/
  "Flat plots
  Reports with large numbers of samples may contain flat plots. These are rendered when the MultiQC report is generated using MatPlotLib and are non-interactive (flat) images within the report. The reason for generating these is that large sample numbers can make MultiQC reports very data-intensive and unresponsive (crashing people's browsers in extreme cases). Plotting data in flat images is scalable to any number of samples, however.
  Flat plots in MultiQC have been designed to look as similar to their interactive versions as possible. They are also copied to multiqc_data/multiqc_plots"
* I noticed the FS19all_multiqc_report.html report showed samples 95 and 96 having huge number of reads and per sequence GC content had several large peaks.
* To see if the discrepancies are due to samples 95 and 96 having so many reads (as a result of re-sequencing them) and are therefore skewing the report stats, I moved samples 95 and 96 fastqc.zip to Fastqc_Sample95_96 directory, and ran multiqc to generate FS19_1-94_multiqc_report.html and FS19_1-94_data directory. 
* I also looked at the fastqc reports for samples 95 and 96 individually.
* Files generated:
 * FS19all_multiqc_report.html
 * FS19all_multiqc_data directory
 * FS19_1-94_multiqc_report.html
 * FS19_1-94_multiqc_data directory
