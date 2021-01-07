# Methods

Written summary of methods performed in this repo. This is the methods write up for the paper.

## Sample Collection

## Sequence Analysis

###QC
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
