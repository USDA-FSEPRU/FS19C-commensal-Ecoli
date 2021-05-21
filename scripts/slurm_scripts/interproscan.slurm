#!/bin/bash

#SBATCH --job-name="interproscan"                    # name of the job submitted
#SBATCH -p short                                          # name of the queue you are submitting to
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --mem=60G                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 40:00:00                                    # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=kathy.mou@usda.gov             # will receive an email when job starts, ends or fails
#SBATCH --mail-type=ALL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru

module load interproscan


interproscan.sh -b stec-interproscan -iprlookup --goterms --pathways -i pan_genome_reference.fa -cpu 32
