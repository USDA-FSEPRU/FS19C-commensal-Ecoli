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
roary -f ./roary_95isolates_5referencestrains_output -e -n -v -p 16 *.gff
#query_pan_genome  -o pan_genome_results_union -v -a union *.gff
#query_pan_genome  -o pan_genome_results_core -v -a intersection *.gff
#query_pan_genome  -o pan_genome_results_accessory -v -a complement *.gff
#End of file
