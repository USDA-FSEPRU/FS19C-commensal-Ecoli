#!/bin/bash

#SBATCH --job-name=Expt1                          # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 72                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out with %j = job number, %N = node name, %x = job name, %u = username
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints out standard error
#SBATCH --mem=64G   #
#SBATCH --mail-type=ALL
#SBATCH --mail=type=begin   # email when job begins
#SBATCH --mail=type=end     # email when job ends
#SBATCH --mail-user=kathy.mou@usda.gov

# ENTER COMMANDS HERE:
