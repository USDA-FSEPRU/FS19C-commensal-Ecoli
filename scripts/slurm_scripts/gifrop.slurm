#!/bin/bash
#SBATCH --job-name=panpipe                           # name of the job submitted
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



#source activate /project/fsepru/conda_envs/gifrop2
source activate /project/fsepru/conda_envs/gifrop

# you should run this from a folder containing all your assemblies in fasta format with the suffix ".fna"

# change "REFERENCE.GBK" to whatever the name of your desired reference genome is.

# will annotate all with prokka in parallel (will do 24 genomes at a time, each with 1 thread)
# run roary and generate a core genome alignment
# run gifrop and extract classify and cluster genomic islands

# Ecoli K-12 gbk for prokka
#pan_pipe --prokka_args "--proteins Ecoli_K12_MG1655.gbk --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"

# specifying EDL933 gbk for prokka
# pan_pipe --prokka_args "--proteins *.gbk --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"

# without specifying any Ecoli gbk for prokka. See gifropEcoli.slurm on Ceres for this specific command
pan_pipe --prokka_args "--genus Escherichia --species coli --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
