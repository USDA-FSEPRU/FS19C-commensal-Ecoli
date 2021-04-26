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

## DRAM database setup, make sure set #SBATCH --ntasks-per-core=16
# DRAM-setup.py prepare_databases --output_dir DRAM_data3 --threads 16

## Use this if manually downloaded uniref90.fasta.gz and want to skip downloading it again when running database setup
# DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 16 --uniref_loc /path/to/uniref90.fasta.gz

DRAM.py annotate -i /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/*.fna -o annotation --threads 20
