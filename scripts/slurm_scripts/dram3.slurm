#!/bin/bash
#SBATCH --job-name=dram3                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 32                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH --ntasks-per-core=20            # didn't run this step for `DRAM.py annotate`
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

## DRAM database setup, make sure set #SBATCH --ntasks-per-core=16
# DRAM-setup.py prepare_databases --output_dir DRAM_data3 --threads 16

## From Chris Anderson: Use this if manually downloaded uniref90.fasta.gz and want to skip downloading it again when running database setup
# DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 16 --uniref_loc /path/to/uniref90.fasta.gz

## DRAM annotate
# DRAM.py annotate -i '/project/fsepru/kmou/FS19C/stecandcommensalEcoli_gifrop/seconddramannotation/*.fna' -o annotation_v4 --threads 32

## DRAM distill for annotation_v4
# DRAM.py distill -i annotation_v4/annotations.tsv -o genome_summaries_annotation_v4 --trna_path annotation_v4/trnas.tsv --rrna_path annotation_v4/rrnas.tsv

## DRAM distill for annotation_v3
DRAM.py distill -i annotation_v3_dramfirstrun/mergedannotation_dramfirstrun.tsv -o genome_summaries_annotation_v3
