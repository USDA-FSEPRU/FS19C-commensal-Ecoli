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
