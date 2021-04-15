#!/bin/bash

#SBATCH --job-name=virgene_BLAST                          # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 72                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name, %x adds job name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --mem=64G   #
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov

# ENTER COMMANDS HERE:

module load blast+
module load parallel


#mkdir metal_res

#find ./passing_QC -name "*fasta" | parallel 'blastn -db ~/blastdb/METAL_ISLAND.fasta -query {} -out ./metal_res/{/.}.metal -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore"'

#cd metal_res

# collects and combines all blast results into one file
#find . -name "*metal" | xargs -n 100 cat >> ALL_metal.results

# removes all files ending with 'metal'
#find . -name "*metal" | xargs rm

#End of file

mkdir virulencegenes

find /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates6refgenomes/renamed_contigs/gifropWithEcoliAnnotation/ -name "*.fna" | parallel 'blastn -db ~/blastdb/METAL_ISLAND.fasta -query {} -out ./virulencegenes/{/.}.virulence -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore"'

cd virulencegenes

# collects and combines all blast results into one file
find . -name "*virulence" | xargs -n 100 cat >> ALL_virulence.results

# removes all files ending with 'metal'
find . -name "*virulence" | xargs rm

#End of file
