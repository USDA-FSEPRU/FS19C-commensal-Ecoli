# 04_Pan-genome Analysis
* Summary: Run gifrop to annotate genomes with Prokka and do pan-genome analysis with Roary.
* Platform: Ceres
  * `/project/fsepru/kmou/FS19C/STECgenomes/`

1. For fasta files, you'll need to rename suffix of `.fasta` files to `.fna`. Also need to modify contig IDs in fasta files for gifrop to work properly (able to call virulence genes, etc. from various databases). Jules showed me his `rename_contigs` script

  <details><summary>rename_contigs.sh</summary>

  ```
  rename_contigs() {
          FILE=$1
          BASE=$2
          awk -v basev="$BASE" '/^>/{print ">"basev"_"++i;next}{print}' "$FILE" > "$BASE"_rename.fasta
          rm $FILE
          mv "$BASE"_rename.fasta $FILE
  }
  export -f rename_contigs
  ```
  </details>

2. Add `rename_contigs.sh` script to your `.bashrc` profile in home directory on Ceres. Do `source ~/.bashrc` from wherever you are to update bashrc profile.

3. Run the bash script `~/scripts/renamecontigs_forloop.sh` in desired directory where fasta files are. This script runs a for-loop of `rename_contigs` on all fasta files and also adds `_pol` to name of file. See `BashScriptLesson.md` for more details.

4. Run `gifrop.slurm`

  <details><summary>gifrop.slurm script</summary>

  ```
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

  source activate /project/fsepru/conda_envs/gifrop

  # without specifying any Ecoli gbk for prokka.
  pan_pipe --prokka_args "--genus Escherichia --species coli --cpus 1 --centre X --compliant" --roary_args "-p 24 -e -n -z -v" --gifrop_args "--threads 24"
  ```
  </details>

#### Files generated:
* *_pol/ or Ecoli_*/
  * *_pol.err
  * *_pol.faa
  * *_pol.ffn
  * *_pol.fna
  * *_pol.fsa
  * *_pol.gbk
  * *_pol.gff
  * *_pol.log
  * *_pol.sqn
  * *_pol.tbl
  * *_pol.tsv
  * *_pol.txt
  * proteins.faa
  * proteins.pdb
  * proteins.pot
  * proteins.ptf
  * proteins.pto
* pan/
  * *.gff
  * accessory_binary_genes.fa
  * accessory_binary_genes.fa.newick
  * _accessory_clusters
  * _accessory_clusters.clstr
  * accessory_graph.dot
  * accessory.header.embl
  * accessory.tab
  * blast_identity_frequency.Rtab
  * _blast_results
  * _clustered
  * _clustered.clstr
  * clustered_proteins
  * _combined_files
  * _combined_files.groups
  * core_accessory_graph.dot
  * core_accessory.header.embl
  * core_accessory.tab
  * core_alignment_header.embl
  * core_gene_alignment.aln
  * core_gene_alignment.aln.reduced
  * gene_presence_absence.csv
  * gene_presence_absence.Rtab
  * gifrop_out/
    * clustered_island_info.csv
    * figures/
      * island_length_histogram.png
      * islands_per_isolate_no_unknowns.png
      * islands_per_isolate.png
      * Number_of_occurances.png
      * Number_of_occurances_secondary.png
    * gifrop.log
    * islands_pangenome_gff.csv
    * my_islands/
      * abricate/
        * All_islands.megares2
        * All_islands.ncbi
        * All_islands.plasmidfinder
        * All_islands.vfdb
        * All_islands.viroseqs
      * island_info.csv
      * All_islands.fasta
    * pan_only_islands.csv
    * pan_with_island_info.csv
    * sequence_data/
      * *.fna
      * *_short.gff
  * _inflated_mcl_groups
  * _inflated_unsplit_mcl_groups
  * _labeled_mcl_groups
  * M7lUUryBzC/
    * *.gff.proteome.faa
  * number_of_conserved_genes.Rtab
  * number_of_genes_in_pan_genome.Rtab
  * number_of_new_genes.Rtab
  * number_of_unique_genes.Rtab
  * pan_genome_reference.fa
  * pan_genome_sequences/
  * summary_statistics.txt
  * _uninflated_mcl_groups
