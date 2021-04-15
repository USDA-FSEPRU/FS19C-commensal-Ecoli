polishedgenomesprokka_95isolates6refgenomes/# Files generated for FS19C project and their location

## Raw data location
### Illumina NovaSeq Reads
**Q:/_TempTransfer/DBayles/mou**
**I:/Mou**
* 1_33298_01_1-A01-1-428RN3A_HTFVy_1801.tar
  * Samples 1-96 from Sept. 15, 2020
* 1_33298_01_1-A01-1-428RN3A_HVHJT_1839.tar
  * Re-sequenced samples 95 and 96 from 16Dec2020 because there were no sequences in the first NovaSeq run

## Ceres /project/fsepru/kmou/ directory setup
```
|_project/
    |_fsepru/
        |_kmou/
            |_programs/
              |_bin  
              |_exonerate  
              |_gapseq
              |_share
            |_dot_files/
                |_.conda
                |_software/
                  |_adapt_polish.sh
                  |_bbmap/
                  |_good_contig_names.R
                  |_SPAdes-3.14.1-Linux/
            |_for_hannah_fs9/
            |_FS9/
            |_FS9_R/
            |_FS19C/
                |_lane1/
                |_polished_genomes_100X/
                    |_polishedgenomesforprokka_95isolates6refgenomes/

```

## Data generated from sequence analysis
**(Ceres) /project/fsepru/kmou/FS19C/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_FS19C_4Samples100X/
                |_FS19C_4Samples250X/
```
* FS19C_4Samples100X/
  * Files for isolates 1, 20, 94, and 96
* FS19C_4Samples250X
  * Files for isolates 1, 20, 94, and 96
* 1_33298_01_1-A01-1-428RN3A_HTFVy_1801.tar

**(Ceres) /project/fsepru/kmou/FS19C/lane1/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_lane1/
                    |_*.fastq.gz
                    |_badsequencedata/
                    |_fastqc/
                    |_linked/
```
* *.fastq.gz sequence data for samples 1-96
* badsequencedata/
  * 1-B08-20-427FEC_S27_L001_R*_001.fastq.gz
  * 1-H12-96-441FEC_S103_L001_R*_001.fastq.gz
* fastqc/
  * *fastqc.zip
  * *fastqc.html
  * Fastqc_Samples95_96/
    * *.html
    * *.zip
  * FS19_1-94_multiqc_data/
    * *.txt
    * *.log
    * *.jason
  * FS19_1-94_multiqc_report.html
  * FS19all_multiqc_data/
    * *.txt
    * *.log
    * *.jason
  * FS19all_multiqc_report.html
  * stderr.*.fastqc
  * stdout.*.fastqc
  * stderr.*.multiqc
  * stdout.*.multiqc
* linked/
  * *_covstats.txt
  * *.fasta
  * *.names
  * *_pol.fasta
  * *.slurm
  * stderr.*
  * stdout.*

### FastQC-related files of importance
* *fastqc.zip
* *fastqc.html

### MultiQC-related files of importance
* FS19all_multiqc_report.html
* FS19all_multiqc_data/
* FS19_1-94_multiqc_report.html
* FS19_1-94_multiqc_data/

### Mash-related files of importance
* distances_thirdrun.tab

### FastANI-related files of importance
* fs19cfastanioutput2.out.tab
* FS19CfastANIoutput2.xlsx

### Genome Assembly (BBMap, spades, mash, fastani)
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_fs19cfastanioutput.out2
                    |_referencegenomes/
                    |_polishedgenomesforprokka_95isolates6refgenomes/
                    |_querylist2.txt
                    |_mash_all/
```
* fs19cfastanioutput.out2
* referencegenomes/
  * Ecoli_HS.fasta  
  * Ecoli_K-12_MG1655.fasta  
  * Ecoli_NADC6564.fasta  
  * Ecoli_Nissle1917.fasta  
  * Ecoli_O157H7_EDL933.fasta
  * Ecoli_TW14588.fasta
  * Cjejuni_11168.fasta
  * Clostridium_N1-4.fasta
  * Styphimurium_LT2.fasta
* querylist2.txt
* mash_all/
  * *_pol.fasta
  * *.fasta
  * distances_1strun.tab
  * distances_secondrun.tab
  * distances_thirdrun.tab

### Genome Annotation, Pan-genome Analysis, Phylogenetic Tree, Genomic Island ID (prokka, roary, raxml, gifrop)
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_polishedgenomesforprokka_95isolates6refgenomes/
                        |_*.fna
                        |_renamed_contigs/
```
* *.fna
* */
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
  * proteins.*
    * .faa
    * .pdb
    * .pot
    * .ptf
    * .pto
* pan/
  * *.gff
  ```
  accessory_binary_genes.fa
accessory_binary_genes.fa.newick
_accessory_clusters
_accessory_clusters.clstr
accessory_graph.dot
accessory.header.embl
accessory.tab
blast_identity_frequency.Rtab
_blast_results
_clustered
_clustered.clstr
clustered_proteins
_combined_files
_combined_files.groups
core_accessory_graph.dot
core_accessory.header.embl
core_accessory.tab
core_alignment_header.embl
core_gene_alignment.aln
gene_presence_absence.csv
gene_presence_absence.Rtab
_inflated_mcl_groups
_inflated_unsplit_mcl_groups
_labeled_mcl_groups
number_of_conserved_genes.Rtab
number_of_genes_in_pan_genome.Rtab
number_of_new_genes.Rtab
number_of_unique_genes.Rtab
pan_genome_reference.fa
RAxML_bestTree.core_genome_tree_1
RAxML_bipartitionsBranchLabels.core_genome_tree_1
RAxML_bipartitions.core_genome_tree_1
RAxML_bootstrap.core_genome_tree_1
RAxML_info.core_genome_tree_1
summary_statistics.txt
_uninflated_mcl_groups
  ```
  * pan_genome_sequences/
    * *.aln
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
* panpipe_logs/
  * gifrop.log  
  * prokka_logs.txt  
  * roary.log
* prokka_cmds.txt

### Pan-genome Analysis with PPanGGOLiN
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/renamed_contigs/prokka_gbk/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_polishedgenomesforprokka_95isolates6refgenomes/
                        |_renamed_contigs/
                            |_prokka_gbk/
                                |_*_pol.gbk
                                |_ppanggolin_output_DATE2021-02-24_HOUR09.59.47_PID10989/
```
* *.gbk
* Ecoligbkpath.txt
* Ecoligbklist.txt
* Ecoligbk.txt
* ppanggolin_output_DATE2021-02-24_HOUR09.59.47_PID10989/
  * gene_presence_absence.Rtab       
  * organisms_statistics.tsv  
  * pangenomeGraph_light.gexf  
  * projection/
    * *.gbk.tsv
  * matrix.csv                       
  * pangenomeGraph.gexf       
  * pangenome.h5               
  * tile_plot.html
  * mean_persistent_duplication.tsv  
  * pangenomeGraph.json       
  * partitions/
    ```
    cloud.txt            persistent.txt  S3.txt  shell.txt           undefined.txt
    exact_accessory.txt  S1.txt          S4.txt  soft_accessory.txt
    exact_core.txt       S2.txt          S5.txt  soft_core.txt
    ```
  * Ushaped_plot.html


### Pan-genome Analysis with gapseq
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesforprokka_95isolates6refgenomes/gapseq**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_polishedgenomesforprokka_95isolates6refgenomes/
                       |_gapseq/
```
* -all-*Pathways.tbl
* -all-*Reactions.tbl
* stderr.*.gapseq
* stdout.*.gapseq







## Files in Files directory
* FS19C_metadata.xlsx
* 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
* 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html
* 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
* 1-H11-95-440FED_S1_L002_R1_001_fastqc.html
* distances_thirdrun.tab
* FS19_1-94_multiqc_report.html
* FS19_outline_09.23.19.docx
* FS19all_multiqc_report.html
* FS19C 96 S-S+ E. coli gDNA gels.pdf
* FS19C Samples 1-96 Final Data.xlsx
* FS19C_metadata.xlsx
* fs19cfastanioutput2.out.tab
* FS19CfastANIoutput2.xlsx
* gene_presence_absence.Rtab
* Hannah Sorbitol-positive isolates - MALDI, list for sequencing.xlsx
* KathyMou_NovaSeq_Submission_Form_8June2020.xlsx
* Sorbitol-negative isolates - agglutination, MALDI, list for sequencing.xlsx






## Files I need to add in Files directory
* querylist.txt (genome assembly)
* referencelist.txt (genome assembly)
* polished genomes as a tar file, place on Ceres fsepru directory somewhere easy to find
* FS19C slurm progress.xlsx
