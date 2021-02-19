# Files generated for FS19C project and their location

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
                    |_polishedgenomesprokka_95isolates5refgenomes/

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
* FS19Cmashdistances.xlsx

### FastANI-related files of importance
* fs19cfastanioutput2.out.tab
* FS19CfastANIoutput2.xlsx

### Genome Assembly (BBMap, spades, mash, fastani)
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_fs19cfastanioutput.out2
                    |_referencegenomes/
                    |_polishedgenomesprokka_95isolates5refgenomes/
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

### Genome Annotation, Pan-genome, Genomic Island ID
**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka_95isolates5refgenomes/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_polishedgenomesprokka_95isolates5refgenomes/
                        |_*_pol.fasta
                        |_*_pol.fasta%.fasta_prokka/
```
* *_pol.fna
* *_pol/
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
Ecoli_HS.gff
Ecoli_K-12_MG1655.gff
Ecoli_NADC6564.gff
Ecoli_Nissle1917.gff
Ecoli_O157H7_EDL933.gff
Ecoli_TW14588.gff
gene_presence_absence.csv
gene_presence_absence.Rtab
_inflated_mcl_groups
_inflated_unsplit_mcl_groups
_labeled_mcl_groups
number_of_conserved_genes.Rtab
number_of_genes_in_pan_genome.Rtab
number_of_new_genes.Rtab
number_of_unique_genes.Rtab
P7mObRUGMX
pan_genome_reference.fa

summary_statistics.txt
_uninflated_mcl_groups
  ```
  * pan_genome_sequences/
    * *.aln
  * gifrop_out/
    ```
    clustered_island_info.csv  my_islands
    figures                    pan_only_islands.csv
    gifrop.log                 pan_with_island_info.csv
    islands_pangenome_gff.csv  sequence_data
    ```
* panpipe_logs/
  * gifrop.log  
  * prokka_logs.txt  
  * roary.log
* prokka_cmds.txt

**(Ceres) /project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka/prokka_gbk**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_polished_genomes_100X/
                    |_polishedgenomesprokka_95isolates5refgenomes
                        |_prokka_gbk/
                            |_*_pol.fasta%.fasta.gbk
                            |_ppanggolin_output_DATE/
```
* *_pol.fasta%.fasta.gbk
* Ecoligbkpath.txt
* Ecoligbk.txt
* ppanggolin_output_DAT
  * gene_presence_absence.Rtab       
  * organisms_statistics.tsv  
  * pangenomeGraph_light.gexf  
  * projection/
  * matrix.csv                       
  * pangenomeGraph.gexf       
  * pangenome.h5               
  * tile_plot.html
  * mean_persistent_duplication.tsv  
  * pangenomeGraph.json       
  * partitions/
  * Ushaped_plot.html

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
* FS19C Samples 1-96 Final Data.xlsx
* FS19C_metadata.xlsx
* fs19cfastanioutput2.out.tab
* FS19CfastANIoutput2.xlsx
* FS19Cmashdistances.xlsx
* Hannah Sorbitol-positive isolates - MALDI, list for sequencing.xlsx
* KathyMou_NovaSeq_Submission_Form_8June2020.xlsx
* Sorbitol-negative isolates - agglutination, MALDI, list for sequencing.xlsx
* FS19C 96 S-S+ E. coli gDNA gels.pdf




## Files I need to add in Files directory
* querylist.txt (genome assembly)
* referencelist.txt (genome assembly)
* polished genomes as a tar file, place on Ceres fsepru directory somewhere easy to find
* FS19C slurm progress.xlsx
