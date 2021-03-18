# Results
Written summary of results obtained in this repo. This is the results write up for the paper and includes figures.

## 02a_Methods.md
After testing bbmap and spades with 4 of the E. coli isolates (#1, #20, #94, #96), I decided to run assemblies at 100X coverage, skip isolate #20. Modified slurm file from Jules and renamed as `SRAassemblyPipeline.FS19C.SLURM_TEMPLATE`.

## 02_Methods.md

### BBmap + spades
Assembly and qc check was successful for all isolates except for #20 (skip).
##### Files/directories generated (for each isolate if indicated with a '*'):
* polishedfasta.txt
* *_1.fastq.gz
* *_2.fastq.gz
* *_pol.fasta
* *.slurm
* *.names
* *_covstats.txt
* *.fasta
* *_spades_out/

### fastQC
Ran all samples through fastqc to use for multiqc. I also looked at fastqc reports for re-sequenced isolates #95 and #96 - they look good.
##### Files generated (for each isolate):
  * *fastqc.zip
  * *fastqc.html

### MultiQC
Ran multiqc for isolates #1-94 and report shows sequences are good. Looked at fastqc reports for re-sequenced isolates #95 and #96 and confirmed they are good.
##### Files generated:
 * FS19all_multiqc_report.html
 * FS19all_multiqc_data directory
 * FS19_1-94_multiqc_report.html
 * FS19_1-94_multiqc_data directory
 * 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
 * 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
 * 1-H11-95-440FED_S1_L002_R1_001_fastqc.html

### MDS from fastANI and mash
#### fastANI
##### Files generated:
  * FS19CfastANIoutput.xlsx
  * FS19CfastANIoutput2.xlsx
  * fs19cfastanioutput.out
  * fs19cfastanioutput2.out

#### mash
##### Files generated:
  * distances_thirdrun.tab
  * FS19Cmashdistances.xlsx

#### MDS
There weren't any obvious outliers, like what Jules had hinted when he ran the code on his end. At microbe meeting (22Jan2021), Crystal pointed out that it's good we're seeing some differences of the commensal isolates (like isolates collected from EDL933 group cluster separately from EDL933 isolate) from the STEC isolates in this plot. Regardless of whether we do find any differences in metabolic genes or not, there are other differences we can explore too.

Also ran MDS of mash-generated `distances_thirdrun.tab` of 95 isolates + 6 reference strains + 3 non-*E. coli* isolates. The 3 non-*E. coli* isolates were either in the same family, different family, or non-proteobacteria: Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, Campylobacter jejuni subsp. jejuni NCTC 11168, Clostridium saccharoperbutylacetonicum N1-4(HMT), respectively. We got what we expected where non-*E. coli* isolates were very distant from *E. coli* isolates and reference strains.

##### Files generated:
* fastANImashMDSheatmaps.pptx
* FS19C_fastaniMDS.tiff
* FS19C_mashMDS.tiff
* qc_mds.R

### gifrop
##### Files generated:
* **_pol/ or Ecoli_*/
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

### PPanGGOLiN
##### Files generated:
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

### RAxML
##### Files generated:
* RAxML_bestTree.core_genome_tree_1
* RAxML_bipartitionsBranchLabels.core_genome_tree_1
* RAxML_bipartitions.core_genome_tree_1
* RAxML_bootstrap.core_genome_tree_1
* RAxML_info.core_genome_tree_1

### gapseq
