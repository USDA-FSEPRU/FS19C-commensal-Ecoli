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

### prokka
Ran prokka without specifying a specific reference genome to use for annotation. I used their default lists of databases for annotation. I did specify Escherichia genus.
#### Files generated (for each isolate):
* *_pol.fasta%.fasta_prokka/
  * *_pol.fasta%.fasta.err
  * *_pol.fasta%.fasta.faa
  * *_pol.fasta%.fasta.ffn
  * *_pol.fasta%.fasta.fna
  * *_pol.fasta%.fasta.fsa
  * *_pol.fasta%.fasta.gbk
  * *_pol.fasta%.fasta.gff
  * *_pol.fasta%.fasta.log
  * *_pol.fasta%.fasta.sqn
  * *_pol.fasta%.fasta.tbl
  * *_pol.fasta%.fasta.tsv
  * *_pol.fasta%.fasta.txt

### roary
#### Files generated:
* accessory.header.embl			
* core_alignment_header.embl
* accessory.tab				
* core_gene_alignment.aln
* accessory_binary_genes.fa		
* gene_presence_absence.Rtab
* accessory_binary_genes.fa.newick
* gene_presence_absence.csv
* accessory_graph.dot			
* number_of_conserved_genes.Rtab
* blast_identity_frequency.Rtab		
* number_of_genes_in_pan_genome.Rtab
* clustered_proteins			
* number_of_new_genes.Rtab
* core_accessory.header.embl		
* number_of_unique_genes.Rtab
* core_accessory.tab			
* pan_genome_reference.fa
* core_accessory_graph.dot		
* summary_statistics.txt

### roary

#### R pheatmap package to create heatmap of select genes and select isolates from `gene_presence_absence.csv`

### gifrop

### PPanGGOLiN

### RAxML
