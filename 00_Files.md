# Files generated for FS19C project and their location


## Raw data location

### Illumina NovaSeq Reads
**Q:/_TempTransfer/DBayles/mou**
* 1_33298_01_1-A01-1-428RN3A_HTFVy_1801.tar
  * Samples 1-96 from Sept. 15, 2020
* 1_33298_01_1-A01-1-428RN3A_HVHJT_1839.tar
  * Re-sequenced samples 95 and 96 from 16Dec2020 because there were no sequences in the first NovaSeq run

## Data generated from sequence analysis

**(Ceres) /project/fsepru/kmou/FS19C/lane1/**
* *.fastq.gz sequence data for samples 1-96
* badsequencedata/
  * 1-B08-20-427FEC_S27_L001_R*_001.fastq.gz
  * 1-H12-96-441FEC_S103_L001_R*_001.fastq.gz
* fastqc/
* linked/
  * *_covstats.txt
  * *.fasta
  * *.names
  * *_pol.fasta
  * *.slurm
  * stderr.*
  * stdout.*
* renamefiles.batch

### FastQC
* *fastqc.zip
* *fastqc.html

### MultiQC
* FS19all_multiqc_report.html
* FS19all_multiqc_data/
* FS19_1-94_multiqc_report.html
* FS19_1-94_multiqc_data/

### Mash
* distances.tab
* FS19Cmashdistances.xlsx

### FastANI
* fs19cfastanioutput2.out.tab
* FS19CfastANIoutput2.xlsx

### Genome Assembly (BBMap, spades)
**(Ceres)
/project/fsepru/kmou/FS19C/polished_genomes_100X**
* *_pol.fasta
* fs19cfastanioutput.out
* referencegenomes/
  * Ecoli_HS.fasta  
  * Ecoli_K-12_MG1655.fasta  
  * Ecoli_NADC6564.fasta  
  * Ecoli_Nissle1917.fasta  
  * Ecoli_O157H7_EDL933.fasta
* querylist.txt
* referencelist.txt

### Genome Annotation (prokka)
**(Ceres)
/project/fsepru/kmou/FS19C/polished_genomes_100X/polishedgenomesprokka**
* *_pol.fasta
* *_pol.fasta%.fasta_prokka
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

## Files in Files directory
* FS19C_metadata.xlsx
* 1-H12-96-441FEC_S2_L002_R2_001_fastqc.html
* 1-H12-96-441FEC_S2_L002_R1_001_fastqc.html
* 1-H11-95-440FED_S1_L002_R2_001_fastqc.html
* 1-H11-95-440FED_S1_L002_R1_001_fastqc.html
* distances.tab
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
