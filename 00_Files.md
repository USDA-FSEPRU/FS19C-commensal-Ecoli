# Files generated for FS19C project and their location

## (01) Raw data location
### (01A) Commensal E. coli Illumina NovaSeq Reads
**Q:/_TempTransfer/DBayles/mou**
**I:/Mou/CommensalEcoliIsolates1-96NovaSeq**
* `1_33298_01_1-A01-1-428RN3A_HTFVy_1801-CommensalEcoliIsolates1-96.tar`
  * Samples 1-96 from Sept. 15, 2020
* `1_33298_01_1-A01-1-428RN3A_HVHJT_1839-CommensalEcoliIsolates95_96_resequence.tar`
  * Re-sequenced samples 95 and 96 from 16Dec2020 because there were no sequences in the first NovaSeq run

### (01B) Daniel's O157 challenge strain sequence data
**Q:\FSEP_DataTransfer\Trachsel\O157_challenge_strains**
* Includes TW14588, EDL933, FRIK1989, RM6067

## (02) Ceres /project/fsepru/kmou/ master directory setup
```
|_project/
    |_fsepru/
        |_kmou/
            |_programs/
                |_adapt_polish.sh
                |_bbmap/
                |_bin/
                |_exonerate/
                |_gapseq/
                |_good_contig_names.R
                |_share/
                |_SPAdes-3.14.1-Linux/              
            |_conda_envs/
                |_annotation_v3_dramfirstrun/         
                |_annotation_v4_dramsecondrun/
                |_DRAM/
                |_dram3.slurm
                |_DRAM_data3/
                |_environment.yaml
                |_genome_summaries_annotation_v3/
                |_genome_summaries_annotation_v4/
                |_pan_genome_reference.fa
                |_prokka_env/
            |_my_pkg_cache/
            |_from_jules/
            |_for_hannah_fs9/
            |_FS9/
            |_FS19C/
                |_fs19c_sequences/
                  |_assemblywithbbtoolsandspades/
                  |_badsequencedata/
                  |_polished_genomes_100X/
                  |_renamefiles.batch
                  |_touchfilenames.batch
                |_SRAassemblyPipeline.FS19C.SLURM_TEMPLATE  
                |_stecandcommensalEcoli_gifrop/
                    |_gbkforppanggolin_prokka/
                    |_gfftest/
                    |_pan/
                    |_panpipe_logs/
                    |_prokka_cmds.txt
                    |_secondramannotation/
                    |_stecftp.txt
                    |_stec.slurm

```

## (03) Ceres /KEEP/fsepru/kathy.mou/ directory setup
```
|_KEEP/
    |_fsepru/
        |_kathy.mou/
            |_10EcoliIsolates/
```

## (04) Ceres: Sequences, sequence assembly, mash
**/project/fsepru/kmou/FS19C/fs19c_sequences/**
```
|_project/
  |_fsepru/
    |_kmou/
      |_FS19C/
        |_fs19c_sequences/
          |_*.fastq.gz
          |_assemblywithbbtoolsandspades/
          |_badsequencedata/
          |_polished_genomes_100X/
            |_mash/
            |_polishedgenomesforprokka_95isolates6refgenomes/
              |_renamed_contigs
            |_referencegenomes/
          |_renamefiles.batch
          |_touchfilenames.batch
```
* *.fastq.gz sequence data for samples 1-96
* assemblywithbbtoolsandspades/
  * *_covstats.txt
  * *.fasta
  * *.names
  * *_pol.fasta (final polished assembly)
  * *.slurm
  * stderr.*
  * stdout.*
  * *_spades_out/
* badsequencedata/
  * 1-B08-20-427FEC_S27_L001_R*_001.fastq.gz
  * 1-H12-96-441FEC_S103_L001_R*_001.fastq.gz
* mash/
  * *.fasta
  * distances_thirdrun.tab
* polishedgenomesforprokka_95isolates6refgenomes/
  * fs19cpolishedgenomes.tar.gz (all polished assemblies in tar file)  
  * renamed_contigs/
    * *.fna (same as fasta files in assemblywithbbtoolsandspades/ - these are polished)
* referencegenomes/
  * *.fasta (reference genomes used for mash/fastANI, and/or gifrop)

### (05) Ceres: gifrop (genome annotation, pan-genome analysis)
**/project/fsepru/kmou/FS19C/stecandcommensalEcoli_gifrop/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_FS19C/
                |_stecandcommensalEcoli_gifrop/
                    |_*_pol/ or *_genomic/ (for each isolate)
                    |_*.fna
                    |_gbkforppanggolin_prokka/
                    |_gfftest/
                    |_pan/
                    |_panpipe_logs/
                    |_prokka_cmds.txt
                    |_secondramannotation/
                    |_stecftp.txt
                    |_stec.slurm/
```
* gbkforppanggolin_prokka/
  * has output from ppanggolin (used prokka annotated gbk files in ppanggolin)
* pan/
  * has output from roary
* gfftest/
  * Attempt to rename gff files with sample name


### (06) Ceres: DRAM (genome annotation, metabolic pathways analysis)
**/project/fsepru/kmou/conda_envs/**
```
|_project/
    |_fsepru/
        |_kmou/
            |_conda_envs/
                |_annotation_v3_dramfirstrun/      
                    |_working_dir/
                    |_mergedannotation_dramfirstrun.tsv
                |_annotation_v4_dramsecondrun/
                |_DRAM/
                |_dram3.slurm
                |_DRAM_data3/
                |_environment.yaml
                |_genome_summaries_annotation_v3/
                |_genome_summaries_annotation_v4/
                |_pan_genome_reference.fa
                |_prokka_env/
```
* annotation_v3_dramfirstrun/working_dir/
  * 195 of 231 STEC and commensal E. coli annotations generated from `DRAM.py annotate`
  * each isolate's directory has: *.gbk, genes.annotated.faa, genes.annotated.gff, scaffolds.annotated.fa, annotations.tsv, genes.annotated.fna, rrnas.tsv, trnas.tsv
* annotation_v4_dramsecondrun/
  * remaining 18 of 231 STEC and commensal E. coli annotations generated from `DRAM.py annotate`
  * annotations.tsv, genbank, genes.faa, genes.fna, genes.gff, rrnas.tsv, scaffolds.fna, trnas.tsv
* _DRAM_data3/
  * databases compiled by `DRAM-setup.py`
* genome_summaries_annotation_v3/
  * `DRAM.py distill` output of `annotation_v3_dramfirstrun/mergedannotation_dramfirstrun.tsv`
* genome_summaries_annotation_v4/
  * `DRAM.py distill` output of `annotation_v4_dramsecondrun/annotations.tsv`
* _prokka_env/
  * conda environment that has fastani, multiqc, mash prokka, PPanGGOLiN installed

## (07) Files in Files directory
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
