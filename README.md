# **FS19C project**

This is a lab notebook describing the comparative genomics project of commensal E. coli isolates and development of a pipeline to identify potential probiotic isolates that could be tested for competitive exclusion against STEC. \n
96 isolates (sorbitol-positive and sorbitol-negative) collected from FS19 STEC colonization study in dairy calves (Fecal and RAMS samples) were tested in this pipeline.

## **Table of contents**
| File Name  | Description |
| -- | -- |
| 00_Files.md | List of relevant files used and generated from analysis |
| 00a_Metadata.md | Important metadata related to the project and files in 00_Files.md and Files directory|
| 01_Background.md | Notes taken from OSQR plan, references listed in OSQR that are related to this project |
| 02a_Pre-Methods.md | Lists entries from OneNote FS19C lab notebook on practicing using tools (bbmap and spades) on 4 FS19C isolates' sequence data |
| 02b_Methods.md | Lists entries from OneNote FS19C lab notebook of relevance, conda environments created for sequence analysis, description and scripts for each method|
| 02c_MethodsSummary.md| Summary of sequence analyses methods performed. Only shows the final corrected code + procedures for each of the tools used.|
| 02d_MethodsDescription.md | Ideas for how to write methods in manuscript |
| 03_Results.md | Notes to include in paper |
| 04_Introduction.md | Notes to include in paper |
| 05_Discussion.md | Notes to include in paper |
| 06_AuthorInfo.md | List of authors, contact info, contribution |
| BashScriptLesson.md | Contains the script for `rename.contigs` and the for-loop script `rename_contigs.sh` to run `rename.contigs` on all fasta files in a directory. This page also explains each step of `rename_contigs` and `rename_contigs.sh` |
| Notes.md | To-do list |
| [Files](https://github.com/k39ajdM2/Notebook/tree/main/Files) | Relevant data files or files in scripts |
| [scripts](https://github.com/k39ajdM2/Notebook/tree/main/scripts) | Scripts used in data analysis on HPC (Ceres) or locally (RStudio) |

## **General Pipeline**
1. Assemble short read sequences (Illumina) with BBMap and spades
2. Determine how closely related isolates are using mash and MDS (in R)
3. Run gifrop to annotate with prokka, pan-genome with roary, and run through several databases to identify genomic islands (virulence gene focus)
4. Use the pangenome matrix file produced by roary (`gene_presence_absence.Rtab`)
5. Group genomes to whatever groups you are interested in comparing.
6. Run Fisher's exact test with FDR correction to see which genes are more significantly enriched in each group.
7. Identify isolates with interested metabolic pathways (DRAM, BioCyc, blast).
8. Of those isolates, exclude those that have virulence genes (blast).
