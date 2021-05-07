# **FS19C project**

This is a lab notebook describing the comparative genomics project of commensal E. coli isolates and development of a pipeline to identify potential probiotic isolates that could be tested for competitive exclusion against STEC. /n
96 isolates (sorbitol-positive and sorbitol-negative) collected from FS19 STEC colonization study in dairy calves (Fecal and RAMS samples) were tested in this pipeline.

## **Table of contents**
| File Name  | Description |
| -- | -- |
| [00_Files.md](https://github.com/k39ajdM2/Notebook/tree/main/00_Files.md) | List of relevant files used and generated from analysis |
| [00a_Metadata.md](https://github.com/k39ajdM2/Notebook/tree/main/00a_Metadata.md) | Important metadata related to the project and files in 00_Files.md and Files directory|
| [01_Background.md](https://github.com/k39ajdM2/Notebook/tree/main/01_Background.md) | Notes taken from OSQR plan, references listed in OSQR that are related to this project |
| [02_Methods](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods) | Lists entries from OneNote FS19C lab notebook on practicing using tools (bbmap and spades) on 4 FS19C isolates' sequence data; entries from OneNote FS19C lab notebook of relevance, conda environments created for sequence analysis, description and scripts for each method; each individual method as its own step; summary of sequence analyses methods performed (will eventually remove MethodsSummary.md) |
| [02_Methods-Manuscript.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods-Manuscript.md) | Ideas for how to write methods in manuscript |
| [03_Results.md](https://github.com/k39ajdM2/Notebook/tree/main/03_Results.md)| Notes to include in paper |
| [04_Introduction.md](https://github.com/k39ajdM2/Notebook/tree/main/04_Introduction.md) | Notes to include in paper |
| [05_Discussion.md](https://github.com/k39ajdM2/Notebook/tree/main/05_Discussion.md) | Notes to include in paper |
| [06_AuthorInfo.md](https://github.com/k39ajdM2/Notebook/tree/main/06_AuthorInfo.md) | List of authors, contact info, contribution |
| [BashScriptLesson.md](https://github.com/k39ajdM2/Notebook/tree/main/BashScriptLesson.md) | Contains the script for `rename.contigs` and the for-loop script `rename_contigs.sh` to run `rename.contigs` on all fasta files in a directory. This page also explains each step of `rename_contigs` and `rename_contigs.sh` |
| [Notes.md](https://github.com/k39ajdM2/Notebook/tree/main/Notes.md) | To-do list |
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
