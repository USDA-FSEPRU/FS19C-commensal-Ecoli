# **FS19C project**

This is Kathy Mou's lab notebook describing the comparative genomics project of commensal E. coli isolates and development of a pipeline to identify potential probiotic isolates that could be tested for competitive exclusion against STEC. /n
96 isolates (sorbitol-positive and sorbitol-negative) collected from FS19 STEC colonization study in dairy calves (Fecal and RAMS samples) were tested in this pipeline.

## **Table of contents**
| File Name  | Description |
| -- | -- |
| [00_Files.md](https://github.com/k39ajdM2/Notebook/tree/main/00_Files.md) | List of relevant files used and generated from analysis, locations of files on network drives or Ceres fsepru project folder |
| [00a_Metadata.md](https://github.com/k39ajdM2/Notebook/tree/main/00a_Metadata.md) | Important metadata related to the project and files in 00_Files.md and Files directory|
| [01_Background.md](https://github.com/k39ajdM2/Notebook/tree/main/01_Background.md) | Notes taken from OSQR plan, references listed in OSQR that are related to this project |
| [02_Methods](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods) | Lists entries from OneNote FS19C lab notebook on practicing using tools (bbmap and spades) on 4 FS19C isolates' sequence data; entries from OneNote FS19C lab notebook of relevance, conda environments created for sequence analysis, description and scripts for each method; each individual method as its own step. |
| [02_Methods-Manuscript.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods-Manuscript.md) | Ideas for how to write methods in manuscript |
| [03_Results.md](https://github.com/k39ajdM2/Notebook/tree/main/03_Results.md)| Notes to include in paper |
| [04_Introduction.md](https://github.com/k39ajdM2/Notebook/tree/main/04_Introduction.md) | Notes to include in paper |
| [05_Discussion.md](https://github.com/k39ajdM2/Notebook/tree/main/05_Discussion.md) | Notes to include in paper |
| [06_AuthorInfo.md](https://github.com/k39ajdM2/Notebook/tree/main/06_AuthorInfo.md) | List of authors, contact info, contributions |
| [BashScriptLesson.md](https://github.com/k39ajdM2/Notebook/tree/main/BashScriptLesson.md) | Contains the script for `rename.contigs` and the for-loop script `rename_contigs.sh` to run `rename.contigs` on all fasta files in a directory. This page also explains each step of `rename_contigs` and `rename_contigs.sh` |
| [Notes.md](https://github.com/k39ajdM2/Notebook/tree/main/Notes.md) | To-do list and other miscellaneous items |
| [Files](https://github.com/k39ajdM2/Notebook/tree/main/Files) | Relevant data files or files in scripts |
| [scripts](https://github.com/k39ajdM2/Notebook/tree/main/scripts) | Scripts used in data analysis on HPC (Ceres) or locally (RStudio) |

## **General Pipeline**
1. Assemble short read sequences (Illumina) with BBMap and spades
2. Fetch STEC genomes from NCBI (RefSeq)
3. Determine how closely related isolates are using mash and MDS (in R)
4. Run gifrop to annotate with prokka and analyze pan-genomes of commensal E. coli and STEC with roary.
5. Identify metabolic pathways of interest in STECs that are also present in commensal E. coli (DRAM, BioCyc, blast) to find candidate commensal E. coli isolates.
6. Interproscan: GO term enrichment using roary output
