# **Methods_README**
1. These are the finalized version of each step in the pipeline.
2. Notes on what order to run some of these methods after steps 1-3 (01_SequenceAssembly, 02_FetchEcoliGenomes, 03_QC):
  * Steps 4 (04_Pan-genomeAnalysis_gifrop) and 5 (05_DRAM) can be run simultaneously and are not dependent on each other.
  * Step 6 (06_Interproscan) can be run only after completing step 4 or 5.
  * Step 7 (07_ppanggolin) must be run only after step 4 is completed.

## **Table of contents**
| File Name | Description |
| -- | -- |
| [01_SequenceAssembly.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/01_SequenceAssembly.md) | Assembly sequences with bbmap and spades. |
| [02_FetchEcoliGenomes.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/02_FetchEcoliGenomes.md) | Import STEC and non-STEC E. coli genomes from RefSeq. |
| [03_QC.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/03_QC.md) | Run Mash and MDS in R. |
| [04_Pan-genomeAnalysis_gifrop.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/04_Pan-genomeAnalysis_gifrop.md) | Annotation, pangenome analysis with gifrop |
| [05_DRAM.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/05_DRAM.md) | Database setup, annotation, distill to get biological pathways present/absent in each isolate |
| [06_Interproscan.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/06_Interproscan.md) | GO term enrichment of pan-genome results from roary, DRAM (?) |
| [07_ppanggolin.md](https://github.com/k39ajdM2/Notebook/tree/main/02_Methods/07_ppanggolin.md) | Another method of pan-genome analysis |
