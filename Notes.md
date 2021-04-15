# Notes, Ideas and Thoughts

#Skeleton
### Method
* Summary:
* Began on:
* Completed on:
* Platform:

#To do

## FS19C
* new approach:
    * run roary of commensals and STEC separate but make sure they have the same annotation approach?
    * How to do Fisher exact test?
    * DRAM
* current approach:
    * run prokka with Escherichia coli of commensals + STECs via gifrop
    * parallel blast to make sure virulence genes exist in the strains
      * custom database of virulence genes: track which source (EDL933, Sakai, etc.)


* read papers from Vijay (eut operon)
* Read DRAM publication, go through github site to see how to download
* highlight isolates in raxml tree once can narrow down which isolates have LEE operon, other virulence genes. Maybe gifrop results can help narrow that down?
* check out anvi'o: https://merenlab.org/2016/11/08/pangenomics-v2/
* troubleshoot gapseq results
* update 02b_Methods.md (gifrop, gapseq), 03_Results.md (gifrop, gapseq), 02c_MethodsSummary.md (gapseq)
* Run isolates through Daniel's github tools
  * https://github.com/nielsend/O157LineageAssignment
    `./LSPA6Long.sh [output CSV File Name] [Fasta File(s)]`
  * https://github.com/nielsend/ABRicateSequenceExtraction
    * Not as sure how applicable the output from this will be... pull out specific genes of interest and search for polymorphisms, etc...
* lit review: read https://mbio.asm.org/content/3/3/e00050-12.long


## Ceres
* You can change your default Slurm account using running slurm-account-selector.sh on the login node.
  * To see all your Slurm accounts at any time, use “sacctmgr -Pns show user format=account”
* reduce space used on Ceres: /project/fsepru/kmou/FS19C/lane1/**
  * *.fastq.gz sequence data for samples 1-96
  * tar these files once finish running assembly slurm scripts and get covstats.txt
  * /project/fsepru/kmou/FS19C/polished_genomes: gzip fasta files from polished_genomes_100X, mash_all/, polishedgenomesprokka/, referencegenomes/
  * remove intermediate files
  * make sure raw data somehwere!

# Journal Target
* [Journal submission guidelines](link)

# Random stuff
<details><summary>ClIcK iF yOu DaRe</summary>
Boo!! Congratulations
</details>
