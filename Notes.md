# Notes, Ideas and Thoughts

#Skeleton
### Method
* Summary:
* Began on:
* Completed on:
* Platform:

#To do

## FS19C
* go through gifrop files to find isolates that are hemolysin-, LEE operon negative
  * read papers from Vijay (lee paper, eut operon)
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
* List from Vijay
```
You can search sequenced genomes of E. coli isolates for the genes/operons involved in assimilation of following carbon/sugar substrates :
glucose, sucrose, galactose, arabinose, lactose, fucose, maltose, hexuornate, mannose, ribose, N-acetylglucosamine, N-acetylgalactosamine N-acetylneuraminate, sialic acid and D-gluconate.  E. coli strains use these sugar substrates in various combinations/hierarchy for growth in the animal intestine:
E. coli O157, unlike most human commensal E. coli, can use galactose, hexuronate, mannose and ribose for in vivo growth.

In addition you can screen the sequenced genomes for:
Genes constituting eut operon (for ethanolamine utilization)
Genes encoding bacteriocins and microcins (especially microcins M and H47)
```

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
