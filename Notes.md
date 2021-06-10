# Notes, Ideas and Thoughts

# To do
## FS19C
* GO term enrichment of roary and/or DRAM output?
  * Roary output: translate the `pan_genome_reference.fa` fasta file to AA sequence, and then run through interproscan.
* Run `pan_genome_reference.fa` through `DRAM.py distill` but make sure to convert the fasta file to AA sequence (because running `pan_genome_reference.fa` as is produced strange results for `product.html`. I sadly already deleted this output. Could try running it again through `DRAM.py distill` and examining the results closer, compare with `genome_summaries_annotation_v3` results).
* ppangolin: get fasta file of gene families or protein families, find out which ones are present/absent in STEC, commensals - blast those fasta files


## Ceres
* You can change your default Slurm account using running slurm-account-selector.sh on the login node.
  * To see all your Slurm accounts at any time, use “sacctmgr -Pns show user format=account”
