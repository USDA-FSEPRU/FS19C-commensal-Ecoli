# Notes, Ideas and Thoughts

# To do
## FS19C
* GO term enrichment of roary and/or DRAM output?
  * Roary output: translate the `pan_genome_reference.fa` fasta file to AA sequence, and then run through interproscan.
* Discussions with Jules and Crystal: converting presence/absence data from `metabolism_summary.xlsx` into an ordination to see which commensals fall closer to STEC and could be candidates for further study. Can also color code and include virulence genes.
* ppangolin: get fasta file of gene families or protein families, find out which ones are present/absent in STEC, commensals - blast those fasta files


## Ceres
* You can change your default Slurm account using running slurm-account-selector.sh on the login node.
  * To see all your Slurm accounts at any time, use “sacctmgr -Pns show user format=account”
