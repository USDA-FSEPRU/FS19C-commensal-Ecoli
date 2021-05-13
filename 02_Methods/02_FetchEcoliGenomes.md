# 02_Fetch E. coli genomes
* Summary: Fetch STEC and non-STEC *Escherichia coli* genomes from RefSeq.
* Platform: Ceres
  * `/project/fsepru/kmou/FS19C/STECgenomes`

1. Copy commensal E. coli polished fasta files to `/project/fsepru/kmou/FS19C/STECgenomes`.

2. Download fasta files of the following *E. coli* reference genomes, add with your polished commensal *E. coli* fasta files.
  * E. coli MG1655 https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2
    * Saved as Ecoli_K12_MG1655_GCF_000005845.2_ASM584v2_genomic.gff.gz
  * E. coli HS https://www.ncbi.nlm.nih.gov/assembly/GCF_000017765.1
    * Saved as Ecoli_HS_GCF_000017765.1_ASM1776v1_genomic.gff.gz
  * E. coli Nissle 1917 https://www.ncbi.nlm.nih.gov/assembly/GCF_000714595.1
    * Saved as Ecoli_Nissle1917_GCF_000714595.1_ASM71459v1_genomic.gff.gz
  * E. coli O157:H7 str. NADC 6564 https://www.ncbi.nlm.nih.gov/assembly/GCF_001806285.1
    * Saved as Ecoli_O157H7_NADC_6564_GCF_001806285.1_ASM180628v1_genomic.gff.gz
  * (This is already available on RefSeq so will ignore this version) E. coli O157:H7 EDL933 https://www.ncbi.nlm.nih.gov/assembly/GCF_000732965.1
    * Saved as EcoliO157H7_EDL933_GCF_000732965.1_ASM73296v1_genomic.gff.gz
  * TW14588 https://www.ncbi.nlm.nih.gov/assembly/GCF_000155125.1/
    * Saved as Ecoli_TW14588_GCF_000155125.1_ASM15512v1_genomic.gff.gz

3. Look for STEC genomes on NCBI (RefSeq).
  * Read https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/README.txt and followed directions on how to grab assembled sequence data of STEC strains from NCBI.
  * Looked for latest assembly versions: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/README.txt
  ```
  all_assembly_versions & latest_assembly_versions are not listed for any species with more than 1,000 assemblies.
  Use the assembly_summary.txt file in the species directory to find assemblies of interest
  and then access the data using the ftp_path from column #20 of the file.
  ```

4. Download `assembly_summary.txt` via `wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/assembly_summary.txt`. Save as `assembly_summary.xlsx`.

5. Make the following sheets in `assembly_summary.xlsx`:
  * `assembly_summary`: all assembly data.
  * `Complete Genome`: filtered from `assembly_summary` sheet, only includes complete genomes based on column L `assembly_level`.
  * `O157H7`: filtered from `Complete Genome` sheet, only includes  `Escherichia coli O157:H7` based on column H `organism_name`.
    * total of 133 STEC genomes

6. Within `O157H7` sheet in `assembly_summary.xlsx`, create columns `desired_file_from_directory` and `new_ftp_path`.

7. Under column `desired_file_from_directory`, fill each cell with `_genomic.fna.gz`.

8. Under `new_ftp_path`, combine the following columns using excel combine function: `# assembly_accession` (column A), `asm_name` (column P), `ftp_path` (column T), and `desired_file_from_directory` (column W).
  ```
  =T3&"/"&A3&"_"&P3&W3
  ```

9. Copy new paths under `new_ftp_path` (column W) to text file and save as `stecftp.txt`. Upload file to Ceres.

10. Run `stec.slurm` to copy `*genomic.fna.gz` from RefSeq to `/project/fsepru/kmou/FS19C/STECgenomes`.

  <details><summary>stec.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=stec                             # name of the job submitted
  #SBATCH -p short                                 # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 2                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  #Enter commands here:
  wget -i stecftp.txt
  ```
  </details>

11. Unzip fetched files via `gzip -d *genomic.fna.gz`.
