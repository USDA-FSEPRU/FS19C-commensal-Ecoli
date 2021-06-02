# 07_ppanggolin
* Summary: Run pan-genome analysis with ppanggolin and prokka-annotated gbk files. Identify gene families involving sugar utilization pathways that are present in both STEC and commensal E. coli. This generates a nice summary plot of genes families present and absent in pan-isolates.
* Platform: Ceres, conda
  * `/project/fsepru/kmou/FS19C/stecandcommensalEcoli_gifrop/gbkforppanggolin_prokka`

1. Make directory `gbkforppanggolin_prokka` and copy gbk files from prokka annotation to this directory.

2. Remove Ecoli_TW14588.gbk, Ecoli_NADC6564.gbk, Ecoli_O157H7_EDL933.gbk from directory.

3. Create text file listing all filenames and their paths in `gbkforppanggolin_prokka/` directory and downloaded text file `Ecoligbkpath.txt`.
```
ls -d -1 "$PWD/"*.gbk > Ecoligbkpath.txt
```

4. Open `Ecoligbkpath.txt` in Excel. Goal is to have 1st column containing unique organism name and second column as path to location of gbk file. In the first column, delete the path name `/project/fsepru/kmou/FS19C/stecandcommensalEcoli_gifrop/gbkforppanggolin_prokka/`, `_pol.gbk`, `.gbk`. Save and upload modified text file to Ceres.

5. Create `ppanggolin.slurm` script on Ceres and run script on slurm.

  <details><summary>ppanggolin.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=ppanggolin                             # name of the job submitted
  #SBATCH -p short                                    # name of the queue you are submitting to
  #SBATCH -N 1                                            # number of nodes in this job
  #SBATCH -n 16                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
  #SBATCH -t 48:00:00                                      # time allocated for this job hours:mins:seconds
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --mem=32G   # memory
  #SBATCH --account fsepru
  #SBATCH --mail-user=kathy.mou@usda.gov
  #SBATCH --mail-type=ALL
  #Enter commands here:
  set -e
  set -u
  set +eu
  module load miniconda
  source activate /project/fsepru/conda_envs/ppanggolin
  ppanggolin workflow --anno Ecoligbk.txt
  ```
  </details>

6. Download `ppanggolin_output_DATE*` and look at output, such as `tile_plot.html`. This plot has a nice summary heatmap of gene families present or absent in each isolate of pan-genome.

7. On Ceres, paste the following command in `ppanggolin.slurm`. Move slurm script to `ppanggolin_output_DATE*` directory and submit slurm job.
```
ppanggolin write -p pangenome.h5 --families_tsv --output gene_families
```

8. `gene_families.tsv` shows all the genes (non-descriptive) for each gene family. Will have to blast fasta files. Can obtain fasta file for entire pangenome of genes, gene families, or protein families:
```
ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes all
ppanggolin fasta -p pangenome.h5 --output MY_GENES --gene_families all
ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families all
```
