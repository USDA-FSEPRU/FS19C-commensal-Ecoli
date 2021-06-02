# 08_DRAM
* Summary: Identify sugar utilization pathways in STECs that are also present in commensal E. coli to identify candidate commensal E. colis for further study.
* Platform: Ceres, conda
  * `/project/fsepru/kmou/FS19C/STECgenomes/`

1. Install DRAM via miniconda: https://github.com/shafferm/DRAM/wiki/2.-How-to-Install-and-Set-Up-DRAM in:
```
cd /project/fsepru/kmou/conda_envs_dram_analysis
wget https://raw.githubusercontent.com/shafferm/DRAM/master/environment.yaml
```

2. Edit `environment.yaml` to the following:
```
pkgs_dirs:
  - /project/fsepru/kmou/my_pkg_cache
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.*
  - pandas
  - pytest
  - scikit-bio
  - prodigal
  - mmseqs2!=10.6d92c
  - hmmer!=3.3.1
  - trnascan-se >=2
  - sqlalchemy
  - barrnap
  - altair >=4
  - openpyxl
  - networkx
  - ruby
  - parallel
  - dram
environment.yaml
```

4. Run `DRAM-setup.py print_config` in slurm script `dram3.slurm`.

<details><summary>dram3.slurm script</summary>

```
#!/bin/bash
#SBATCH --job-name=dram3                            # name of the job submitted
#SBATCH -p mem                                    # name of the queue you are submitting to
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=16
#SBATCH --mem=550gb
#SBATCH -t 96:00:00                                      # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
#SBATCH --account fsepru
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kathy.mou@usda.gov
#Enter commands here:
set -e
set -u
set +eu
module load miniconda
source activate /project/fsepru/kmou/conda_envs_dram_analysis/DRAM
DRAM-setup.py prepare_databases --output_dir DRAM_data3 --threads 16
```
</details>

5. Set up directory that has all STEC and commensal E. coli fasta files (`*.fna`). Make sure to exclude the following E. coli STEC reference genomes (used in previous steps) because they are already included in the RefSeq collection: `Ecoli_TW14588.fna`, `Ecoli_NADC6564.fna`, `Ecoli_O157H7_EDL933.fna`.

6. We have a total of 133 STEC genomes + 95 commensal E. coli + 3 commensal reference strains = 231 genomes. DRAM could only annotate 213 of the 231 genomes in one job.
  * DRAM used 20 processors for 80 MAGS and completed in ~17h. I have 231 genomes. If I use 32 cores, with 16 GB per core = 512 GB memory. Set memory at 550GB.

7. Run `DRAM.py annotate` on 231 E. coli genomes in 2 jobs.
<details><summary>dram3.slurm script</summary>

  ```
  #!/bin/bash
  #SBATCH --job-name=dram3                            # name of the job submitted
  #SBATCH -p mem                                    # name of the queue you are submitting to
  #SBATCH -N 1
  #SBATCH -n 32
  #SBATCH --mem=550gb
  #SBATCH -o "stdout.%j.%N.%x"                               # standard out %j adds job number to outputfile name and %N adds the node name
  #SBATCH -e "stderr.%j.%N.%x"                               # optional but it prints our standard error
  #SBATCH --account fsepru
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=kathy.mou@usda.gov

  #Enter commands here:
  set -e
  set -u
  set +eu

  module load miniconda
  source activate /project/fsepru/kmou/conda_envs_dram_analysis/DRAM
  DRAM.py annotate -i '/project/fsepru/kmou/FS19C/stecandcommensalEcoli_gifrop/*.fna' -o annotation_v3 --threads 32
  ```
</details>
