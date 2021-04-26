#to test a slurm script, allocate a debug node to make sure it can start running

salloc -N 1 -p debug -t 01:00:00
module load miniconda      #<= there's a way to check if this is the right module name
source activate prokka_env
ppanggolin workflow --anno Ecoligbkpath.txt
