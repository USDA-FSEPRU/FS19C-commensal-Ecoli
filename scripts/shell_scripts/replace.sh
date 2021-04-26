#usr/bin/bash

## This code replaces all "REPLACE" names with list of sample names in samples.txt and generate separate slurm file for each sample.

while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
