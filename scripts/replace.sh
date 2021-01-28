#usr/bin/bash


while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed "s/REPLACE/$line/g" > "$line".slurm
done < samples.txt
