#!/bin/bash

## This code adds _pol to *.fna fasta file name. Also run rename_contigs.sh as a for-loop

#66-440RED_pol.fna

for file in ./*.fna
  do
   file2=$(basename ${file} | sed 's/\.fna$//g' | sed -r 's/\_pol//g')
   rename_contigs ${file} ${file2}
done
