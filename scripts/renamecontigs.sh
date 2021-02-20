#!/bin/bash

#66-440RED_pol.fna

for file in ./*.fna
  do
   file2=$(basename ${file} | sed 's/\.fna$//g' | sed -r 's/\_pol//g')
   rename_contigs ${file} ${file2}
done
