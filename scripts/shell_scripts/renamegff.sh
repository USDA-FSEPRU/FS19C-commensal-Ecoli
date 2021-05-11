#!/bin/bash

# This script renames the `genes.annotated.gff` from each E. coli annotated genome with directory name.gff
# genes.annotated.gff -> 62-438RED_pol.gff

for pathname in */genes.annotated.gff; do
  cp "$pathname" "$( basename "$( dirname "$pathname" )" ).gff"
done
