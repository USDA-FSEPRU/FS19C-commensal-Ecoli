#!/bin/bash

## This code, from Jules, renames the contig ID (first line) of fasta file

rename_contigs() {
        FILE=$1
        BASE=$2
        awk -v basev="$BASE" '/^>/{print ">"basev"_"++i;next}{print}' "$FILE" > "$BASE"_rename.fasta
        rm $FILE
        mv "$BASE"_rename.fasta $FILE
}
export -f rename_contigs
