# Fancy bash scripts
* Source of bash scripts: https://github.com/Jtrachsel/bash_fun

##rename_contigs
### How to execute command:
`rename_contigs 27-433FED_pol.fna 27-433FED`

### rename_contig script:
```
rename_contigs() {
        FILE=$1
        BASE=$2
        awk -v awkv="$BASE" '/^>/{print ">"awkv"_"++i;next}{print}' "$FILE" > "$BASE"_rename.fasta
        rm $FILE
        mv "$BASE"_rename.fasta $FILE
}
export -f rename_contigs
```

### Breakdown of rename_contigs
```
rename_contigs() {
        FILE=$1
        BASE=$2

###

awk -v basev="$BASE"      #<=declare variable basev as $2
 '
  /^>/      #<=find a line that starts with >
  {
    print ">"   basev   "_"  ++i      #<=replace that line with >basev_i, where i is an index that increments each match. i can be called whatever, except certain restricted words (like index)
    ; #separates two commands
    next      #<=go to the next line, executes next iteration of loop
  }
 {print}      #<=reproduces line it finds (only gets executed if line doesn't start with >)
 '
 "$FILE"      #<=tells awk to operate on this file
 > "$BASE"_rename.fasta      #<=rename $FILE to $BASE_rename.fasta

###

 rm $FILE #<=remove $FILE
        mv "$BASE"_rename.fasta $FILE      #<=move "$BASE"_rename.fasta to where $FILE was
}
export -f rename_contigs      #<=execute this function in ~/.bashrc so you it will recognize and run rename_contigs in the future
```

### for-loop to run rename_contigs on multiple fasta files in a directory
```
for file in ./*.fna
  do
   file2=$(basename ${file} | sed 's/\.fna$//g' | sed -r 's/\_pol//g')
   rename_contigs ${file} ${file2}
done
```

### Breakdown of for-loop
```
for file in ./*.fna      #<= read one *.fna as 'file'
  do
   file2=$(basename ${file} | sed 's/\.fna$//g' | sed -r 's/\_pol//g')      #<= basename removes "./", sed replaces .fna with nothing, then replaces _pol with nothing
   rename_contigs ${file} ${file2}      #<=run rename_contigs with file and file2
done
```

Example filename to model for-loop off of: 66-440RED_pol.fna

### In pseudocode that would look like:
for each line in file {
  if (line starts with >) {print ">basev_index"; goto next item in for loop}
  print line #only gets executed if "if" statement fails
}
