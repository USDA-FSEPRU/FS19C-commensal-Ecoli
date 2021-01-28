#usr/bin/bash

mkdir linked
for myfile in ./*.fastq.gz; do
        target=$(echo $myfile|sed -r 's/1\-[A-H][0-9]+\-([0-9]+\-[0-9]+[A-Z]+[0-9]?[A-Z]?\_)[A-Z]+[0-9]+\_[A-Z]+[0-9]+\_R([0-9])\_[0-9]+/\1\2/g')
        echo "$myfile is being renamed to $target"
        ln -sr $myfile linked/$target
done
