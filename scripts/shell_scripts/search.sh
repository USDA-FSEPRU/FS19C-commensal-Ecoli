#!/bin/bash

## This code helps to find E. coli virulence genes that are listed in Ecolivirulencegene.txt and search through gene_presence_absence_2.txt and generate pa.isolates.txt output showing search results

echo "Beginning of file" > pa.isolates.txt
for line in `cat Ecolivirulencegene.txt`;
do
	echo $line >> pa.isolates.txt
	grep -n "$line" gene_presence_absence_2.txt >> pa.isolates.txt
	echo "END_OF_GENE" >> pa.isolates.txt

done
echo "End of file" >> pa.isolates.txt


# to be able to branch and whatnot, would have to do something like this
# for x in file.txt:
# 	gene = x;
# 	output = grep gene in O157h7.txt
# 	n = count number of ":"
# 	line numbers = regex{\d+:}
# 	if n==0:
#		print(gene, "not found")
# 	else:
#		print("found",n,"times on lines:", output)
