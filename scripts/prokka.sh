#usr/bin/bash

for file in *.fasta; do tag=$file%.fasta; prokka -prefix "$tag" -locustag "$tag" -genus Escherichia -strain "$tag" -outdir "$tag"_prokka -force -addgenes "$file" -centre X -compliant; done
