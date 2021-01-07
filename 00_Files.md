# Files generated for FS19C project and their location


## Raw data location

### Illumina Reads
Q drive

### QC reports

#### MultiQC
* Completed on Jan. 7, 2020
* Platform: Ceres, conda environment
```
salloc
module load miniconda
source activate fastanienv
conda install -c bioconda multiqc
multiqc *.fastqc.zip
```
* File generated:  multiqc_report.html

* Why are the plots flat plot?
From: https://multiqc.info/docs/
"Flat plots
Reports with large numbers of samples may contain flat plots. These are rendered when the MultiQC report is generated using MatPlotLib and are non-interactive (flat) images within the report. The reason for generating these is that large sample numbers can make MultiQC reports very data-intensive and unresponsive (crashing people's browsers in extreme cases). Plotting data in flat images is scalable to any number of samples, however.
Flat plots in MultiQC have been designed to look as similar to their interactive versions as possible. They are also copied to multiqc_data/multiqc_plots"


## Genome Assembly
Ceres


## Genome Annotation



## Final Files of Importance
