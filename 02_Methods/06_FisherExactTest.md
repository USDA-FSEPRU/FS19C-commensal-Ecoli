# 06 Fisher Exact Test
* Summary: Use output from Roary to run through R script.
* Platform: R

1. Modify `gene_presence_absence.Rtab` similar to `https://github.com/k39ajdM2/Notebook/tree/main/Files/BraddHaley_SampleMatrix.csv` and save as `Matrix.csv`. Modifications include:
  * Rename 1st column heading to `ID`.
  * Add a second row `Rstatus` and label respective group name for each strain. Try to group strains of the same group together (helps make it easier to call specific columns for a group in R script)

2. Run `Matrix.csv` through `FisherExactTest.R`. You will need to modify the code to fit with number of columns and rows for each group and gene in `Matrix.csv`

  <details><summary>FisherExactTest.R script</summary>



  </details>
