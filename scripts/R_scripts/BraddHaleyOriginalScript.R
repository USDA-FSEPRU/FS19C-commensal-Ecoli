
# Set working directory where data resides.

setwd("g:/client projects/haley/20009")

# Read the row of MDR and SUSC labels from the CSV spreadsheet, to identify the range of columns for each label.

lbltype <- read.csv("SampleMatrix.csv",header=FALSE,skip=1,nrows=1)
mdr <- lbltype[]=="MDR"
mdr
dim(mdr)
mdr[1:101]  # Columns   2:100 (i.e.,  99) are MDR
mdr[99:217] # Columns 101:217 (i.e., 117) are SUSC

# Skip column headings and MDR/SUSC labels; read 1/0 presence/absence data values.

genomes <- read.csv("SampleMatrix.csv",header=FALSE,skip=2)

# Assign genome names from the first column to be the rownames for the data object.

rownames(genomes) <- genomes[,1]

# Create an empty list.  Each element of this list will store Fisher Exact Test results one genome.

lst.exact2x2 <- c()


# VERY IMPORTANT:  Make certain that when you begin your R or Rstudio session, you right-click on the icon
#                  and select "Run as administrator".  Otherwise, package installations will be unsuccessful.

# Install and load the exact2x2 package.
# exact2x2 produces p-values, odds ratios, and odds ratio confidence intervals.

# Note:  Need install the package only once; not everytime you begin a new R session.
# install.packages("exact2x2")
library(exact2x2)

# Loop over each genome to sum the present and absent counts for MDR and for SUSC

for (i in 1:nrow(genomes))
{
  # Count of MDR Samples
  genomes[i,217] <- sum(genomes[i,2:100])

  # Count of Not MDR Samples
  genomes[i,218] <- 99-genomes[i,217]

  # Count of SUSC Samples
  genomes[i,219] <- sum(genomes[i,101:216])

  # Count of Not SUSC Samples
  genomes[i,220] <- 116-genomes[i,219]

  # Create a 2x2 matrix containing presence/absense counts to be used by exact2x2
  xi <- matrix(c(genomes[i,217],genomes[i,218],genomes[i,219],genomes[i,220]),nrow=2,ncol=2)

  # Save exact test results for each genome into an element of the list.
  lst.exact2x2[[i]] <- exact2x2(xi,tsmethod="central")

  # Assign the Genome ID to an element of the list, for ID and labeling purposes.
  lst.exact2x2[[i]]$genomeID <- as.character(genomes[i,1])
}

# Count Sums and proportions

count.sums <- data.frame(genomes[,217],genomes[,218],genomes[,219],genomes[,220])
colnames(count.sums) <- c("MDR.yes","MDR.no","SUSC.yes","SUSC.no")

prop.mdr <- genomes[,217]/(genomes[,217]+genomes[,218])
prop.susc <- genomes[,219]/(genomes[,219]+genomes[,220])

count.sums <- data.frame(count.sums, prop.mdr, prop.susc)


# Create a data frame from the content of the exact2x2 list.

p.values <- sapply(lst.exact2x2, function(x){as.numeric(x[1])})

orci <- sapply(lst.exact2x2, function(x){as.numeric(x[2][[1]])})
str(orci)
orlcl <- orci[1,]
orucl <- orci[2,]

or <- sapply(lst.exact2x2, function(x){as.numeric(x[3])})

gid <- sapply(lst.exact2x2, function(x){as.character(x[8])})


# "exact.results" contains all results from the exact2x2 Fisher's Exact tests, for each genome.
# BUT, these p-values should be adjusted to protect against false significance caused by conducting
# X individual tests, each with 5% probablity of false significance.

exact.results <- data.frame(gid, p.values, count.sums, or, orlcl, orucl)

original.order <- seq(1:nrow(exact.results))


# The following code conducts an FDR adjustment (i.e., calculates q-values) for each genome.

# Install these 2 packages, required for subsequently installing and using the qvalue package.

# Note:  Need install the package only once; not everytime you begin a new R session.
# install.packages("devtools","BiocManager")

# Load the packages.
library("devtools","BiocManager")

# Install the qvalue package, to adjust the exact2x2 p-values for false discovery rate (FDR).
# Note:  Need install the package only once; not everytime you begin a new R session.
BiocManager::install("qvalue")

# Load the qvalue package.
library(qvalue)

# The vignette can be viewed by typing:
# browseVignettes(package = "qvalue")

q.obj <- qvalue(p=exact.results$p.values)

# The summary function indicates the number of genomes exhibiting significance based on each
# of the 3 criteria:  p-value, q-value, local FDR;  where column heading indicates alpha level.

summary(q.obj)

#plot(q.obj)

q.obj$pi0

qv <- q.obj$qvalues

q.exact.results <- data.frame(original.order, exact.results, qv)

sigq.order <- order(q.exact.results$qv)

ord.q.exact.results <- q.exact.results[sigq.order,]

ord.q.exact.results$sig.order <- seq(1:nrow(ord.q.exact.results))

# Write the results to a csv file.
write.csv(ord.q.exact.results, "Genome 2x2 Exact Test q-values - Sorted from Most to Least Significant.csv")
