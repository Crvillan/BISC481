## Question 5
##Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret")

## Initialization
library(DNAshapeR)
library(caret)
workingPath <- "/Users/ceciliavillanueva/Desktop/BISC481/CTCF/"

# Extract sample sequences
fn <- system.file("extdata", "CGRsample.fa", package = "DNAshapeR")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots- We generate ensemble plots for each of the parameters that we want to observe.
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)
