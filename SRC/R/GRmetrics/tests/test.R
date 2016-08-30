# Case C (example 4) test
#install.packages('readr')      # un-comment and install these packages if necessary
#install.packages('devtools')   # un-comment and install these packages if necessary
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("SummarizedExperiment")

# Load GRmetrics functions
library(GRmetrics)
# Load Case C (example 4) input
data("inputCaseC")
# Run GRfit function with case = "C"
output4 = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "C")

# Load Case A (example 1) input
data("inputCaseA")
# Run GRfit function with case = "A"
output1 = GRfit(inputData = inputCaseA, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "A")

# change type integer to numeric for the sake of testing
metadata(output1)[[1]]$replicate = as.numeric(metadata(output1)[[1]]$replicate)
all.equal(output1, output4)
#[1] TRUE
# Test passed - output from Case C matches output from Case A
