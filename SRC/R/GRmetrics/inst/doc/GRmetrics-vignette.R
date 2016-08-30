## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ---- include = FALSE----------------------------------------------------
## Case A (DRC examples)
library(GRmetrics)

## ------------------------------------------------------------------------
data(inputCaseA)

## ---- include = FALSE----------------------------------------------------
inputCaseA = as.data.frame(inputCaseA)

## ------------------------------------------------------------------------
head(inputCaseA)

## ---- include = FALSE----------------------------------------------------
drc_output = GRfit(inputCaseA, groupingVariables = c('cell_line','agent'))

## ------------------------------------------------------------------------
drc_output

## ------------------------------------------------------------------------
assay(drc_output)

## ------------------------------------------------------------------------
colData(drc_output)

## ------------------------------------------------------------------------
head(metadata(drc_output)[[1]])

## ------------------------------------------------------------------------
metadata(drc_output)[[2]]


## ------------------------------------------------------------------------
# Draw dose-response curves
GRdrawDRC(drc_output)
GRdrawDRC(drc_output, experiments = c('BT20 drugA', 'MCF10A drugA', 
                                      'MCF7 drugA'))
GRdrawDRC(drc_output, experiments = c('BT20 drugA', 'MCF10A drugA', 
                                      'MCF7 drugA'), 
          min = 10^(-4), max = 10^2)
GRdrawDRC(drc_output, plotly = FALSE)

## ------------------------------------------------------------------------
## Case C (scatterplot and boxplot examples)
data(inputCaseC)

## ---- include = FALSE----------------------------------------------------
inputCaseC = as.data.frame(inputCaseC)

## ------------------------------------------------------------------------
head(inputCaseC)

## ---- include = FALSE----------------------------------------------------
output1 = GRfit(inputData = inputCaseC, groupingVariables = 
                  c('cell_line','agent', 'perturbation', 'replicate', 'time'), 
                case = "C")

## ------------------------------------------------------------------------
# Draw scatterplots
GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB')
GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB', 
          plotly = FALSE)

# Draw boxplots
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', 
      pointColor = 'agent')
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', 
      pointColor = 'agent',
      factors = c('BT20', 'MCF10A'))
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', 
      pointColor = 'agent',
      factors = c('BT20', 'MCF10A'), plotly = FALSE)

