# GRmetrics

An R package for calculating Growth Rate Inhibition Metrics 

## Installation

Install the latest development version (on GitHub) via devtools:

```r
install.packages("devtools")
devtools::install_github("uc-bd2k/GRmetrics")
```

####References:
Hafner, M., Niepel, M. Chung, M. and Sorger, P.K., *Growth Rate Inhibition Metrics Correct For Confounders In Measuring Sensitivity To Cancer Drugs*. Nature Methods 13.6 (2016): 521-527.
(http://dx.doi.org/10.1038/nmeth.3853)

Corresponding MATLAB and Python scripts available on repo https://github.com/sorgerlab/gr50_tools 

Browser interface and online tools: http://www.grcalculator.org

Much of the description below has been adapted from https://github.com/sorgerlab/gr50_tools/blob/master/README.md. See there for additional details on corresponding MATLAB and Python scripts.

###Input data
The main function of the package is **GRfit**, which takes in a data frame containing information about concentration, cell counts over time, and additional grouping variables for a dose-response assay and calculates growth-rate inhibition (GR) metrics for each experiment in the dataset.

There are two cases of data input accepted by the **GRfit** function. They are described in detail below. Case "A" is the default option.

**The mandatory inputs for Case "A" are:**

- **inputData**   the name of the input data frame
- **concentration**		the name of the column with concentration values (not log transformed) of the perturbagen on which dose-response curves will be evaluated
- **cell_count**		the name of the column with the measure of cell number (or a surrogate of cell number) after treatment
- **cell_count\_\_time0**   the name of the column with Time 0 cell counts - the measure of cell number in untreated wells grown in parallel until the time of treatment
- **cell_count\_\_ctrl**    the name of the column with the Control cell count - the measure of cell number in control (e.g. untreated or DMSO-treated) wells from the same plate

All other columns will be treated as additional keys on which the data will be grouped (e.g. *cell_line*, *drug*, *time*, *replicate*)

**The default values for the column variables (concentration, cell_count, cell_count\_\_time0, and cell_count\_\_ctrl) are equal to their variable names, so you do not have to specify these arguments unless your columns are named differently.**

See the following link for an example input data file:

- https://github.com/sorgerlab/gr50_tools/blob/master/INPUT/toy_example_input1.tsv

This data frame is stored in the package as "inputCaseA":
```r
library(GRmetrics)
data(inputCaseA)
```

**The mandatory inputs for Case "C" are:**

- **inputData**   the name of the input data frame
- **concentration**		the name of the column with concentration values (not log transformed) of the perturbagen on which dose-response curves will be evaluated
- **cell_count**		the name of the column with the measure of cell number (or a surrogate of cell number)
- **time** 			the name of the column with the time at which a cell count is observed

All other columns will be treated as additional keys on which the data will be grouped (e.g. *cell_line*, *drug*, *time*, *replicate*)

**The default values for the column variables (concentration, cell_count, time) are equal to their variable names, so you do not have to specify these arguments unless your columns are named differently.**

See the following link for an example input data file:

- https://github.com/sorgerlab/gr50_tools/blob/master/INPUT/toy_example_input4.tsv

This data frame is stored in the package as "inputCaseC":
```r
library(GRmetrics)
data(inputCaseC)
```

###Visualization functions
The package contains 3 visualization functions: **GRdrawDRC**, **GRscatter**, and **GRbox**.

All of these functions take in an object created by **GRfit** as well as additional arguments. The results can be viewed in a static ggplot image or an interactive plotly (turned on/off by the *plotly* parameter).

- **GRdrawDRC**   this function draws the (growth-rate inhibition) dose-response curve using the parameters calculated by the **GRfit** function. If *points* is set to TRUE, then it will also plot the points used to fit each curve.
- **GRscatter**   this function draws a scatterplot of a given GR metric (GR50, GRmax, etc.) with the *xaxis* value(s) plotted against the *yaxis* value(s).
- **GRbox**   this function draws boxplots of a given GR metric (GR50, GRmax, etc.) for values of a given grouping variable. It overlays the points used to make these boxplots and can color them according to another grouping variable.

### Examples
```r
## Case A (DRC examples)
library(GRmetrics)
data(inputCaseA) #load the data
head(inputCaseA) #review the data
# Calculate GR values and solve for GR metrics parameters (i.e. fit curves)
drc_output = GRfit(inputCaseA, groupingVariables = c('cell_line','agent'))
head(drc_output$gr_table) #review output table of GR values
head(drc_output$parameter_table) #review output table of GR metrics parameters

# Draw dose-response curves
GRdrawDRC(drc_output)
GRdrawDRC(drc_output, experiments = c('BT20 drugA', 'MCF10A drugA', 'MCF7 drugA'))
GRdrawDRC(drc_output, experiments = c('BT20 drugA', 'MCF10A drugA', 'MCF7 drugA'), min = 10^(-4), max = 10^2)
GRdrawDRC(drc_output, plotly = F)

## Case C (scatterplot and boxplot examples)
data(inputCaseC)
head(inputCaseC)
output1 = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "C")

# Draw scatterplots
GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB')
GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB', plotly = F)

# Draw boxplots
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent')
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'))
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'), plotly = F)

# Input with different column names
data(inputCaseC)
colnames(inputCaseC)
#[1] "cell_line"     "agent"         "perturbation"  "replicate"     "time"          "concentration" "cell_count"
colnames(inputCaseC)[6] = "conc"
colnames(inputCaseC)[7] = "count"
outputC = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent'), case = "C", concentration = 'conc', cell_count = 'count')
GRdrawDRC(outputC)

```

### Input data formats

####Case A: a single file with control values assigned to treated measurements
The control values (both control and time 0 cell counts) are pre-computed by the user and assigned to each treatment (row) in appropriate columns in the input file. Control cell counts should be in a column labeled *cell_count\_\_ctrl* and the time 0 cell counts in a column labeled *cell_count\_\_time0*.

This case corresponds to the toy example 1 in the GitHub folder.
https://github.com/sorgerlab/gr50_tools/blob/master/INPUT/toy_example_input1.tsv

This data frame is stored in the package as "inputCaseA":
```r
library(GRmetrics)
data(inputCaseA)
```
####Case B: three files with control values labelled with a key
Implemented in MATLAB code. Not yet implemented in R package or online calculator.

####Case C: a single file with control values stacked with treated measurements
In the most general case, the control cell counts are in the same file and format as the treated cell counts. Control cell counts will be averaged (using a 50%-trimmed mean) and automatically matched to the treated cell counts based on the keys (columns in the data file). The control cell count values must have a value of 0 for *concentration* and a value for *time* that matches the treated measurements. The time 0 cell count values must have value of 0 for *time*. If the structure of the data is complex, the provided scripts may inappropriately match control and treated cell counts, so users instead should format their data as described in case A or B. 

Case C corresponds to the toy example 4 in the GitHub folder.
https://github.com/sorgerlab/gr50_tools/blob/master/INPUT/toy_example_input4.tsv

This data frame is stored in the package as "inputCaseC":
```r
library(GRmetrics)
data(inputCaseC)
```

##General approach
We have developed scripts to calculate normalized growth rate inhibition (GR) values and corresponding metrics (GR_50, GR_max, ...) based on cell counts measured in dose-response experiments. Users provide a tab-separated data file in which each row represents a separate treatment condition and the columns specify the keys that define the treatment condition (e.g. cell line, drug or other perturbagen, perturbagen concentration, treatment time, replicate) and the measured cell counts (or surrogate). The experimentally measured cell counts that are required for GR metric calculation are as follows: 
- measured cell counts after perturbagen treatment (*cell_count*, *x(c)*)
- measured cell counts of control (e.g. untreated or DMSO-treated) wells on the same plate (*control_cell\_\_count*, *x_ctrl*)
- measured cell counts from an untreated sample grown in parallel until the time of treatment (*time0_cell\_\_count*, *x_0*)

The provided GR scripts compute over the user’s data to calculate GR values individually for each treatment condition (cell line, time, drug, concentration, ...) using the formula:

    GR(c) = 2 ^ ( log2(x(c)/x_0) / log2(x_ctrl/x_0) ) - 1

Based on a set of GR values across a range of concentrations, the data are fitted with a sigmoidal curve:

    GR(c) = GR_inf + (1-GR_inf)/(1 + (c/(GEC_50))^h_GR )

The following GR metrics are calculated: 
- **GR_inf** = GR(c->inf), which reflects asymptotic drug efficacy. 
- Hill coefficient of the sigmoidal curve (**h_GR**), which reflects how steep the dose-response curve is.
- **GEC_50**, the drug concentration at half-maximal effect, which reflects the potency of the drug.
- **GR_50**, the concentration at which the effect reaches a GR value of 0.5 based on interpolation of the fitted curve.
- **GR_AOC**, the area over the dose-response curve, which is the integral of *1-GR(c)* over the range of concentrations tested, normalized by the range of concentration. 
- **GR_max**, the effect at the highest tested concentration. Note that *GR_max* can differ from *GR_inf* if the dose-response does not reach its plateau value.

In addition, the scripts report the r-squared of the fit and evaluate the significance of the sigmoidal fit based on an F-test. If the fit is not significant (p<0.05, or any arbitrary value), the sigmoidal fit is replaced by a constant value (flat fit). The cutoff value for the p-value can be set above 1 for bypassing the F-test. Additional information and considerations are described in the supplemental material of the manuscript referred above.

