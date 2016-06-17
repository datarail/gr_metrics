# Case C (example 4) test
#install.packages('readr')      # un-comment and install these packages if necessary
#install.packages('devtools')   # un-comment and install these packages if necessary

# Load GRmetrics functions
devtools::load_all()
# Load Case C (example 4) input
data("inputCaseC")
# Run GRfit function with case = "C"
output4 = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "C")

# Load Case A (example 1) input
data("inputCaseA")
# Run GRfit function with case = "A"
output1 = GRfit(inputData = inputCaseA, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "A")

# change type integer to numeric for the sake of testing
output1[[1]]$replicate = as.numeric(output1[[1]]$replicate)
all.equal(output1, output4)
#[1] TRUE
# Test passed - output from Case C matches output from Case A

colnames(inputCaseC)
#[1] "cell_line"     "agent"         "perturbation"  "replicate"     "time"          "concentration" "cell_count"
colnames(inputCaseC)[6] = "conc"
colnames(inputCaseC)[7] = "count"
outputC = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "C", concentration = 'conc', cell_count = 'count')

test = readr::read_tsv('data/toy_example_input1.tsv')
colnames(test)
# [1] "cell_line"         "agent"             "perturbation"      "replicate"         "time"
# [6] "concentration"     "cell_count"        "cell_count__ctrl"  "cell_count__time0"

# changing column names just to test it out
colnames(test)[6] = "conc"
colnames(test)[7] = "count_end"
colnames(test)[8] = "count_ctrl"
colnames(test)[9] = "count_start"

test_output = GRfit(inputData = test, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), GRtable = 'both', case = "A", concentration = 'conc', cell_count = 'count_end', cell_count__time0 = 'count_start', cell_count__ctrl = 'count_ctrl')

test = readr::read_tsv('data/toy_example_input1.tsv')
colnames(test)[6] = "conc"
colnames(test)[7] = "count_end"
test_output = GRfit(inputData = test, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), GRtable = 'both', case = "A", concentration = 'conc', cell_count = 'count_end')


# change type integer to numeric for the sake of testing
test_output[[1]]$replicate = as.numeric(test_output[[1]]$replicate)
all.equal(test_output, output1)
#[1] TRUE
# Test passed - output from case A work when column names are changed

drc_output = GRfit(inputCaseA, groupingVariables = c('cell_line','agent'), case = "A")
GRdrawDRC(drc_output, experiments = c('BT20 drugA', 'MCF10A drugA', 'MCF7 drugA'), min = 10^(-4), max = 10^2)
GRdrawDRC(drc_output)
GRdrawDRC(drc_output, plotly = F)

GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB')
GRscatter(output1, 'GR50', 'agent', c('drugA','drugD'), 'drugB', plotly = F)

GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'))
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'), plotly = F)
GRbox(output1, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'), plotly = F)

