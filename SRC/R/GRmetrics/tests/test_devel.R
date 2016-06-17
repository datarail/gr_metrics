# development version test

#install.packages('readr')
library(readr)
test = read_tsv('data/toy_example_input1.tsv')
colnames(test)
# [1] "cell_line"         "agent"             "perturbation"      "replicate"         "time"
# [6] "concentration"     "cell_count"        "cell_count__ctrl"  "cell_count__time0"

# changing column names just to test it out
colnames(test)[6] = "conc"
colnames(test)[7] = "count_end"
colnames(test)[8] = "count_ctrl"
colnames(test)[9] = "count_start"

#install.packages('devtools')
library(devtools)
load_all()
gr_table = GRcalculate(test, 'count_start', 'count_end', 'count_ctrl')
GRfitData = GRlogisticFit(gr_table, c('cell_line','agent', 'perturbation','replicate', 'time'), 'conc')
#GRfitData = GRfit(test, c('cell_line','agent', 'perturbation','replicate', 'time'))
GRscatter(GRfitData, 'GR50', 'agent', c('drugA','drugD'), 'drugB', plotly = F)
# considering deleting GRscatterAdd
#GRscatterAdd(GRfitData, 'GR50', 'agent', 'drugB', 'drugB', plotly = F)
GRbox(GRfitData, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = "all", plotly = T)

GRbox(GRfitData, GRmetric ='GRinf', groupVariable = 'cell_line', pointColor = 'agent' , factors = c('BT20', 'MCF10A'), plotly = F)

########
gr_table = GRcalculate(test, 'count_start', 'count_end', 'count_ctrl')
GRfitData = GRlogisticFit(gr_table, c('cell_line','agent'), 'conc')
GRdrawDRC(gr_table, GRfitData, c('cell_line','agent'), 'conc', min = "auto", max = "auto", points = T, curves = T, plotly = T)
