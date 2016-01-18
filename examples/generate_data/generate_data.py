import math as ma
import random as rn
import pandas as pd
import numpy as np

cap_error = 1.1;
std_error = 0.05;

# sigmoidal function for generating the data (with some noise)
def cellcount(conc, max_count, EC50, Hill, min_count):
        cc = min_count + (max_count - min_count)/(1. + ((10.**conc)/(10.**EC50))**Hill)        
        # capped random distribution 
        cc_rand = rn.normalvariate(cc, cc*std_error)
        return max(cc/cap_error, min(cc*cap_error, cc_rand)) # (capped to 10% error)

def trimmean(vals):
        return np.mean(sorted(vals)[int(len(vals)*.25):int(len(vals)*.75)])
        

# step for concentrations
log10_concentrations = [i/2. for i in range(-6,3)]
# durg names (generic)
agents = [u'drugA', u'drugB', u'drugC', u'drugD']

# initialization
ctrl_tag = 0
time0_tag = 0
output = [];
Drug_metrics = [];

# inputs are: cell line name, division rate in different conditions, cell_count__time0
for iC,r in enumerate((('MCF10A', (1, .7), 300),
        ('MCF7', (.6, .55), 400),
        ('BT20', (.6, .4), 500))   ):
    
    # initialize the seed for consistency
    rn.seed(4*iC-1)
    
    # define the ideal properties of the drug
    agents_log10_EC50 = [rn.uniform(-1.8, .5)  for i in range(len(agents))]
    agents_GRinf = [rn.uniform(-.9, .8)  for i in range(len(agents))]
    agents_Hill = [rn.uniform(1, 3)  for i in range(len(agents))]

    if iC==2:
        # this is an outlier case for triggering the F-test
        agents_GRinf[1] = .98;
    
    for ia,agent in enumerate(agents):
        # save the drug results in the table
        Drug_metrics.append((r[0], agent, 10**agents_log10_EC50[ia], agents_GRinf[ia], agents_Hill[ia]))
        
    for d in range(len(r[1])):
        # time0 cell count: define 20 time_0_count with some noise                        
        time0s = [ int(max(r[2]/cap_error, min(r[2]*cap_error,rn.normalvariate(r[2], r[2]*std_error)))) for i in range(20) ]
        # takes the trimmed mean
        time_0_count = trimmean(time0s)
        time0_tag += 1

        # save results in the table
        for time0 in time0s:
            output.append((r[0], u'-', d, u'NaN', 0, 0, 
                time0, u'NaN', time_0_count, 'c_%i' % ctrl_tag, 't0_%i' % time0_tag, u'NaN'))

                            
        for time in (2.,3.):
            # two time points
            for rep in range(1,3+(iC>1)):
                    # maximum number of cells
                    maxcount = int(time_0_count*2.**(r[1][d]*time))
                    
                    # define 8 control_count with some noise                        
                    ctrls = [ int(max(maxcount/cap_error, min(maxcount*cap_error, rn.normalvariate(maxcount, maxcount*std_error)))) for i in range(8) ]
                    
                    # takes the trimmed mean
                    control_count = trimmean(ctrls)
                    ctrl_tag += 1
                    
                    # save results in the table
                    for ctrl in ctrls:
                            output.append((r[0], u'-', d, rep, time*24., 0, 
                                ctrl, control_count, time_0_count, 'c_%i' % ctrl_tag, 't0_%i' % time0_tag, u'NaN'))
                    
                    # rg: relative growth
                    rg_ctrl = 1.*control_count/time_0_count
                    log2rg_ctrl = ma.log(rg_ctrl, 2)

                    for ia,agent in enumerate(agents):
                            # properties of the drug response with some noise
                            mincount = round(time_0_count * ((control_count/time_0_count)**(ma.log(agents_GRinf[ia]+1,2))))
                            
                            for c in log10_concentrations:
                                    # results of the sigmoidal function (with some noise)
                                    count = int( cellcount(c, maxcount, agents_log10_EC50[ia], agents_Hill[ia], mincount) )
                                    rg = 1.*count/time_0_count
                                    log2rg = ma.log(rg, 2)

                                    # GRval: normalized growth-rate inhibition
                                    GRval = '%.6g' % (2.**(log2rg/log2rg_ctrl) - 1)
                                            
                                    # save results in the table
                                    output.append((r[0], agent, d, rep, time*24, '%.4g' % (10.**c), 
                                        count, control_count, time_0_count, 'c_%i' % ctrl_tag, 't0_%i' % time0_tag, GRval))


# compile the data in a table
headers = ('cell_line agent perturbation replicate time concentration cell_count cell_count__ctrl cell_count__time0 ctrl_tag time0_tag GRvalue'.split())
df_data = pd.DataFrame(output, columns=headers)

headers = ('cell_line agent EC50 GRinf Hill'.split())
df_drugs = pd.DataFrame(Drug_metrics, columns=headers)

## printing out all the files with data to evaluate GR values

# print out the case 1 example (already assigned controls)
df_1 = df_data.loc[df_data.agent != '-', :'cell_count__time0']
df_1.to_csv('../../INPUT/toy_example_input1.tsv', '\t', float_format='%.4g', index=False)

# print out the case 2 example (collapsed controls in separate files)
#       and the case 3 example (all controls in separate files)
df_trt = pd.concat((df_data.loc[df_data.agent != '-', :'cell_count'], df_data.loc[df_data.agent != '-', ('ctrl_tag', 'time0_tag') ]), axis=1)
df_trt.to_csv('../../INPUT/toy_example_input2_data.tsv', '\t', float_format='%.4g', index=False)
df_trt.to_csv('../../INPUT/toy_example_input3_data.tsv', '\t', float_format='%.4g', index=False)

df_t0 = df_data.loc[df_data.time == 0, ('cell_line', 'perturbation', 'cell_count', 'time0_tag')]
df_t0.to_csv('../../INPUT/toy_example_input3_time0.tsv', '\t', float_format='%.4g')
df_t0 = df_t0.groupby(['cell_line', 'perturbation', 'time0_tag']).aggregate(trimmean)
df_t0.to_csv('../../INPUT/toy_example_input2_time0.tsv', '\t', float_format='%.4g')
 
df_ctrl = df_data.loc[(df_data.agent=='-') & (df_data.time>0), ('cell_line', 'perturbation', 'replicate', 'time', 'cell_count', 'time0_tag', 'ctrl_tag')]
df_ctrl.to_csv('../../INPUT/toy_example_input3_ctrl.tsv', '\t', float_format='%.4g')
df_ctrl = df_ctrl.groupby(['cell_line', 'perturbation', 'replicate', 'time', 'time0_tag', 'ctrl_tag']).aggregate(trimmean)
df_ctrl.to_csv('../../INPUT/toy_example_input2_ctrl.tsv', '\t', float_format='%.4g')
        
# print out the case 4 example (all in a file)
df_4 = df_data.loc[:, :'cell_count']
df_4.to_csv('../../INPUT/toy_example_input4.tsv', '\t', float_format='%.4g', index=False)

# print out the output (already assigned controls)
df_out = pd.concat((df_data.loc[df_data.agent != '-', :'cell_count'], df_data.loc[df_data.agent != '-', ('GRvalue',)]), axis=1)
df_out.to_csv('../../OUTPUT/toy_example_output.tsv', '\t', float_format='%.4g', index=False)

df_drugs.to_csv('../../OUTPUT/toy_example_DrugParameters.tsv', '\t', float_format='%.4g', index=False)