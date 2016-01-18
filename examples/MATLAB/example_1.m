addpath('../../SRC/MATLAB/')

%% case where the controls are in separate files
toy1_input_data = '../../INPUT/toy_example_input1.tsv';

%% evaluate the GR value for the data
[t_GRvalues, t_GRmetrics] = GRmetrics([], toy1_input_data);

%% plotting the results

plot_GR_results(t_GRvalues, t_GRmetrics);