addpath('../../SRC/MATLAB/')

%% case where the controls are in separate files
toy3_input_data = '../../INPUT/toy_example_input3_data.tsv';
toy3_input_ctrl = '../../INPUT/toy_example_input3_ctrl.tsv';
toy3_input_time0 = '../../INPUT/toy_example_input3_time0.tsv';

output_tag = '../../OUTPUT/MATLAB_results3';

%% evaluate the GR value for the data
[t_GRvalues, t_GRmetrics] = GRmetrics(output_tag, toy3_input_data, ...
    toy3_input_ctrl, toy3_input_time0);

%% plotting the results

plot_GR_results(t_GRvalues, t_GRmetrics);