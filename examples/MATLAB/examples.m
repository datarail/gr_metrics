addpath('../../SRC/MATLAB/')

%% reference example where controls are assgined
toy1_input = '../../INPUT/toy_example_input1.tsv';
% evaluate the GR value for the data
[t_GRvalues_1, t_GRmetrics_1] = GRmetrics([], toy1_input);

%% case where the controls are in separate files
toy3_input_data = '../../INPUT/toy_example_input3_data.tsv';
toy3_input_ctrl = '../../INPUT/toy_example_input3_ctrl.tsv';
toy3_input_time0 = '../../INPUT/toy_example_input3_time0.tsv';

output_tag = '../../OUTPUT/MATLAB_results_input3.tsv';

% evaluate the GR value for the data
[t_GRvalues_3, t_GRmetrics_3] = GRmetrics(output_tag, toy3_input_data, toy3_input_ctrl, toy3_input_time0);

%% plot the results

t_GRvalues_ref = readtable('../../OUTPUT/toy_example_output.tsv', ...
    'filetype','text','delimiter','\t');

t_GRvalues = t_GRvalues_1;
t_GRmetrics = t_GRmetrics_1;
cell_lines = unique(t_GRvalues.cell_line);
xc = 10.^(-3.5:.01:1.5);
colors = [.8 .1 0; .7 1 .1; .2 0 .8];
for d = 'A':'D'
    figure(d-64);clf;
    set(gcf,'position', [(d-64)*30 (d-64)*30 800 700])
    
    drug = ['drug' d];
for ip = [0 1]
    for it = [48 72]
        for ir = 1:3
            subplot(3,4, (ir-1)*4+ip*2+(it-24)/24)
            title(sprintf('t=%i, perturb=%i', it, ip))
            hold on
            h = [];
            for iC = 1:length(cell_lines)
                t_ = t_GRvalues_ref(t_GRvalues_ref.perturbation==ip & ...
                    t_GRvalues_ref.time==it & t_GRvalues_ref.replicate==ir & ...
                    strcmp(t_GRvalues_ref.cell_line, cell_lines{iC}) & ...
                    strcmp(t_GRvalues_ref.agent,drug),:);
                if isempty(t_), continue, end
                h(iC) = plot(log10(t_.concentration), t_.GRvalue, 'x', 'color', ...
                    colors(iC,:));
                
                t_ = t_GRvalues(t_GRvalues.perturbation==ip & ...
                    t_GRvalues.time==it & t_GRvalues.replicate==ir & ...
                    strcmp(t_GRvalues.cell_line, cell_lines{iC}) & ...
                    strcmp(t_GRvalues.agent,drug),:);
                if isempty(t_), continue, end
                h(iC) = plot(log10(t_.concentration), t_.GRvalue, 'o', 'color', ...
                    colors(iC,:));
                
                
                t_ = t_GRmetrics(t_GRmetrics.perturbation==ip & ...
                    t_GRmetrics.time==it & t_GRmetrics.replicate==ir & ...
                    t_GRmetrics.cell_line==cell_lines{iC} & ...
                    t_GRmetrics.agent==drug,:);
                
                plot(log10(xc), t_.GRinf + (1-t_.GRinf)./ ...
                    ( 1 + (xc/t_.EC50).^t_.Hill), '-', 'color', ...
                    colors(iC,:))
            end
            xlim([-3.5 1.5])
            ylim([-1 1.3])
            
            if ((ir-1)*4+ip*2+(it-24)/24)==1
                hl = legend(h, cell_lines, 'orientation', 'horizontal');
                set(hl,'position',[.1 .01 .8 .04])
            end
        end
    end
end
end