addpath('../../SRC/MATLAB/')

t_data = readtable('../../INPUT/toy_example_input1.tsv','filetype','text','delimiter','\t');
t_GRvalues = evaluate_GRvalue(t_data(strcmp(t_data.agent,'drugA'),2:end));

%%
t_GRmetrics = evaluate_GRmetrics(t_GRvalues);


%%
figure(1);clf
cell_lines = unique(t_GRvalues.cell_line);
xc = 10.^(-3.5:.01:1.5);
for ip = [0 1]
    for it = [48 72]
        for ir = 1:3
            subplot(3,4, (ir-1)*4+ip*2+(it-24)/24)
            hold on
            for iC = 1:length(cell_lines)
                t_ = t_GRvalues(t_GRvalues.perturbation==ip & ...
                    t_GRvalues.time==it & t_GRvalues.replicate==ir & ...
                    strcmp(t_GRvalues.cell_line, cell_lines{iC}) & ...
                    strcmp(t_GRvalues.agent,'drugA'),:);
                plot(log10(t_.concentration), t_.GRvalue, 'o', 'color', ...
                    [iC/length(cell_lines) 0 1-(iC/length(cell_lines))])
                
                
                t_ = t_GRmetrics(t_GRmetrics.perturbation==ip & ...
                    t_GRmetrics.time==it & t_GRmetrics.replicate==ir & ...
                    t_GRmetrics.cell_line==cell_lines{iC} & ...
                    t_GRmetrics.agent=='drugA',:);
                
                plot(log10(xc), t_.GRinf + (1-t_.GRinf)./ ...
                    ( 1 + (xc/t_.EC50).^t_.Hill), '-', 'color', ...
                    [iC/length(cell_lines) 0 1-(iC/length(cell_lines))])
            end
            xlim([-3.5 1.5])
            ylim([-1 1])
        end
    end
end