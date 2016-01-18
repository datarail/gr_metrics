function [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, input_ctrl, input_time0)
% [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, input_ctrl, input_time0)


%% load the data
t_data = readtable(input_data,'filetype','text','delimiter','\t');

if any(~ismember({'ctrl_tag' 'time0_tag'}, t_data.Properties.VariableNames))
    % need to assign the controls to each condition
    
    if nargin>2
        % case of multiple input files
        
        % endpoint controls
        t_ctrl = readtable(input_ctrl,'filetype','text','delimiter','\t');
        t_ctrl = t_ctrl(:, {'cell_count' 'ctrl_tag'});
        t_ctrl = grpstats(t_ctrl,'ctrl_tag',@(x)trimmean(x,50));
        t_ctrl.GroupCount = [];
        t_ctrl.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__ctrl';
        
        % time0 controls
        t_time0 = readtable(input_time0,'filetype','text','delimiter','\t');
        t_time0 = t_time0(:, {'cell_count' 'time0_tag'});
        t_time0 = grpstats(t_time0,'time0_tag',@(x)trimmean(x,50));
        t_time0.GroupCount = [];
        t_time0.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__time0';
        
        % merge the controls
        t_data = outerjoin(t_data, t_ctrl, 'type', 'left', ...
            'MergeKeys', true, 'keys', 'ctrl_tag');
        t_data.ctrl_tag = [];
        
        t_data = outerjoin(t_data, t_time0, 'type', 'left', ...
            'MergeKeys', true, 'keys', 'time0_tag');
        t_data.time0_tag = [];
        
    else
        % case of one long with with the controls
        
        
    end
    
end

%% evaluate the GR value for the data
t_GRvalues = evaluate_GRvalue(t_data);


%% calculate the GR metrics
t_GRmetrics = evaluate_GRmetrics(t_GRvalues);


%% write the output files
if ~isempty(output)
    writetable(t_GRvalues, [output '_GRvalues.tsv'], ...
        'filetype','text' , 'delimiter', '\t');
    
    writetable(t_GRmetrics, [output '_GRmetrics.tsv'], ...
        'filetype','text' , 'delimiter', '\t');
end

if nargout==0
    clear t_GRvalues t_GRmetrics
end
