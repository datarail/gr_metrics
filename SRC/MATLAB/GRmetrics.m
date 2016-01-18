function [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, input_ctrl, input_time0)
% [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, input_ctrl, input_time0)
% input variables:
%   - output:       folder and tag for the output files (empty means no 
%                       output file will be written)
%   - input_data:   file name for the data on which to compute GR values
%   - input_ctrl:   file name for the control data (optional, see below)
%   - input_time0:  file name for the time0 data (optional, see below)
%
% input files:
%   the files are tab-separated (.tsv) files with the following columns:
%       - 'concentration'
%       - 'cell_count'
%       - 'time' (mandatory only for the input type 'C', see below)
%       - controls values (time0 and untreated control); see below
%       - any other column will be considered a key on which the data will
%           be grouped. E.g. 'cellline', 'drug', 'time', 'replicate', ...
%
% input type for controls:
%   A) one file with the controls matching each treated measure. The
%       untreated controls should be in a column 'cell_count__ctrl'; the
%       time 0 data should be in a column 'cell_count__time0'
%   B) one file ('input_data') with the all treated measures and a key to 
%       match the columns. The keys in the 'input_data' files should be: 
%       'ctrl_tag' for untreated controls and 'time0_tag' for time 0 data.
%       The files 'input_ctrl' and 'input_time0' should each contain a
%       column with 'cell_count' and the key 'ctrl_tag' and 'time0_tag',
%       respectively. Multiple measures for the same key will be averaged
%       (50%-trimmed mean).
%   C) one file will with the controls and all treated measures. The
%       controls will be automatically matched to the treated measures
%       based on the keys. The untreated controls must have
%       'concentration=0' and 'time' matching the treated measures. The
%       time 0 data must have 'time=0'.
%
% output variables are tables:
%   - t_GRvalues contains the GR values for all treated measures.
%   - t_GRmetrics contains the results of the sigmoidal fit for all set of 
%       keys found in the treated measures. The columns of the t_GRmetrics
%       are the keys and the fitted parameters and values: 
%       'GR50' 'GRmax' 'GR_AUC' 'EC50' 'GRinf' 'Hill' 'r2' 'pval'
%       obtained from the sigmoidal fit. 
%   NOTE: The quality of the fit is tested against a flat fit and if the 
%       sigmoidal fit is not significant (p>0.05), the flat fit is prefered.
%
%
%% load the data
t_data = readtable(input_data,'filetype','text','delimiter','\t');

if any(~ismember({'cell_count__ctrl' 'cell_count__time0'}, t_data.Properties.VariableNames))
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
t_GRvalues = sortrows(t_GRvalues);

%% calculate the GR metrics
t_GRmetrics = evaluate_GRmetrics(t_GRvalues);
t_GRmetrics = sortrows(t_GRmetrics);

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
