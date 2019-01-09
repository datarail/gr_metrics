function [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, varargin)
% [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, [input_ctrl, input_time0], [varargin])
%
% input variables:
%   - output:       folder and tag for the output files (empty means no
%                       output file will be written)
%   - input_data:   file name for the data on which to compute GR values
%   - input_ctrl:   file name for the control data (optional, see below)
%   - input_time0:  file name for the time0 data (optional, see below)
%
% optional input parameters (property/value pairs):
%   - 'pcutoff':    cutoff value for the flat fit using the F-test
%                       (default is pcutoff=0.05)
%   - 'collapseKeys': column header (key) to average the data on
%                       (default is none)
%
% input files:
%   the files are tab-separated (.tsv) files with the following columns:
%       - 'concentration'
%       - 'cell_count'
%       - 'time' (mandatory only for the input type 'C', see below)
%       - control values (time0 and untreated control); see below
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
%    NOTE: for the case C), some assumptions are made regarding the keys
%       and how controls are assigned. If the structure of the data is
%       complex, it is better the use cases A) or B)
%
% output variables are tables:
%   - t_GRvalues contains the GR values for all treated measures.
%   - t_GRmetrics contains the results of the sigmoidal fit for all set of
%       keys found in the treated measures. The columns of the t_GRmetrics
%       are the keys and the fitted parameters and values:
%       'GR50' 'GRmax' 'GR_AOC' 'GEC50' 'GRinf' 'h_GR' 'r2_GR' 'pval_GR'
%       obtained from the sigmoidal fit.
%   NOTE: The quality of the fit is tested against a flat fit and if the
%       sigmoidal fit is not significant (p>0.05), the flat fit is prefered.
%
%

% further improvements to implement (MH 1/18/16):
%   - check on the input data and output descriptive error messages:
%       - cases B/C: all treated measures have corresponding controls
%       - case C: input for time 0 have concentration=0
%       - case C: input for untreated controls have concentration=0
%   - assess quality of the data
%       - variablity of the controls
%
% code improvements to implement (MH 2/25/16):
%   - proper handling of errors and tests (remove assert functions)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% assign the optional parameters
assert(mod(nargin,2)==0, 'Expecting a even number of arguments - see help')

caseB = false;
if length(varargin)>=2
    % handle the case B (overloading of input arguments)
    if exist(varargin{1}, 'file')==2 && exist(varargin{2}, 'file')==2
        input_ctrl = varargin{1};
        input_time0 = varargin{2};
        varargin = varargin(3:end);
        caseB = true;
    end
end

% get the parameters from the input arguments
p = inputParser;
addParameter(p,'pcutoff', .05, @isnumeric);
addParameter(p,'collapseKeys', '', @(x) ischar(x) || iscellstr(x));
parse(p,varargin{:});
p = p.Results;

%% load the data
t_data = readtable(input_data,'filetype','text','delimiter','\t');

% --> change to a real error handling MH 16/1/21
assert(all(ismember({'cell_count', 'concentration'}, t_data.Properties.VariableNames)), ...
    'Need the columns ''cell_count'', ''concentration'' in the data')

if any(~ismember({'cell_count__ctrl' 'cell_count__time0'}, t_data.Properties.VariableNames))
    % need to assign the controls to each condition

    if caseB
        % case of multiple input files

        % endpoint controls
        t_ctrl = readtable(input_ctrl,'filetype','text','delimiter','\t');
        % time0 controls
        t_time0 = readtable(input_time0,'filetype','text','delimiter','\t');

        % match the controls with the data
        t_data = add_controls(t_data, t_ctrl, t_time0);

    else
        % case of one long with with the controls:
        %   compute the controls and assign them to the measured data
        t_data = assign_controls(t_data);
    end

end

%% changed to averaging GR values after calculation instead of
%% averaging cell counts before GR value calculation - NC 07/20/2017
%%% collapse the data (if neeeded)
%if ~isempty(p.collapseKeys)
%    t_data = collapsed_data(t_data, p.collapseKeys);
%end

%% evaluate the GR value for the data
t_GRvalues = evaluate_GRvalue(t_data);
t_GRvalues = sortrows(t_GRvalues); % for consistent output

%% calculate the GR metrics
if ~isempty(output) || nargout>1
    % skipped if not needed for output
    t_GRmetrics = evaluate_GRmetrics(t_GRvalues, p.pcutoff, p.collapseKeys);
    t_GRmetrics = sortrows(t_GRmetrics); % for consistent output
end
%% write the output files
if ~isempty(output)
    writetable(t_GRvalues, strcat(output, '_GRvalues.tsv'), ...
        'filetype','text' , 'delimiter', '\t');
    writetable(t_GRmetrics, strcat(output, '_GRmetrics.tsv'), ...
        'filetype','text' , 'delimiter', '\t');
end

if nargout==0
    clear t_GRvalues t_GRmetrics
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> write the subfunctions in separate files MH 16/1/21
function t_data = assign_controls(t_all)
% split the data in measured values and controls based on keys
%   and assign the proper tags placing the controls right

keys = setdiff(t_all.Properties.VariableNames, {'concentration' 'cell_count'}, ...
    'stable');

% find the time 0:
t_time0 = t_all(t_all.time==0, [keys, 'cell_count']);

% find the untreated controls:
t_ctrl = t_all(t_all.concentration==0 & t_all.time>0, [keys, 'cell_count']);

% filter the data
t_data = t_all(t_all.concentration~=0 & t_all.time>0,:);

% find which columns are keys for the different controls (to exclude
% treatements)
time0_keys = false(1, length(keys));
ctrl_keys = false(1, length(keys));
for i = 1:length(keys)
    time0_keys(i) = ~isempty(intersect(t_time0.(keys{i}), t_data.(keys{i})));
    ctrl_keys(i) = ~isempty(intersect(t_ctrl.(keys{i}), t_data.(keys{i})));
end

% compute the trimmed mean for time0 data
t_time0_mean = grpstats(t_time0(:,[keys(time0_keys) 'cell_count']), ...
    keys(time0_keys),@(x)trimmean(x,50));
% clean the table
t_time0_mean.GroupCount = [];
t_time0_mean.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__time0';

% compute the trimmed mean for untreated data
t_ctrl_mean = grpstats(t_ctrl(:,[keys(ctrl_keys) 'cell_count']), ...
    keys(ctrl_keys),@(x)trimmean(x,50));
% clean the table
t_ctrl_mean.GroupCount = [];
t_ctrl_mean.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__ctrl';

% join all the tables
% --- > need to check for proper join (no empty rows)  MH 16/1/21
t_data = outerjoin(t_data, t_time0_mean, 'type', 'left', ...
    'MergeKeys', true, 'keys', keys(time0_keys));
t_data = outerjoin(t_data, t_ctrl_mean, 'type', 'left', ...
    'MergeKeys', true, 'keys', keys(ctrl_keys));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> write the subfunctions in separate files MH 16/1/21

function t_data = add_controls(t_data, t_ctrl, t_time0)
% join the tables of controls with the table of data

% endpoint controls (clean and perform trimmed mean)
t_ctrl = t_ctrl(:, {'cell_count' 'ctrl_tag'});
t_ctrl = grpstats(t_ctrl,'ctrl_tag',@(x)trimmean(x,50));
t_ctrl.GroupCount = [];
t_ctrl.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__ctrl';

% time0 controls (clean and perform trimmed mean)
t_time0 = t_time0(:, {'cell_count' 'time0_tag'});
t_time0 = grpstats(t_time0,'time0_tag',@(x)trimmean(x,50));
t_time0.GroupCount = [];
t_time0.Properties.VariableNames{'Fun1_cell_count'} = 'cell_count__time0';

% merge the controls with the measured data
t_data = outerjoin(t_data, t_ctrl, 'type', 'left', ...
    'MergeKeys', true, 'keys', 'ctrl_tag');
t_data = outerjoin(t_data, t_time0, 'type', 'left', ...
    'MergeKeys', true, 'keys', 'time0_tag');

% clean the flags
t_data.ctrl_tag = [];
t_data.time0_tag = [];
end


% changed to averaging GR values after calculation instead of
% averaging cell counts before calculating GR values - NC 07/20/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function t_data_collapsed = collapsed_data(t_data, collapseKeys)
%% average the results (cell count) by removing one key
%if ischar(collapseKeys), collapseKeys = {collapseKeys}; end
%assert(all(ismember(collapseKeys, t_data.Properties.VariableNames)))
%
%cell_count_columns = {'cell_count' 'cell_count__ctrl' 'cell_count__time0'};
%
%[t_data_collapsed, ~, idx] = unique(t_data(:,setdiff(t_data.Properties.VariableNames, ...
%    [collapseKeys cell_count_columns], 'stable')), 'stable');
%
%t_data_collapsed = [t_data_collapsed array2table(NaN(height(t_data_collapsed),3), ...
%    'variableNames', cell_count_columns)];
%for i=1:height(t_data_collapsed)
%    for j=1:length(cell_count_columns)
%        t_data_collapsed.(cell_count_columns{j})(i) = ...
%            mean(t_data.(cell_count_columns{j})(idx==i));
%    end
%end
%end
