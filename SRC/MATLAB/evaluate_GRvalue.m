function t_out = evaluate_GRvalue(t_in)
% t_out = evaluate_GRvalue(t_in)
%   evalute GRvalue on a table with columns: cell_count, cell_count__time0 and cell_count__ctrl

% --> change to a real error handling MH 16/1/21
assert(all(ismember({'cell_count', 'cell_count__time0' 'cell_count__ctrl'}, ...
    t_in.Properties.VariableNames)), ...
    'Need the columns ''cell_count'', ''cell_count__time0'', ''cell_count__ctrl'' in the data')

t_out = t_in;
t_out.GRvalue = 2.^( log2(t_in.cell_count./t_in.cell_count__time0) ./ ...
    log2(t_in.cell_count__ctrl./t_in.cell_count__time0) ) -1;

end
