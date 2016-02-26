addpath('../../SRC/MATLAB/')

headers = {'concentration' 'cell_count' 'cell_count__ctrl' 'cell_count__time0'};

merged_data = [ 10.^[-3:.5:1]' [1;1;1;.95;.9;.6;.3;.25;.2] ones(9,1) .3*ones(9,1)];
t_test = [table(repmat({'A'},9,1), 'variablenames', {'agent'}) ...
    array2table(merged_data, 'variablenames', headers)];

shift = [0 .05 .05 .05];
test_data_1 = [
    merged_data+repmat(shift,9,1);
    merged_data-repmat(shift,9,1)];
t_test_1 = [
    table([repmat({'a'},9,1);repmat({'b'},9,1)], 'variablenames', {'condition'}), ...
    repmat(t_test(:,'agent'),2,1) ...
    array2table(test_data_1, 'variablenames', headers)];

shift = [0 .025 .025 .025];
test_data_2 = [
    test_data_1+repmat(shift,18,1);
    test_data_1-repmat(shift,18,1)];
t_test_2 = [
    table([ones(18,1);2*ones(18,1)], 'variablenames', {'replicate'}), ...
    repmat(t_test_1(:, {'condition' 'agent'}),2,1) ...
    array2table(test_data_2, 'variablenames', headers)];


%%
table2tsv(t_test,'temp.tsv')
t_GRv = GRmetrics([], 'temp.tsv');

table2tsv(t_test_1,'temp.tsv')
t_GRv_1 = GRmetrics([], 'temp.tsv');
t_GRv_a = GRmetrics([], 'temp.tsv', 'collapseKey', 'condition');

table2tsv(t_test_2,'temp.tsv')
t_GRv_1b = GRmetrics([], 'temp.tsv', 'collapseKey', 'replicate');
t_GRv_b = GRmetrics([], 'temp.tsv', 'collapseKey', {'replicate' 'condition'});

delete temp.tsv

%%
assert(all(abs(t_GRv.GRvalue-t_GRv_a.GRvalue)<1e-10))
assert(all(abs(t_GRv.GRvalue-t_GRv_b.GRvalue)<1e-10))
assert(all(abs(t_GRv_1.GRvalue-t_GRv_1b.GRvalue)<1e-10))
