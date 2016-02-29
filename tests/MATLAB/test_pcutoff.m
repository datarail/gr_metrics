
% load some data
toy1_input = '../../INPUT/toy_example_input1.tsv';
t_test = tsv2table(toy1_input);
t_test = t_test(strcmp(t_test.cell_line,'BT20') & ismember(t_test.agent, {'drugA' 'drugB'}),:);

table2tsv(t_test,'temp.tsv')
%% evaluate the GR value for the data
[t_GRvalues, t_GRmetrics] = GRmetrics([], 'temp.tsv');
[t_GRvalues_1, t_GRmetrics_1] = GRmetrics([], 'temp.tsv', 'pcutoff',2);

delete temp.tsv
%%
assert(all(t_GRmetrics.GR50==t_GRmetrics_1.GR50))
assert(all(t_GRmetrics.GRmax==t_GRmetrics_1.GRmax))
assert(all(t_GRmetrics.EC50==t_GRmetrics_1.EC50 | t_GRmetrics.agent=='drugB'))
assert(all(t_GRmetrics.EC50~=t_GRmetrics_1.EC50 | t_GRmetrics.agent=='drugA'))

assert(all(t_GRmetrics.Hill(t_GRmetrics.agent=='drugB')==.01))
assert(all(t_GRmetrics_1.Hill(t_GRmetrics.Hill==.01)~=.01))
assert(all(t_GRmetrics.EC50==0 | t_GRmetrics.agent=='drugA'))
