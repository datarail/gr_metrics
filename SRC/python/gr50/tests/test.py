import sys, os, os.path
import pandas as pd
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import gr50
import numpy as np

def output_test():
    curfilePath = os.path.abspath(__file__)
    Dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(curfilePath)))))
    Dirfile = os.path.abspath(os.path.join(Dir,'INPUT/toy_example_input1.tsv'))
    df = pd.read_csv(Dirfile, delimiter='\t')
    # Read the matlab output.
    Dirfile = os.path.abspath(os.path.join(Dir,'OUTPUT/matlab_input1_GRmetrics.tsv'))
    metrics_matlab = pd.read_csv(Dirfile, delimiter='\t')

    # Compute the GR metrics from the data.
    gr_values = gr50.compute_gr(df)
    metrics_python = gr50.gr_metrics(gr_values)

    print(metrics_python.columns)
    print(metrics_matlab.columns)

    # Determine key columns and sort both dataframes by those columns. Since
    # this script is expected to work on arbitrary inputs, we need to "sniff"
    # out which columns are the keys (otherwise we'd need that as an input too).
    # Both the matlab and python implementations put the keys first and the
    # metrics second (starting with GR50) so we will use that to our advantage.
    # Of course they don't both sort the rows the same way which is why we need
    # to do that here.
    first_data_col_index = list(metrics_python).index('GR50')
    keys_python = metrics_python.columns[:first_data_col_index]
    keys_matlab = metrics_matlab.columns[:first_data_col_index]
    assert sorted(keys_python) == sorted(keys_matlab), "Key column mismatch"
    keys = list(keys_python)
    for df in metrics_python, metrics_matlab:
        # Perform the sort.
        df.sort_values(keys, inplace=True)
        df.reset_index(drop=True, inplace=True)
        # Compute log10 of GR50 and drop original column.
        df.insert(first_data_col_index, 'log10_GR50', np.log10(df['GR50']))
        del df['GR50']
        df.insert(first_data_col_index, 'log10_GEC50', np.log10(df['GEC50']))
        del df['GEC50']

    metrics = [c for c in metrics_python if c not in keys]
    errs = ['err_' + m for m in metrics]
    # ignore pval and r^2 for now
    errs = errs[:-2]
    error = metrics_python[keys].copy()
    for m, e in zip(metrics, errs):
        err = (1 - metrics_python[m] / metrics_matlab[m])
        error[e] = err.fillna(0)
    rejects = np.isinf(metrics_python.log10_GR50)
    rejects2 = np.isinf(metrics_python.log10_GEC50)
    # Sanity check to make sure both versions pass/reject the same things.
    assert (rejects == np.isinf(metrics_matlab.log10_GR50)).all()
    assert (rejects2 == np.isinf(metrics_matlab.log10_GEC50)).all()

    # Keep pandas from wrapping.
    pd.set_option('display.width', None)

    print
    print ("Records with > 0.1% error between python and matlab values of any "
           "metric:")
    test1 = error[(abs(error[errs]) > 0.001).any(axis=1)]
    print(test1)

    print
    print ("Records with > 10% error between python and matlab values of any "
           "metric:")
    test2 = error[(abs(error[errs]) > 0.1).any(axis=1)]
    print(test2)

    # Insert rows at the top with median and IQR of each error column.
    median = error[errs].median()
    iqr = np.subtract(*np.percentile(error[errs], [75, 25], axis=0))
    error.loc[-2] = ['<median>', -1, -1, 'N/A', -1] + list(median)
    error.loc[-1] = ['<IQR>', -1, -1, 'N/A', -1] + list(iqr)
    error.sort_index(inplace=True)
    # Passes if none of the metrics differ by 10%
    # Just an example of a passing test.
    assert test2.shape[0] == 0
