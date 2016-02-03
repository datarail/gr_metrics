import sys
import textwrap
import os.path as path
import numpy as np
import pandas as pd

def main():
    '''
    Usage: compare_python_matlab.py python_metrics.tsv matlab_metrics.tsv

    Compare GR50 values calculated by python and matlab implementations.

    The two input files must have the same key columns and the same sets of
    values in those columns, but the ordering does not have to be the same.
    '''

    if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) != 3:
        print textwrap.dedent(main.__doc__)
        return

    metrics_python = pd.read_csv(sys.argv[1], delimiter='\t')
    metrics_matlab = pd.read_csv(sys.argv[2], delimiter='\t')

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

    metrics = [c for c in metrics_python if c not in keys]
    errs = ['err_' + m for m in metrics]
    error = metrics_python[keys].copy()
    for m, e in zip(metrics, errs):
        err = (1 - metrics_python[m] / metrics_matlab[m])
        error[e] = err.fillna(0)
    rejects = np.isinf(metrics_python.log10_GR50)
    # Sanity check to make sure both versions pass/reject the same things.
    assert (rejects == np.isinf(metrics_matlab.log10_GR50)).all()

    # Keep pandas from wrapping.
    pd.set_option('display.width', None)

    print
    print ("Records with > 0.1% error between python and matlab values of any "
           "metric:")
    print
    print error[(abs(error[errs]) > 0.001).any(axis=1)]

    # Insert rows at the top with median and IQR of each error column.
    median = error[errs].median()
    iqr = np.subtract(*np.percentile(error[errs], [75, 25], axis=0))
    error.loc[-2] = ['<median>', -1, -1, 'N/A', -1] + list(median)
    error.loc[-1] = ['<IQR>', -1, -1, 'N/A', -1] + list(iqr)
    error.sort_index(inplace=True)

    base_path = path.join(path.dirname(path.abspath(__file__)),
                          '..', '..', '..')
    output_path = path.join(base_path, 'OUTPUT', 'python_matlab_comparison.tsv')
    print
    print "Saving full report to %s" % path.relpath(output_path)
    print
    error.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    main()
