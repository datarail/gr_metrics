import sys
import textwrap
import fileinput
import pandas as pd
import gr50
from gr50.linereader import LineReader

def main():
    '''
    Usage: compute_gr_metrics.py input.tsv > output.tsv

    Compute Growth Response metrics for precomputed GR dose-response data in
    TSV format.

    The input tsv file must meet the following requirements:

    - The first row contains column names.
    - The 'concentration' column contains the numeric dose of a perturbagen.
    - The 'GRvalue' column contains the GR value for that condition.

    The input file may also have other "key" columns to distinguish multiple
    separate dose-response experiments within the same file.

    The GR value may be computed from raw cell counts using the
    `add_gr_column.py` script. The columns required as input for that script
    will be ignored by this one, and as such the output of `add_gr_column.py`
    may be piped directly to the input of this script.

    The output tsv will have all key columns in addition to a column for each GR
    metric. (See the documentation of gr50.gr_metrics for details)
    '''

    if '-h' in sys.argv or '--help' in sys.argv:
        print textwrap.dedent(main.__doc__)
        return
    reader = LineReader(fileinput.input(mode='rb'))
    data = pd.read_csv(reader, delimiter='\t')
    metrics = gr50.gr_metrics(data)
    print metrics.to_csv(sep='\t', index=False)

if __name__ == '__main__':
    main()
