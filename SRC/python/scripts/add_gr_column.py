import sys
import textwrap
import fileinput
import csv
import collections
import decimal
import gr50

# ------------------------------------------------------------------------------

def print_augmented_row(row, last_col):
    print '%s\t%s' % ('\t'.join(row), last_col)

# ------------------------------------------------------------------------------

def main():
    '''
    Usage: add_gr_column.py input.tsv > output.tsv

    Compute Growth Response value from cell count data in TSV format.

    The input tsv file must meet the following requirements:

    - The first row contains column names.
    - The 'cell_count' column contains measured cell counts under a given
      perturbation at some time.
    - The 'cell_count__time0' column contains cell counts for the corresponding
      time=0 controls.
    - The 'cell_count__ctrl' column contains cell counts for the corresponding
      no-perturbation controls.

    All other columns are ignored and passed through untouched.

    The output tsv will have all columns from the input, as well as a new
    'GRvalue' column at the end containing the GR values.
    '''

    if '-h' in sys.argv or '--help' in sys.argv:
        print textwrap.dedent(main.__doc__)
        return

    reader = csv.reader(fileinput.input(mode='rb'), delimiter='\t')
    headers = next(reader)
    print_augmented_row(row=headers, last_col='GRvalue')

    record = collections.namedtuple(typename='record',
                                    field_names=headers,
                                    rename=True)
    for r in (record(*row) for row in reader):
        gr = gr50.compute_gr_single(r)
        print_augmented_row(row=r, last_col=gr)

if __name__ == '__main__':
    main()
