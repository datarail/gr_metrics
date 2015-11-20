#!/usr/bin/env python

import csv
import fileinput as fi
import collections as co
import math as ma

# ------------------------------------------------------------------------------
log2 = lambda x: ma.log(x, 2)

def normalize_log2(n, n_0_0):
    assert isinstance(n, float) or isinstance(n_0_0, float)
    normalized = n/n_0_0
    return log2(normalized)

def print_augmented_row(row, last_col):
    print '%s\t%s' % ('\t'.join(row), last_col)

# ------------------------------------------------------------------------------

def main():
    '''
    Usage:

    add_ngri_column_0.py input.tsv > output.tsv

    input.tsv must

    - have column headers in the first row
    - contain numeric columns with headers n_t_c, n_t_0, n_0_0

    '''
    reader = csv.reader(fi.input(mode='rb'), delimiter='\t')

    headers = next(reader)
    to_record = co.namedtuple(typename='record',
                              field_names=headers,
                              rename=True)

    print_augmented_row(row=headers, last_col='ngri')

    for r in (to_record(*row) for row in reader):

        n_0_0 = float(r.n_0_0)
        log2nn, log2nn_ctrl = (normalize_log2(float(n), n_0_0)
                               for n in (r.n_t_c, r.n_t_0))

        # ngri: normalized growth-rate inhibition
        ngri = '%.6g' % (2**(log2nn/log2nn_ctrl) - 1)

        print_augmented_row(row=r, last_col=ngri)

main()
