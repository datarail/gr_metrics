#!/usr/bin/env python

import csv
import fileinput as fi
import collections as co
import math as ma
import decimal

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

    add_gr_column_0.py input.tsv > output.tsv

    input.tsv must

    - have column headers in the first row
    - contain numeric columns with headers cell_count, cell_count__ctrl, and
      cell_count__time0

    '''
    reader = csv.reader(fi.input(mode='rb'), delimiter='\t')

    headers = next(reader)
    to_record = co.namedtuple(typename='record',
                              field_names=headers,
                              rename=True)

    def adjust_sigfigs(val, _context=decimal.Context(prec=3)):
        return _context.create_decimal(val)

    for r in (to_record(*row) for row in reader):
        cell_count__time0 = float(r.cell_count__time0)
        log2nn, log2nn_ctrl = (normalize_log2(float(n), cell_count__time0)
                               for n in (r.cell_count, r.cell_count__ctrl))

        gr = 2**(log2nn/log2nn_ctrl) - 1

        print_augmented_row(row=r, last_col=adjust_sigfigs(gr))

main()
