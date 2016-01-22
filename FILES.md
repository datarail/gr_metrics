* name: GR50_tools
* start date: Fri Nov 20 12:05:53 EST 2015
* description
* files
** GR_tools_tutorial.pdf
** GR_tools_tutorial.md
** input/
    `toy_example_*.tsv
        generated with /examples/generate_data/generate_data.py
** output/
    `toy_example_*.tsv
        generated with /examples/generate_data/generate_data.py

** NOTES.org

* software
** SRC/python/add_gr_column.py

    Usage:
    add_gr_column.py input.tsv > output.tsv

    input.tsv must (case A)
    - have column headers in the first row
    - have keys on the first columns
    - contain numeric columns with headers cell_count, cell_count__ctrl, and
      cell_count__time0

**  SRC/MATLAB/GRmetrics.m

    Usage: 
    [t_GRvalues, t_GRmetrics] = GRmetrics(output, input_data, input_ctrl, input_time0)
    
    all input variables are file names. See help of function for details on the type 
    of inputs and the output format (MATLAB table)
    
* tests
** MATLAB/*
    test only the consistency of the different cases. No error testing implemented yet