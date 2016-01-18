* [151120F]
** implemented add_ngri_column_0.py
*** column names are provisional
** to do
*** finalize column names
* [151124T]
** fixed column names
** changed log10c column to concentration
*** changed numeric format to 3 significant figures
* [151125W]
** changed abbreviation ngri to gr throughout
** changed numeric format of gr column to 3 significant figures
** updated PROTOCOL
* [151127F]
** fixed numeric format of gr column
   Previous version only approximated a 3-significant figure format
   using a string formatting template; it worked fine for numbers like
   0.123 but rendered as, e.g., 0.7 what should have been 0.700; the
   current version gets that right.
** eliminated some abbreviations in the code
   ...for the sake of clarity
* [151130M]
** restored line accidentally deleted from add_gr_column_0.py
** reordered columns of input/toy.tsv
   primarily so that the `agent` and `concentration` columns are
   adjacent
* [160111M]
** done
   - fix .gitnit -> .gitignore
   - 17AAG -> 17-AAG in ./**/toy*
   - input/ -> INPUT/, &c
   - root_pipeline -> STANDARDIZATION
   - fix add_gr_column_0.py
     - rename to add_gr_column.py
     - run main conditionally
   - update warmup/PROTOCOL.org
** state
*** tree
    % tree -I .git $HOME/_/projects/gr50_browser
    /home/berriz/_/projects/gr50_browser
    .
    |-- .gitignore
    |-- STANDARDIZATION/
    |   |-- INPUT/
    |   |-- OUTPUT/
    |   |   `-- LOG/
    |   |-- PROTOCOL.org
    |   `-- SRC/
    |       `-- TEST/
    `-- warmup/
        |-- INPUT/
        |   |-- toy.tsv
        |   `-- toy.txt
        |-- NOTES.org
        |-- OUTPUT/
        |   |-- LOG/
        |   `-- toy1.tsv
        |-- PROTOCOL.org
        `-- SRC/
            |-- add_gr_column.py*
            |-- makefake.py
            `-- TEST/
                `-- test_format_sig_figs.py
*** git ls-files
    .gitignore
    warmup/INPUT/toy.tsv
    warmup/NOTES.org
    warmup/OUTPUT/toy1.tsv
    warmup/PROTOCOL.org
    warmup/SRC/TEST/test_format_sig_figs.py
    warmup/SRC/add_gr_column.py

** todo
   - update and commit warmup/SRC/TEST

* [180118M] and prior days
** done
    - MATLAB implementation
    - Python scripts to generate examples
    - tests and examples for the MATLAB scripts
    - documentation (Word tutorial)
    
** state
    major changes on the structure of the repo for publication
    
** to do
    - convert the tutorial in .md (such that it is self-contained in GitHub)
    - additional functionality as described in the MATLAB file GRmetrics.m
    -- priority: testing of the input!
    - implementation of the sigmoidal fit in Python
    - long term: plotting functions