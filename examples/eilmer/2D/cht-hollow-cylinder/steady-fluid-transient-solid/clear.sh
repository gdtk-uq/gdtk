#!/bin/bash

# subdirectories
rm -r config flow grid hist loads plot solid solid-grid residuals limiter-values CellData

# gas/chemistry model files
rm *.gas *.chem *.exch 

# solver files
rm *saved e4-nk.diagnostics.dat *pdf

# general tmp files
rm eilmer_solution.dat log* *e4-nk.diagnostics.dat *~
