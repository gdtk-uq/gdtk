========================================================
 Species database description
========================================================
:Author: Rowan J. Gollan
:Date: 2016-01-19


Overview
--------
The species data is stored in separate files. The names of
the files correspond to the names of the compounds for
simple species names. For example, the data for diatomic
oxygen (O2) is found in the file O2.lua. For more complex
compounds, such as ions, the file name is more verbose.
The compound O2+ for example has its data stored in the
file O2_plus.lua. This is to avoid problems with handling
"+" and "-" characters in the file names. The translation
between filenames and actual compound name is handled in
the file species-list.txt. The primary data is stored as
plain text in Lua tables that allow easy manipulation by 
Lua utilities.

Adding a new species
--------------------

The steps to create a new species are as follows.

1. Create a new species file: species_name.lua
   1a. Use an existing file as a template and fill out
       the appropriate data.
   1b. Any missing data will default to the values found in
       defaults.lua.

2. Add the new species entry to species-list.txt.
   2a. Follow the colon-separated format in species-list.txt.
   2b. Add your new entry in alphabetical order.

3. Update the species-database.lua by running the makefile. 
   3a. > make




