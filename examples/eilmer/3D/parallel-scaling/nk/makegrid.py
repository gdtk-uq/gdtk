"""
Write out an su2 format grid for the timing tests

@author: Nick Gibbons
"""

from subprocess import run

run('e4shared --custom-script --script-file=gengrid.lua'.split(), check=True)
