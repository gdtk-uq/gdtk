#!/usr/bin/env python
# file: fpdiff.py
"""
Compare two files, looking significant differences in numerical values.

Can be used as a stand-alone script or can be imported into a Python
program for more extensive application.

PJ, 28-Sep-02
    23-Jul-2018 Updated to Python3
"""

from __future__ import print_function
import sys, os

def count_differences(filenameA, filenameB, tolerance=1.0e-9, do_print=0):
    """
    Given the file names, count the differences in their numerical data.

    Optionally set the tolerance for testing differences and allow printing
    of the lines with differences.
    """
    fA = open(filenameA, 'r')
    fB = open(filenameB, 'r')

    line_count = 0
    difference_count = 0
    while 1:
        lineA = fA.readline()
        lineB = fB.readline()
        if lineA == '' or lineB == '':
            break
        line_count += 1
        lineA = lineA.strip()
        lineB = lineB.strip()
        numbersA = lineA.split()
        numbersB = lineB.split()
        linediff = 0
        lenA = len(numbersA)
        lenB = len(numbersB)
        if lenA != lenB:
            linediff += 1
        for i in range( min(lenA, lenB) ):
            itemA = numbersA[i]
            itemB = numbersB[i]
            # print '(', itemA, itemB, ')'
            # if itemA[0] == '#' or itemB[0] == '#':
            #     break
            try:
                a = float(itemA)
                b = float(itemB)
            except ValueError as e:
                break
            if abs(a - b)/(0.5 * (abs(a) + abs(b)) + 1.0) > tolerance :
                linediff += 1
        if linediff > 0 and do_print == 1:
            print('Difference at line', line_count)
            print('fileA:', lineA)
            print('fileB:', lineB)
        difference_count += linediff
    return difference_count


def print_usage_message():
    print('Usage: fpdiff <fileA> <fileB> ?tol=<tolerance>? ?summary?')

#-----------------------------------------------------------------------------    
if __name__ == '__main__':
    # We need at least the two file names to start work.
    if len(sys.argv) < 3:
        print_usage_message()
        sys.exit()

    filenameA = sys.argv[1]
    filenameB = sys.argv[2]
    print_summary = 0
    tolerance = 1.0e-9
    
    if os.path.isfile(filenameA) and os.path.isfile(filenameB):
        # Pull apart the rest of the command-line arguments.
        rest_of_argv = sys.argv[3:]
        for arg in rest_of_argv:
            if arg.startswith("tol") or arg.startswith("-tol"):
                try:
                    tolerance = float(arg.split('=')[1])
                except:
                    tolerance = 1.0e-9
            if arg.startswith("sum") or arg.startswith("-sum"):
                print_summary = 1
                
        diff_count = count_differences(filenameA, filenameB, tolerance, 1)
        if print_summary:
            print('Found', diff_count, 'differences with tolerance', tolerance)
    else:
        print('One file (at least) is invalid.')
        print_usage_message()
        
