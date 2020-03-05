#! /usr/bin/env python
#
# Script to compare two flow solutions
# inputs:
# 1. 1D ordered flow solution
# 2. 1D ordered flow solution
#
# note: input 2 may have more data points than input 1
#
# Kyle Damm (01/02/2020)

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class FlowSolution():    
    def __init__(self):
        self.pos = np.array([])
        self.rho = np.array([])
        self.p = np.array([])
        self.v = np.array([])
            
    def read_flow_solution_from_file(self, filename):
        with open(filename, 'r') as f:
            data = f.readlines()
            for line in data:
                dat = line.split()
                self.pos = np.append(self.pos, float(dat[0]))
                self.rho = np.append(self.rho, float(dat[1]))
                self.p = np.append(self.p, float(dat[2]))
                self.v = np.append(self.v, float(dat[3]))

    def sample_from_other_flow_solution(self, sample_pos, other):
        for x in sample_pos:
            idx = binary_search(x, other.pos)
            err_msg = "WARNING: binary_search(x, nums) failed to find x={} in nums.".format(x)
            assert(idx != -1), err_msg
            self.pos = np.append(self.pos, x)
            self.rho = np.append(self.rho, other.rho[idx])
            self.p = np.append(self.p, other.p[idx])
            self.v = np.append(self.v, other.v[idx])

def plot_data(numeric_data, ref_data):

    err_msg = "WARNING: plot_data(numeric_data, ref_data) input data are not of equal length."
    assert (len(numeric_data.pos)==len(ref_data.pos)), err_msg
    x_loc = numeric_data.pos
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_size_inches(16, 4)
    
    ax1.plot(x_loc, numeric_data.p/1000, color='firebrick', label='Num.') # scale output from Pa to kPa
    ax1.plot(x_loc, ref_data.p/1000, color='royalblue',label='Ref.') # scale output from Pa to kPa
    ax1.set(xlabel='x (m)', ylabel='pressure ($kPa$)',
           title='Pressure distribution')
    ax1.grid()
    ax1.legend()
    
    ax2.plot(x_loc, numeric_data.rho, color='firebrick', label='Num.') 
    ax2.plot(x_loc, ref_data.rho, color='royalblue',label='Ref.') 
    ax2.set(xlabel='x (m)', ylabel='density  ($kg/m^3$)',
           title='Density distribution')
    ax2.grid()
    ax2.legend()

    ax3.plot(x_loc, numeric_data.v, color='firebrick', label='Num.') 
    ax3.plot(x_loc, ref_data.v, color='royalblue',label='Ref.') 
    ax3.set(xlabel='x (m)', ylabel='velocity  ($m/s$)',
           title='Velocity distribution')
    ax3.grid()
    ax3.legend()

    fig.savefig("flow_viz.png")
    

def relative_root_mean_square_error(vec1, vec2):
    err_msg = "WARNING: root_mean_square_error(vec1, vec2) input vectors are not of equal length."
    assert (len(vec1)==len(vec2)), err_msg
    n = len(vec1)
    return ( np.sqrt(np.sum(np.square(vec1-vec2))/n) ) / (np.mean(np.append(vec1,vec2)))

def linear_search(x, nums):
    # O(n)
    for i in range(len(nums)):
        if nums[i] == x:         # item found, return the index value
            return i
    return -1                    # loop finished, item was not in list

def binary_search(x, nums):
    # O(log n)
    low = 0
    high = len(nums) - 1
    while low <= high:           # There is still a range to search
        mid = (low + high)/2     # position of middle item
        item = nums[mid]
        if x == item :           # Found it! Return the index
            return mid
        elif x < item:           # x is in lower half of range
            high = mid - 1       #     move top marker down
        else:                    # x is in upper half of range
            low = mid + 1        #     move bottom marker up
    return -1                    # no range left to search,
                                 # x is not there

if __name__=='__main__':

    # process command line agruments
    assert (len(sys.argv) == 3), "{} arguments given, 2 expected.".format(len(sys.argv))
    numerical_filename = sys.argv[1]
    ref_filename = sys.argv[2]

    # solution from QUODAS
    numericSol = FlowSolution()
    numericSol.read_flow_solution_from_file(numerical_filename)

    # reference solution
    refSol = FlowSolution()
    refSol.read_flow_solution_from_file(ref_filename)

    # sample reference solution at numeric solution x positions
    sampledSol = FlowSolution()
    sampledSol.sample_from_other_flow_solution(numericSol.pos, refSol)

    # plot data
    plot_data(numericSol, sampledSol)
    
    # determine error
    err = [("density", relative_root_mean_square_error(sampledSol.rho, numericSol.rho)),
           ("pressure", relative_root_mean_square_error(sampledSol.p, numericSol.p)),
           ("velocity", relative_root_mean_square_error(sampledSol.v, numericSol.v))]
    passed = True
    for name,e in err:
        if (e > 0.1):
            passed = False
            print("FAILED on {} RMSE: {}".format(name, e));

    if (passed):
        print("=====================")
        print("==== TEST PASSED ====")
        print("=====================")
