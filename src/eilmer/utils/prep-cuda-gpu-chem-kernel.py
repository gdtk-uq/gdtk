#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
import sys
import os
#------------------------------------------------------------------------------------
# Python module for generating system of ODEs for finite-rate chemical kinetics
#  - adapted from compact notation method
# author: Kyle Damm, 2015
#------------------------------------------------------------------------------------
'''
Start by finding and storing the correct path to the kernel template
'''

src_dir = os.getenv("DGD")
tmplt_file = src_dir + "/share/alpha_qss_kernel_cuda_tmplt.cu"
#tmplt_file = '/home/kdamm/codes/eilmer4/dgd-gpu-dev/src/eilmer/utils/alpha_qss_kernel_cuda_tmplt.cu'
#------------------------------------------------------------------------------------
# 			      Read in Data 
#------------------------------------------------------------------------------------
''' 
Let's read in the number of reactions, number of species, the stoichiometric 
coefficients and the third body efficiencies from the input file. The input file
is located in the source directory along with the template file.
'''
print "reading in reaction data file"

vji_p = []
MM = []

with open('chem-compact-notation.inp', 'r') as f:
    data = f.readlines()
    count = 0
    for line in data:
        dat = line.split()
        if count==0:                   #reading in no. species/reactions
            num_spec = int(dat[0])
            num_reac = int(dat[1])
            vji_pp = zeros((num_reac, num_spec))
            vji_p = zeros((num_reac, num_spec))
            MM = zeros((num_reac, num_spec+1))
        if count>1 and count<num_reac+2:                  #stoichiometric coeffs (reactants)
            for i in linspace(0, len(dat)-1, len(dat)):
                i = int(i)
                vji_p[count-2][i] = dat[i]
        if count>num_reac+2 and count<num_reac*2+3:       #stoichiometric coeffs (products)
            for i in linspace(0, len(dat)-1, len(dat)):
                i = int(i)
                vji_pp[count-num_reac-3][i] = dat[i]
        if count>num_reac*2+3 and count<num_reac*3+4:    #third body efficiencies
            if int(dat[0])==0:
                MM[count-num_reac*2-4][0]=dat[0]
            if int(dat[0])==1:
                MM[count-num_reac*2-4][0]=dat[0]
                for i in linspace(0, len(dat)-1, len(dat)):
                    i = int(i)
                    MM[count-num_reac*2-4][i] = dat[i]
            elif int(dat[0])>1:
                print "Error: encountered integer above 1"
                sys.exit()

        count += 1

#print vji_pp
#print vji_p
#print MM
print "data read successfully"

#------------------------------------------------------------------------------------
# 					M equations 
#------------------------------------------------------------------------------------
'''
Now let's form the M equations for the third body colliders. There theoretically
could be an M equation for every reaction. This will generate an equation only
for reactions with M in them, and only multiples species concentrations with
non-zero efficiencies to save on floating point operations.
'''

print "setting up M equation"

with open('workfile', 'w') as f:
    for i in linspace(0, num_reac-1, num_reac):
        if MM[i][0] == 0:
            pass
        else:
            i = int(i)
            f.write('MM[idx+numcell*')
            f.write(str(int(i)))
            f.write(']=')
            for j in linspace(0, num_spec-1, num_spec):
                if MM[i][1+j] == 0.0:
                    pass
                else:
                    value = 'Y[idx+numcell*'
                    s = str(value)
                    f.write(s)
                    f.write(str(int(j)))
                    f.write(']*')
                    f.write(str(MM[i][1+j]))
                    if j < num_spec - 1:
                        if j == num_spec - 2 and MM[i][2+j] == 0.0:  #watch out that you don't have a '+' if last efficiency is 0
                            pass
                        else:
                            f.write('+')
                if j >= num_spec - 1:
                    f.write(';')
                    f.write('\n')

#------------------------------------------------------------------------------------
# 					ODEs - Production terms (q)
#------------------------------------------------------------------------------------
    '''
        We will now formulate the set of ODEs, seperated into production terms q and loss terms,
        p, which fit the ODE structure of dy/dt = q - py. This is based off of the compact notation
        method found in Turns' textbook An Introduction to Combustion.
        '''
    print "setting up Production equations (q)"
    
    vji = zeros((num_reac, num_spec))  
    for j in linspace(0, num_spec-1, num_spec):
        for i in linspace(0, num_reac-1, num_reac):
            i = int(i)
            j = int(j)            
            vji[i][j] = vji_pp[i][j] - vji_p[i][j]
    for j in linspace(0, num_spec-1, num_spec):
        zero_count = 0
        value = 'q[idx+numcell*'
        s = str(value)
        f.write(s)
        f.write(str(int(j)))
        f.write(']=')
        for i in linspace(0, num_reac-1, num_reac):
            if vji[i][j] == 0:
                zero_count += 1
                pass
            else:
                if vji[i][j] < 0:			# if this is -ve then the production is from kb
                    value = '+kb[idx+numcell*'
                    s = str(value)						
                    f.write(s)
                    f.write(str(int(i)))
                    f.write(']')
                    f.write('*')
                    f.write(str(abs(vji[i][j])))
                    for l in linspace(0, num_spec-1, num_spec):
                        l = int(l)
                        i = int(i)
                        if vji_pp[i][l] > 0:
                            for pow in linspace(0, vji_pp[i][l]-1, vji_pp[i][l]):
                                value = '*Y[idx+numcell*'	
                                s = str(value)
                                f.write(s)
                                f.write(str(int(l)))
                                f.write(']')
                    if MM[i][0] > 0.0:
                        value = '*MM[idx+numcell*'	
                        s = str(value)
                        f.write(s)
                        f.write(str(int(i)))
                        f.write(']')
                if vji[i][j] > 0:			# if this is +ve then the production is from kf
                    value = '+kf[idx+numcell*'
                    s = str(value)						
                    f.write(s)
                    f.write(str(int(i)))
                    f.write(']')
                    f.write('*')
                    f.write(str(abs(vji[i][j])))
                    for l in linspace(0, num_spec-1, num_spec):
                        l = int(l)
                        i = int(i)
                        if vji_p[i][l] > 0:
                            for pow in linspace(0, vji_p[i][l]-1, vji_p[i][l]):
                                value = '*Y[idx+numcell*'	
                                s = str(value)
                                f.write(s)
                                f.write(str(int(l)))
                                f.write(']')
                    if MM[i][0] > 0.0:
                        value = '*MM[idx+numcell*'	
                        s = str(value)
                        f.write(s)
                        f.write(str(int(i)))
                        f.write(']')
        if (zero_count == num_reac):
            f.write('0.0')
        f.write(';')
        f.write('\n')
#------------------------------------------------------------------------------------
# 					ODEs - Loss terms (p)
#------------------------------------------------------------------------------------
    print "setting up Loss equation (p)"
    
    vji = zeros((num_reac, num_spec))  
    for j in linspace(0, num_spec-1, num_spec):
        for i in linspace(0, num_reac-1, num_reac):
            i = int(i)
            j = int(j)            
            vji[i][j] = vji_pp[i][j] - vji_p[i][j]
    for j in linspace(0, num_spec-1, num_spec):
        zero_count = 0
        value = 'p[idx+numcell*'
        s = str(value)
        f.write(s)
        f.write(str(int(j)))
        f.write(']=')
        for i in linspace(0, num_reac-1, num_reac):
            if vji[i][j] == 0:
                zero_count += 1
                pass
            else:
                if vji[i][j] < 0:			# if this is -ve then the loss is from kf
                    value = '+kf[idx+numcell*'
                    s = str(value)						
                    f.write(s)
                    f.write(str(int(i)))
                    f.write(']')
                    f.write('*')
                    f.write(str(abs(vji[i][j])))
                    for l in linspace(0, num_spec-1, num_spec):
                        l = int(l)
                        i = int(i)
                        if vji_p[i][l] > 0:
                            if j == l:				#check for squared terms to make sure only multiply once (p*y in mott algorithm)
                                if vji_p[i][l] == 1.0:
                                    pass
                                elif vji_p[i][l] > 1.0:
                                    for pow in linspace(0, vji_p[i][l]-2, vji_p[i][l]-1):
                                        value = '*Y[idx+numcell*'	
                                        s = str(value)
                                        f.write(s)
                                        f.write(str(int(l)))
                                        f.write(']')
                            else:
                                for pow in linspace(0, vji_p[i][l]-1, vji_p[i][l]):
                                    value = '*Y[idx+numcell*'	
                                    s = str(value)
                                    f.write(s)
                                    f.write(str(int(l)))
                                    f.write(']')
                    if MM[i][0] > 0.0:
                        value = '*MM[idx+numcell*'	
                        s = str(value)
                        f.write(s)
                        f.write(str(int(i)))
                        f.write(']')
                if vji[i][j] > 0:			# if this is +ve then the loss is from kb
                    value = '+kb[idx+numcell*'
                    s = str(value)						
                    f.write(s)
                    f.write(str(int(i)))
                    f.write(']')
                    f.write('*')
                    f.write(str(abs(vji[i][j])))
                    for l in linspace(0, num_spec-1, num_spec):
                        l = int(l)
                        i = int(i)
                        if vji_pp[i][l] > 0:
                            if j == l:			#check for squared terms to make sure only multiply once (p*y in mott algorithm)
                                if vji_pp[i][l] == 1.0:
                                    pass
                                elif vji_pp[i][l] > 1.0:
                                    for pow in linspace(0, vji_pp[i][l]-2, vji_pp[i][l]-1):
                                        value = '*Y[idx+numcell*'	
                                        s = str(value)
                                        f.write(s)
                                        f.write(str(int(l)))
                                        f.write(']')
                            else:	
                                for pow in linspace(0, vji_pp[i][l]-1, vji_pp[i][l]):
                                    value = '*Y[idx+numcell*'	
                                    s = str(value)
                                    f.write(s)
                                    f.write(str(int(l)))
                                    f.write(']')
                    if MM[i][0] > 0.0:
                        value = '*MM[idx+numcell*'	
                        s = str(value)
                        f.write(s)
                        f.write(str(int(i)))
                        f.write(']')
        if (zero_count == num_reac):
            f.write('0.0')
        f.write(';')
        f.write('\n')
f.close
#------------------------------------------------------------------------------------
# 			        Piece Together Kernel
#------------------------------------------------------------------------------------
'''
Let's now piece together the equations with the kernel template file. 
'''
print "putting Kernel together"

in_file = open('workfile')   #open and read set of equations
indata = in_file.read()    
fin = open(tmplt_file, 'r')  #open template and store contents
template_text = fin.read()
fin.close()          #close template file
update_text = template_text.replace('// <insert-source-terms-here>',indata) #place in equations
fout = open('alpha_qss_kernel.cu', 'w') #open and store kernel for gpu
fout.write(update_text)
fout.close()           #close update file
os.remove("workfile")  #delete temporary equations file

# now build cuda kernel into shared object
cmd = "nvcc --shared -o libcudakernel.so alpha_qss_kernel.cu --compiler-options '-fPIC'"
os.system(cmd)
dest = src_dir + "/lib/"
cmd = "cp libcudakernel.so " + dest
os.system(cmd)

print "Kernel ready for launch."
#------------------------------------------------------------------------------------
# 			       Kernel is ready
#------------------------------------------------------------------------------------
