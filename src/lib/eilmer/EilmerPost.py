import numpy as np
import gzip as gz
import re
import os

class flowState:

    def __init__(self, data_list):

        self.posx = data_list[0]
        self.posy = data_list[1]
        self.posz = data_list[2]
        self.vol = data_list[3]
        self.rho = data_list[4]
        self.velx = data_list[5]
        self.vely = data_list[6]
        self.velz = data_list[7]
        self.p = data_list[8]
        self.a = data_list[9]
        self.mu = data_list[10]
        self.k = data_list[11]
        self.mu_t = data_list[12]
        self.k_t = data_list[13]
        self.S = data_list[14]
        self.massf = data_list[15]
        self.u = data_list[16]
        self.T = data_list[17]

def get_block_size(jobName, tstep, block_id, directory = "."):

    """
    get_block_size(jobName, tstep, block_id):
    Gets the dimensions of the specified data block
    jobName: name of job as string
    tstep: index of time file as integer
    block_id: id of block of interest as integer
    returns ni, nj, nk (dimensions of the data block in terms of number of cells)
    """

    # Build the file name
    if 100 <= block_id < 1000:
        block_str = "0" + str(block_id)
    elif 10 <= block_id < 100:
        block_str = "00" + str(block_id)
    elif 0 <= block_id < 10:
        block_str = "000" + str(block_id)
    else:
        raise ValueError ("Block ID out of range")

    if 100 <= tstep < 1000:
        time_str = "0" + str(tstep)
    elif 10 <= tstep < 100:
        time_str = "00" + str(tstep)
    elif 0 <= tstep < 10:
        time_str = "000" + str(tstep)
    else:
        raise ValueError ("Time ID out of range")
    
    filename = directory + "/flow/t" + time_str + "/" + jobName + ".flow.b" + block_str + ".t" + time_str + ".gz"
    f = gz.open(filename, "rt")

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")

    # Pull out information about the data block
    inc = 0
    for line in f:
        
        inc += 1

        if inc == 3:
            simTime = float(line.split()[1])
        if inc == 7:
            ni = int(line.split()[1])
        if inc == 8:
            nj = int(line.split()[1])
        if inc == 9:
            nk = int(line.split()[1])

    # No useful information past here- stop searching
        if inc == 10:
            break

    f.close()

    return ni, nj, nk

def get_cell_data(jobName, tstep, block_id, i, j = 0, k = 0, directory = "."):

    """
    get_cell_data(jobName, tstep, block_id, i, j = 0, k = 0):
    Get the full set of cell data from a given cell
    jobName: your jobname as a string
    tstep: index of time file as integer
    block_id: data block containing cell of interest as integer
    i, j, k: index of cell in the block

    Returns: [x, y, z, volume, density, v_x, v_y, v_z, p, a, mu, k, mu_t, k_t, S, massf, u (int energy), T], simTime
    """

    # Build the file name
    if 100 <= block_id < 1000:
        block_str = "0" + str(block_id)
    elif 10 <= block_id < 100:
        block_str = "00" + str(block_id)
    elif 0 <= block_id < 10:
        block_str = "000" + str(block_id)
    else:
        raise ValueError ("Block ID out of range")

    if 100 <= tstep < 1000:
        time_str = "0" + str(tstep)
    elif 10 <= tstep < 100:
        time_str = "00" + str(tstep)
    elif 0 <= tstep < 10:
        time_str = "000" + str(tstep)
    else:
        raise ValueError ("Time ID out of range")
    
    filename = directory + "/flow/t" + time_str + "/" + jobName + ".flow.b" + block_str + ".t" + time_str + ".gz"
    f = gz.open(filename, "rt")

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")

    # Pull out information about the data block
    inc = 0
    for line in f:
        
        inc += 1

        if inc == 3:
            simTime = float(line.split()[1])
        if inc == 7:
            ni = int(line.split()[1])
        if inc == 8:
            nj = int(line.split()[1])
        if inc == 9:
            nk = int(line.split()[1])
            break

    # Check the inputs
    if i >= (ni * nj * nk):
        raise ValueError ("i index out of range")
    if j >= nj:
        raise ValueError ("j index out of range")
    if k >= nk:
        raise ValueError ("k index out of range")

    # Now we know the line number we're looking for
    line_number = 10 + (ni + nj) * k + ni * j + i

    # Pull out the data from the selected line- order of variables is: x, y, z, volume, rho, vx, vy, vz, p, a, mu, k, mu_t, k_t, S, massf, u, T
    f.seek(0)
    inc = 0
    cell_data = np.empty(18, dtype = float)

    for line in f:

        inc += 1

        if inc == line_number:
            for indx, element in enumerate(line.split()):
                cell_data[indx] = float(element)

    f.close()
        
    return flowState(cell_data), simTime

def get_vertex(jobName, block_id, i, j = 0, k = 0, tstep = 0, directory = "."):

    """
    get_vertex(jobName, block_id, i, j = 0, k = 0, tstep = 0, directory = ".")
    Gets the coordinates of a given vertex in the block
    jobName: your jobName as a string
    block_id: local block id to search in
    i, j = 0, k = 0: indices of the vertex
    tstep = 0: index of the flow solution, only relevant for moving grid simulations
    directory = ".": directory of the job files

    Returns x, y, z: coordinates of the vertex
    """

    # Build the file name
    if 100 <= block_id < 1000:
        block_str = "0" + str(block_id)
    elif 10 <= block_id < 100:
        block_str = "00" + str(block_id)
    elif 0 <= block_id < 10:
        block_str = "000" + str(block_id)
    else:
        raise ValueError ("Block ID out of range")

    if 100 <= tstep < 1000:
        time_str = "0" + str(tstep)
    elif 10 <= tstep < 100:
        time_str = "00" + str(tstep)
    elif 0 <= tstep < 10:
        time_str = "000" + str(tstep)
    else:
        raise ValueError ("Time ID out of range")
    
    filename = directory + "/grid/t" + time_str + "/" + jobName + ".grid.b" + block_str + ".t" + time_str + ".gz"
    f = gz.open(filename, "rt")

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")

    # Pull out information about the block
    inc = 0
    for line in f:

        inc += 1

        if inc == 4:
            ni = int(line.split()[1])
        if inc == 5:
            nj = int(line.split()[1])
        if inc == 6:
            nk = int(line.split()[1])
            break

    # Check the inputs
    if i >= (ni * nj * nk):
        raise ValueError ("i index out of range")
    if j >= nj:
        raise ValueError ("j index out of range")
    if k >= nk:
        raise ValueError ("k index out of range")

    # What line number are we looking for?
    line_number = 7 + (ni * nj) * nk + ni * j + i

    # Pull out the information about the vertex
    f.seek(0)
    inc = 0
    vector = np.empty(3, dtype = float)

    for line in f:

        inc += 1

        if line_number == inc:
            for indx, element in enumerate(line.split()):
                vector[indx] = float(element)

    f.close()

    return vector





if __name__ == "__main__":

    jobname = "IVA"; tstep = 2; block_id = 2;
    i = 5; j = 2;

    print(get_cell(jobname, tstep, block_id, i, j))
