import numpy as np
import gzip as gz
import re
import os

###------------------------------------------------------------------------------------------------###
# Define some classes to assist in presenting the data
class Vector3:

    def __init__(self, x, y, z = 0):

        self.x = x; self.y = y; self.z = z

class flowState:

    def __init__(self, data_list):

        self.pos = Vector3(data_list[0], data_list[1], data_list[2])
        self.vol = data_list[3]
        self.rho = data_list[4]
        self.vel = Vector3(data_list[5], data_list[6], data_list[7])
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

###------------------------------------------------------------------------------------------------###

# Begin get_block_size

def get_block_size(jobName, block_id, tstep = 0, directory = "."):

    """
    get_block_size(jobName, block_id, tstep = 0, directory = "."):
    Gets the dimensions of the specified data block
    jobName: name of job as string
    block_id: id of block of interest as integer
    tstep: index of time file as integer 
    directory: location of the job
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

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")

    f = gz.open(filename, "rt")

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

#end get_block_size

###------------------------------------------------------------------------------------------------###

#begin get_cell_data

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

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")

    f = gz.open(filename, "rt")

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
        
    return flowState(cell_data)

#end get_cell_data

###------------------------------------------------------------------------------------------------###

#begin get_vertex

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

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step and block_id")
    
    f = gz.open(filename, "rt")

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
    line_number = 7 + (ni * nj) * k + ni * j + i

    # Pull out the information about the vertex
    f.seek(0)
    inc = 0
    vertex_vector = np.empty(3, dtype = float)

    for line in f:

        inc += 1

        if line_number == inc:
            for indx, element in enumerate(line.split()):
                vertex_vector[indx] = float(element)

    f.close()

    return Vector3(vertex_vector[0], vertex_vector[1], vertex_vector[2])

#end get_vertex

###------------------------------------------------------------------------------------------------###

#begin get_sim_time

def get_sim_time(jobName, tstep, block_id = 0, directory = "."):

    """
    get_sim_time(jobName, tstep, block_id = 0, directory = ".")
    jobName: name of job as string
    tstep: index of recorded flow solution to get time at
    block_id: index of block in case blocks running at different time
    directory: directory of job
    
    returns real time of simulation in seconds
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

    # Check that the desired file exists
    if not os.path.isfile(filename):
        raise IOError ("File does not exist- check job name, time step")

    f = gz.open(filename, "rt")

    # Pull out information about the data block
    inc = 0
    for line in f:
        
        inc += 1

        if inc == 3:
            simTime = float(line.split()[1])

    return simTime

#end get_sim_time

###------------------------------------------------------------------------------------------------###

#begin find_nearest_cell_centre

def find_nearest_cell_centre(jobName, nBlocks, x, y, z = 0.0, tstep = 0, directory = ".", quick = False):

    """
    find_nearest_cell_centre(jobName, x, y, z = 0.0, directory = ".", quick = True)
    jobName: your jobName as string
    nBlocks: number of blocks in the simulation
    x, y, z: location of interest
    tstep: index of flow solution, relevant for moving grids
    directory: directory of the job
    quick: either run quick version or slow version

    returns the block_id, i of the nearest cell to (x, y, z)

    The quick version will first attempt to find the block the cell is in by looking at the corners, then
    searching through the block to find the nearest cell. It does this by connecting the block corners with
    straight lines. This can fail if the point is near a block boundary that is not bounded by straight
    lines. The slow version simply checks every cell in the domain.
    """

    # Build the file name

    # Start from block_id of 0 and run until the file no longer exists
    block_id = 0
    finished = False
    dist = 1e6

    while not finished:

        if quick:

            # The quick version- first perform a bounding box check to see if the point is in the current block

            ni, nj, nk = get_block_size(jobName, block_id)

            block_corners = [get_vertex(jobName, block_id, 0, 0), get_vertex(jobName, block_id, 0, nj), get_vertex(jobName, block_id, ni, nj), get_vertex(jobName, block_id, ni, 0)]

            if is_in_bounding_box(Vector3(x, y), block_corners):

                # Open the desired file

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

                # Check that the desired file exists
                if not os.path.isfile(filename):
                    raise IOError ("File does not exist- check job name, time step")

                f = gz.open(filename, 'rt')

                # Run through the file and check distances at each point

                line_number = 0
                for line in f:

                    line_number += 1

                    if line_number >= 10:

                        line_data = line.split()

                        local_dist = np.sqrt((x - float(line_data[0])) ** 2 + (y - float(line_data[1])) ** 2 + (z - float(line_data[2])) ** 2) 

                        if local_dist < dist:
                            dist = local_dist
                            cell_id = line_number - 10
                
                # Return the local block and tell the function we're done
                local_block = block_id
                finished = True
                f.close()

        else:

            # The flow version which simply checks every cell in the domain
            # Open the desired file

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

            # Check that the desired file exists
            if not os.path.isfile(filename):
                raise IOError ("File does not exist- check job name, time step")

            f = gz.open(filename, 'rt')

            # Run through the file and check distances at each point

            line_number = 0
            for line in f:

                line_number += 1

                if line_number >= 10:

                    line_data = line.split()

                    local_dist = np.sqrt((x - float(line_data[0])) ** 2 + (y - float(line_data[1])) ** 2 + (z - float(line_data[2])) ** 2) 

                    if local_dist < dist:
                        dist = local_dist
                        cell_id = line_number - 10
                        local_block = block_id

            f.close()

        block_id += 1

        if block_id == nBlocks:
            finished = True

    return local_block, cell_id

# end find_nearest_cell_centre

###------------------------------------------------------------------------------------------------###

# Begin is_in_bounding_box_function
def is_in_bounding_box(pos, corners):
    
    """
    is_in_bounding_boX(pos, corners)
    pos: Vector3 of the current position
    corners: list of 4 Vector3 objects

    returns boolean stating whether the stated position is in the bounding box defined by the corners

    Order of the corners should be either clockwise or counter-clockwise, precise order does not matter.
    """

    # See c++ code code by Dan Sunday for this winding number implementation

    # Define the isLeft lambda function
    isLeft = lambda P, V0, V1 : ((V1.x - V0.x) * (P.y - V0.y) - (P.x - V0.x) * (V1.y - V0.y))

    # Now perform the winding number algorithm
    wn = 0

    for i in range(len(corners)):
        # Prevent stepping off the end of the list
        if i == 3:
            j = 0
        else:
            j = i + 1

        if corners[i].y <= pos.y:
            if corners[j].y > pos.y:
                # If our point is between the two corners in the y, must have upward crossing
                if isLeft(pos, corners[i], corners[j]) > 0:
                    # Pos is to left
                    wn += 1

        else:
            if corners[j].y <= pos.y:
                # Downward crossing
                if isLeft(pos, corners[i], corners[j]) < 0:
                    wn -= 1

    if wn == 0:
        return False
    else:
        return True

# end is_in_bounding_box

###------------------------------------------------------------------------------------------------###







if __name__ == "__main__":

    jobname = "IVA"; tstep = 2; block_id = 2;
    i = 5; j = 2;

    print(get_cell(jobname, tstep, block_id, i, j))
