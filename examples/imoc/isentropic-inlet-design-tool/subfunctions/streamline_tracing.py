# =============================================================================
# streamline_tracing.py traces the streamline from given start point until
# it reaches the desired turn angle.
# 
# Inputs:
#     start_indx,   int:  index of node at which the streamline is to start
#     comp_type,    str:  compression type; either 'ext' or 'int'
#     filename,     str:  filename (to save the data)
#     theta_turn, float:  turn angle in radians
#     dL,         float:  streamline spacing (step length)
#     dR,         float:  search radius (for more info: unit_process)
#     
# Outputs:
#     stream_indx, list:  list of streamline nodes indices
# 
# =============================================================================



import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit
import numpy as np



def streamline_tracing(start_indx, comp_type, filename,
                       theta_turn=0, dL=200e-6, dR=100.0):
    
    # streamline tracing
    stream_indx = [start_indx]
    kdt = kernel.create_kd_tree()
    
    # external contour
    if comp_type == 'ext':
        while kernel.nodes[stream_indx[-1]].theta < theta_turn:
            new_idx = unit.step_stream_node(stream_indx[-1], dL=dL, dR=dR, 
                                            kdtree=kdt)
            if new_idx is None: break
            stream_indx.append(new_idx)
    
    # internal contour
    elif comp_type == 'int':
        while kernel.nodes[stream_indx[-1]].theta > theta_turn:
            new_idx = unit.step_stream_node(stream_indx[-1], dL=dL, dR=dR, 
                                            kdtree=kdt)
            if new_idx is None: break
            stream_indx.append(new_idx)
            
    else:
        raise ValueError("Make sure to specify the compression type (comp_type)" + \
                   "as either external ('ext') or internal ('int')")
    
    
    
    kernel.register_streamline_start(stream_indx[0])

    # write external streamline points to a file
    stream = np.zeros([len(stream_indx), 2])
    for i in range(len(stream_indx)):
        stream[i][0] = kernel.nodes[stream_indx[i]].x
        stream[i][1] = kernel.nodes[stream_indx[i]].y
    np.savetxt(filename, stream, fmt = '%.8f')
    
    
    
    return stream_indx
    
    
    
