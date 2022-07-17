import numpy as np
import matplotlib.pyplot as plt



def plotNodes(nodes):
    info = {
        "indx_array": np.zeros(len(nodes)),
        "y": np.zeros(len(nodes)),
        "x":  np.zeros(len(nodes))
    }

    for i in range(len(nodes)):
        # print(type(kernel.nodes[i].indx))
        # if kernel.nodes[i].x < xs[-1]+0.01:
            info["indx_array"][i] = int(nodes[i].indx)
            info["y"][i] = nodes[i].y
            info["x"][i] = nodes[i].x

    plt.figure(0, figsize=(60,12))
    #plt.figure(0, figsize=(20,4))
    plt.plot(info["x"], info["y"], markersize=1, linestyle="None", marker='o'
    , color = 'k')
# =============================================================================
#     plt.xlim([0.17, 0.22])
#     plt.ylim([0.03, 0.06])    
#     plt.grid(which='both')
# =============================================================================
    plt.show()