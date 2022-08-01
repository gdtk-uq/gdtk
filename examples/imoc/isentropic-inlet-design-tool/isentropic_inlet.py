# =============================================================================
# isentropic_inlet.py is a tool for inviscid design of an isentropic inlet 
# with a specified cone tip semi-angle. Starting with a conical flowfield  
# (from solution of Taylor-Maccoll equations), it uses Eilmer-embedded method
# of characteristics (IMOC) to model the flowfield. Then the inlet contour 
# is drawn using streamline tracing. Then contour is straightened back
# to horizontal with an arc. Finally, the cowl is drawn using two arcs with
# a specified angle in between them. The output of isentropic_inlet.py are
# 4 contours in form of .dat files.
# 
# Maciej Grybko (maciej.grybko@usq.edu.au)  
# University of Southern Queensland
# =============================================================================





from gdtk.ideal_gas_flow import PM1, PM2, beta_cone2, theta_cone_flowfield
from subfunctions.plot_nodes import plotNodes
from subfunctions.streamline_tracing import streamline_tracing
import numpy as np
import matplotlib.pyplot as plt
import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit





# =============================================================================
# INITIATION AND SETTINGS - feel free to adjust them
# =============================================================================

# freestream properties
M1 = 4.
p1 = 1.
T1 = 1.
R = 287.1
g = 1.4
a1 = np.sqrt(g*R*T1)
V1 = M1*a1

# conical flowfield settings
rays_num = 500              # number of rays emanating from the cone tip,
                            # i.e. characteristic mesh resolution (recom. 5000)
theta_cone = np.deg2rad(10) # semi-angle of the cone (inlet tip)
                              # converted from degrees to radians
dtheta = -1e-5              # angle increment for Taylor Maccoll solution
                              # recommended values:
                              #   -1e-5 for quick solution (2 min @ 5000 rays)
                              #   -1e-7 for refined solution (20 min)
tol = -dtheta/1000          # bisection tolerance

beta = beta_cone2(M1, theta_cone, R, g, tol, dtheta) # conical shock angle

# conical flowfield properties
M, flow_dir, ray_angles, mu = \
    theta_cone_flowfield(V1, p1, T1, beta, theta_cone, 
                         rays_num, R, g, dtheta)

# IMOC parameters
kernel.axisymmetric = True
kernel.g = 1.4

# geometric parameters
y_focal = 0.050             # focal point (the cowl tip)
x_focal = y_focal/np.tan(beta)
kernel.Node(x = x_focal, y = y_focal, nu = PM1(M[0]), mach = M[0], 
            theta = flow_dir[0]) # initialise the focal node
theta_turn = np.deg2rad(30) # turn angle (the degree of compression)

# characteristic mesh parameters
N_theta = int(rays_num/2*theta_turn/np.deg2rad(30)) # number of the flow 
                                                    # turn angle increments
corr_factor = 0.5 # user-defined factor to ensure proximity of the start-line
                    # to the characteristic mesh (must be in range (0, 1])
trim = 0.9        # trim coefficient. It is necessary to trim off the bottom
                    # of the characteristic mesh to prevent numerical errors
                    # when approaching the inlet axis (r=0)
                    # (must be in range (0, 1])
                    
# streamline spacing
dL_stream = 0.04/rays_num

# To ensure reasonably uniform spacing of the nodes, non-uniform
# theta increment is required. The increment is calculated based on
# geometric series. I.e. each next term is multiplied by a constant, q.
### q IS ESSENTIAL IN UNIFORMITY OF THE MESH AND MAY REQUIRE TWEAKING ###
S = theta_turn - kernel.nodes[0].theta*corr_factor
q = 1 + 35/np.rad2deg(theta_cone)/rays_num # should be in range (1, 1.05)

# turning radii scale factors (based on distance between focal point and end
# of the isentropic contour) for the lower and upper arcs respectively
r_scale = 4.2
r_scale2 = r_scale - 0.5

# cowl wedge angle
cowl_angle = np.deg2rad(10)

# =============================================================================
# 
# =============================================================================





# =============================================================================
# ISENTROPIC CONTOUR
# =============================================================================

# tracing the first mach wave from the focal point down to the cone surface
tol = 1e-9
for i in range(rays_num-1):
    
    # initially create a node on desired ray (i.e. desired theta)
    slope = (mu[i]+flow_dir[i] + mu[i+1]+flow_dir[i+1])/2
    new_x = (kernel.nodes[-1].y - kernel.nodes[-1].x*np.tan(slope)) \
        / (np.tan(ray_angles[i+1]) - np.tan(slope))
    new_y = new_x*np.tan(ray_angles[i+1])
    kernel.Node(x = new_x, y = new_y, 
                nu = PM1(M[i+1]), mach = M[i+1], theta = flow_dir[i+1])
    
    # keep correcting node position until it lies on characteristic line
    # (in other words, until it concides with one generated with unit.interior)
    dL = 1                              # pre-define position error
    while dL > tol:
        
        # compare position with IMOC generated node
        if i==0:
            unit.interior(i+1, 0, -1)
        else:
            unit.interior(i+1, i, -1)
        dx = kernel.nodes[-1].x - kernel.nodes[-2].x
        dy = kernel.nodes[-1].y - kernel.nodes[-2].y
        dL = np.sqrt(dx**2 + dy**2)
        
        # correct node position
        slope = (y_focal - kernel.nodes[-1].y) / (x_focal - kernel.nodes[-1].x)
        new_x = (slope*x_focal - y_focal) / (slope - np.tan(ray_angles[i+1]))
        new_y = new_x*np.tan(ray_angles[i+1])
        kernel.nodes[-2].x = new_x
        kernel.nodes[-2].y = new_y
        
        # delete IMOC generated node
        kernel.nodes.pop()



# characteristic mesh
theta_focal = kernel.nodes[0].theta*corr_factor # initial flow angle at cowl
i = len(kernel.nodes) - 2 # pre-define number of points in characteristic line
for j in range(N_theta):
    
    # calculate theta increment (geometric series)
    delta_theta = S*(1-q)/(1-q**N_theta)*q**j
        
    # calculate corresponding Mach num
    theta_focal += delta_theta
    M2 = PM2(PM1(M[0]) - theta_focal)   
    
    # overwrite the focal node
    kernel.Node(x = x_focal, y = y_focal, nu = PM1(M2),
                mach = M2, theta = theta_focal)
    
    # draw the new characteristic line
    k = i + 1
    for i in range(k - 1):
        unit.interior(len(kernel.nodes) - k - 1, len(kernel.nodes) - 1, -1)
                
        # prevent creating invalid nodes (due to numerical errors
        # when approaching the inlet axis)
        if kernel.nodes[-1].y < kernel.nodes[-1].x * np.tan(trim*theta_cone): 
            break

    

# external streamline tracing
ext_stream_indx = rays_num-1
ext_stream_indx = streamline_tracing(ext_stream_indx, 'ext', 
                                     'External_streamline_data.dat',
                                     theta_turn, dL=dL_stream)

# last streamline node
last_stream = kernel.nodes[-1].indx
x_last_stream = kernel.nodes[last_stream].x
y_last_stream = kernel.nodes[last_stream].y

# find the nearest characteristic mesh node to the last streamline node
nearest_node = last_stream
x_nearest = x_last_stream
y_nearest = y_last_stream
while nearest_node == last_stream:
    nearest_node = kernel.find_nodes_near(x_nearest, y_nearest)[0]
    x_nearest += 1e-6
    y_nearest += 1e-6

# find the last focal node
last_focal_indx = nearest_node
while kernel.nodes[last_focal_indx].cplus_down != None:
    last_focal_indx = kernel.nodes[last_focal_indx].cplus_down
    
# find cowl flow direction
cowl_flow_dir = kernel.nodes[last_focal_indx].theta

# =============================================================================
# 
# =============================================================================





# =============================================================================
# INTERNAL ARC
# =============================================================================

# find inlet opening height
h_open = np.sqrt((x_focal - x_last_stream)**2 + (y_focal - y_last_stream)**2)

# find centre of the arc
x_cen1 = x_last_stream + r_scale*h_open*np.sin(theta_turn)
y_cen1 = y_last_stream - r_scale*h_open*np.cos(theta_turn)

# define the arc
arc_num = rays_num # at least 2000 for refined solution
arc = np.zeros([arc_num, 2])
for i in range(arc_num):
    theta = (1 - i/arc_num) * theta_turn
    arc[i][0] = x_cen1 - r_scale*h_open*np.sin(theta)
    arc[i][1] = y_cen1 + r_scale*h_open*np.cos(theta)
    
# save the arc coordinates
np.savetxt("Internal_arc.dat", arc, fmt = '%.8f')

# =============================================================================
# 
# =============================================================================





# =============================================================================
# INTERNAL COWL ARC
# =============================================================================

# find centre of the cowl arc
x_cen2 = x_focal + r_scale2*h_open*np.sin(cowl_flow_dir)
y_cen2 = y_focal - r_scale2*h_open*np.cos(cowl_flow_dir)

# define the cowl arc
int_cowl_arc = np.zeros([arc_num, 2])
for i in range(arc_num):
    theta = (1 - i/arc_num) * cowl_flow_dir
    int_cowl_arc[i][0] = x_cen2 - r_scale2*h_open*np.sin(theta)
    int_cowl_arc[i][1] = y_cen2 + r_scale2*h_open*np.cos(theta)

# save the arc coordinates
np.savetxt("Internal_cowl_arc.dat", int_cowl_arc, fmt = '%.8f')

# =============================================================================
# 
# =============================================================================





# =============================================================================
# EXTERNAL COWL ARC
# =============================================================================

# define the external cowl flow direction 
ext_cowl_flow_dir = cowl_flow_dir + cowl_angle

# find the external cowl arc centre and radius , so that the centre is 
# vertically above the internal cowl arc centre
x_cen3 = x_cen2
y_cen3 = y_focal - (x_cen3 - x_focal)/np.tan(ext_cowl_flow_dir)
r = np.sqrt((x_cen3 - x_focal)**2 + (y_cen3 - y_focal)**2)

# define the cowl arc
ext_cowl_arc = np.zeros([arc_num, 2])
for i in range(arc_num):
    theta = (1 - i/arc_num) * ext_cowl_flow_dir
    ext_cowl_arc[i][0] = x_cen3 - r*np.sin(theta)
    ext_cowl_arc[i][1] = y_cen3 + r*np.cos(theta)

# save the arc coordinates
np.savetxt("External_cowl_arc.dat", ext_cowl_arc, fmt = '%.8f')

# define isolator height
h_iso = int_cowl_arc[-1][1] - arc[-1][1]

# calculate radius to isolator height ratio
R2h_iso_ratio = r_scale*h_open/h_iso

# =============================================================================
# 
# =============================================================================





# =============================================================================
# PLOTTING
# =============================================================================

plotNodes(kernel.nodes)

# plot complete geometry
plt.figure(figsize=(20,4))
plt.plot([kernel.nodes[i].x for i in \
          ext_stream_indx],
         [kernel.nodes[i].y for i in \
          ext_stream_indx])
plt.plot(arc[:, [0]], arc[:, [1]])
plt.plot(int_cowl_arc[:, [0]], int_cowl_arc[:, [1]])
plt.plot(ext_cowl_arc[:, [0]], ext_cowl_arc[:, [1]])
plt.title('Inlet geometry')
plt.legend(['External_streamline_data.dat', 
            'Internal_arc.dat',
            'Internal_cowl_arc.dat',
            'External_cowl_arc.dat'])

# flow properties along the streamline
plt.figure(figsize=(200,40))
plt.plot([kernel.nodes[i].mach for i in ext_stream_indx])
plt.title('Mach number along the isentropic contour')
plt.figure(figsize=(200,40))
plt.plot([np.rad2deg(kernel.nodes[i].theta) for i in ext_stream_indx])
plt.title('Flow direction along the isentropic contour (deg)')










