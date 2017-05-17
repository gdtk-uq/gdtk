Test Case for Sliding Mesh simualtion in Eilmer 4

This example presents the implementation of a sliding mesh itnerface between two 
segments of rotating mesh. From a top view the fluid domain is described by 
three concentric circles defining the inlet, sliding interface, and outlet
respectively. For the simulation the segments are defined with an angle of 
+/-theta degrees. 

To allow testing of more complex rotating mesh interfaces, consisting of multiple
blocks on both sides the outer section of the mesh is split in two blocks in the 
tangential direction and the inner section of the mesh is split in two blocks in 
the axial direction. 

During the simulatuions the outer part is kept stationary and the inner part is
rotated with angular speed OMEGA as far as the interface is concerned 
(Rotational effects acting on the inner blocks can be added by setting the 
z_rotation flag during the block definition). Using the sliding interface the 
area weigthed fluxes of the conserved quantities are exchanged between both 
sides of the interfaces. The areas of overalp between opposing cell faces are 
calculated directly at every timestep based on the rotation angle A = OMEAG * t. 

The accompanying test case Sliding_Mesh_Test.lua and the associated lua scripts 
can be executed using the scripts
$ . prep_run.sh     # to create the mesh and start simualtion
$ . post.sh         # to post-process the simualtion

The current restriction of the sliding mesh interface are:
- upstream side of the interface must be defined by WEST faces
- downstream side of the interface must be defined by SOUTH faces
- both sides of the interface must have a structured mesh
- currently only angle segments of < 180deg have been trialled
- currently only implemented for inviscid simulations. 


The simulation set up is created in two files.
job.lua --> defined the mesh and define sides of the sliding mesh interface. See
        Sliding_Mesh_test.py as an example.
udf-config_e4.lua --> set up file that defines the required data for the user 
        define boundary condition. See inside the file for detailed instuctions. 

The minimum parameters that need to be defined are:
N_BLADE = 16        -- Number of blades
OMEGA = 20000       -- rotational speed in rad/s
theta_min = -xx     -- angle at which interface starts in rad
theta_min = +xx     -- angle at which interface ends in rad
relax_SR = 0.5      -- relaxtion factor applied to forwards flow
relax_RS = 0.5      -- relaxtion factor applied to backwards flow
upstream_row0 = {0,1} --blk0, blk1
UP_row_list = {upstream_row0}       -- list of rows of blocks defining upstream 
                                    -- side of interface
downstream_row0 = {2}  -- blk2
downstream_row1 = {3}  -- blk3
DOWN_row_list = {downstream_row0, downstream_row1}  -- list of rows of blocks 
                                                    -- defining downstream side 
                                                    -- of interface

The output from this process is the data file "CONFIG_tbl.lua" which stores the 
data required for the set up of the sliding interface. 

The actual set-up of the BC is performed by udf-process_e3.lua. This loads the 
input data created by the lua configuration script, performs some cell 
interpolation tasks and stores the requried data in the files "UP_tbl.lua" and 
"DOWN_tbl.lua". 

Finally the boundary conditions are computed by the udf-rotor_in_e3.lua and
udf-stator_out_e3.lua, which are called as user-defined boundary conditions. 


Author: Ingo Jahn
Last Modified: 17/05/2017
