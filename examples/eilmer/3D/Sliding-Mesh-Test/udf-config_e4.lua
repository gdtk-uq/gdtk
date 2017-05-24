-- udf-config_e4.lua
-- File to define Sliding Mesh Interface 

require 'lua_helper'

-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-- +++   The following is used to set up the sliding mesh   +++
-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



print("Hello from the Sliding Mesh Set-up, inside udf-config_e3.lua .")


-- Create Table with Simulation Conditions required by the udfs
N_BLADE = 16
OMEGA = -500/60 * 2 * math.pi -- in rad/s
theta_min = -0.5 * 2*math.pi/N_BLADE -- in rad
theta_max =  0.5 * 2*math.pi/N_BLADE -- in rad
relax_SR = 0.5   -- relaxtion factor applied to forwards flow
relax_RS = 0.5   -- relaxtion factor applied to backwards flow
               -- relax = 1. --> 100% new flow will be imported.  
               -- default is 0.5 <--> 0.5 which give correct inviscid solution 

-- the following are optional parameters for use with a user define outlet boundary  
P_out = 8.e6 -- oulet pressure in Pa
P_out_relax = 1. -- relaxation parameter for outlet pressure 
P_out_RC = 12.e6 -- oulet pressure in Pa
P_out_relax_RC = 1. -- relaxation parameter for outlet pressure 




-- Example 1: --> from Sliding mesh example
-- create lists of blocks that are positioned on upstream side of sliding interface
-- All faces must be WEST faces
-- Looking Radially outwards
--    Col1  Col0 
-- +------+------+
-- | blk1 | blk0 | ^
-- |      |      | | Z,k
-- +------+------+ |
--            <----+
--              T,j
upstream_row0 = {0,1} --blk0, blk1
UP_row_list = {upstream_row0}

-- create lists of blocks that are positioned on downstream side of sliding interface
-- All faces must be SOUTH faces
-- Looking Radially outwards
--      Col0 
-- +-------------+
-- |    blk3     |  
-- |             |  
-- +-------------+     
-- |    blk2     |  ^
-- |             |  | Z,k
-- +-------------+  |   
--            <-----+ 
--              T,i      

downstream_row0 = {2}  -- blk2
downstream_row1 = {3}  -- blk3
DOWN_row_list = {downstream_row0, downstream_row1}


-- Example 2: --> turbine
-- create lists of blocks that are positioned on upstream side of sliding interface
-- All faces must be WEST faces
-- Looking Radially outwards
--   Col4  Col3  Col2  Col1  Col0 
-- +-----+-----+-----+-----+-----+
-- |     |     |     |     |     |
-- |  12 |  13 |  14 |  15 |  16 |  Row0
-- |     |     |     |     |     |  ^
-- +-----+-----+-----+-----+-----+  | Z,k
--                           T <----+
--                               j
--
--upstream_row0 = {16,15,14,13,12} -- E4, E3, E2, E1, E0
--UP_row_list = {upstream_row0}
--
--
-- create lists of blocks that are positioned on downstream side of sliding interface
-- All faces must be SOUTH faces
-- Looking Radially outwards
--   Col1  Col0 
-- +-----+-----+
-- |     |     |
-- | IT0 | IT1 |  Row2
-- |  20 |  21 |  
-- +-----+-----+
-- |     |     |
-- | IC0 | IC1 |  Row1
-- |  23 |  24 |  
-- +-----+-----+
-- |     |     |
-- | IB0 | IB1 |  Row0
-- |  26 |  27 |  ^
-- +-----+-----+  | Z,k
--         T <----+
--              i  
--
--downstream_row0 = {27,26}  -- IB0, IB1
--downstream_row1 = {24,23}  -- IC0, IC1
--downstream_row2 = {21,20}  -- IT0, IT1
--DOWN_row_list = {downstream_row0, downstream_row1, downstream_row2}




-- do some data processing
LIST = {} -- create table containing all blocks that are part of interface
for i,row in ipairs(UP_row_list) do
    for j,b in ipairs(row) do
        table.insert(LIST,b)
    end
end
for i,row in ipairs(DOWN_row_list) do
    for j,b in ipairs(row) do
        table.insert(LIST,b)
    end
end


-- save data in data structure for subsequent processing
DATA = {}
DATA[0] = N_BLADE
DATA[1] = OMEGA
DATA[2] = UP_row_list
DATA[3] = DOWN_row_list
DATA[4] = LIST 
DATA[11] = theta_min
DATA[12] = theta_max   
DATA[13] = relax_RS  -- Relaxation for mapping from Rotor --> Stator
DATA[14] = relax_SR  -- Relaxation for mapping from Stator --> Rotor
DATA[15] = P_out            -- Optional Parameter for setting up rotor outlet
DATA[16] = P_out_relax      -- Optional Parameter for setting up rotor outlet
DATA[17] = P_out_RC         -- Optional Parameter for setting up rotor outlet
DATA[18] = P_out_relax_RC   -- Optional Parameter for setting up rotor outlet

-- write mapping tables to files  
assert( table.save( DATA, "CONFIG_tbl.lua" ) == nil )

print('Sliding Mesh Configuration file: CONFIG_tbl.lua  --> created')





