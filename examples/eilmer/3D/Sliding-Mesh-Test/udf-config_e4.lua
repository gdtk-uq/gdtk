-- udf-config_e3.lua
-- File to define Sliding Mesh Interface 

-- Define Functions required by lua set-up script
--// exportstring( string )
--// returns a "Lua" portable version of the string
local function exportstring( s )
   return string.format("%q", s)
end

--// The Save Function
   function table.save(  tbl,filename )
      local charS,charE = "   ","\n"
      local file,err = io.open( filename, "wb" )
      if err then return err end

      -- initiate variables for save procedure
      local tables,lookup = { tbl },{ [tbl] = 1 }
      file:write( "return {"..charE )

      for idx,t in ipairs( tables ) do
         file:write( "-- Table: {"..idx.."}"..charE )
         file:write( "{"..charE )
         local thandled = {}

         for i,v in ipairs( t ) do
            thandled[i] = true
            local stype = type( v )
            -- only handle value
            if stype == "table" then
               if not lookup[v] then
                  table.insert( tables, v )
                  lookup[v] = #tables
               end
               file:write( charS.."{"..lookup[v].."},"..charE )
            elseif stype == "string" then
               file:write(  charS..exportstring( v )..","..charE )
            elseif stype == "number" then
               file:write(  charS..tostring( v )..","..charE )
            end
         end

         for i,v in pairs( t ) do
            -- escape handled values
            if (not thandled[i]) then
            
               local str = ""
               local stype = type( i )
               -- handle index
               if stype == "table" then
                  if not lookup[i] then
                     table.insert( tables,i )
                     lookup[i] = #tables
                  end
                  str = charS.."[{"..lookup[i].."}]="
               elseif stype == "string" then
                  str = charS.."["..exportstring( i ).."]="
               elseif stype == "number" then
                  str = charS.."["..tostring( i ).."]="
               end
            
               if str ~= "" then
                  stype = type( v )
                  -- handle value
                  if stype == "table" then
                     if not lookup[v] then
                        table.insert( tables,v )
                        lookup[v] = #tables
                     end
                     file:write( str.."{"..lookup[v].."},"..charE )
                  elseif stype == "string" then
                     file:write( str..exportstring( v )..","..charE )
                  elseif stype == "number" then
                     file:write( str..tostring( v )..","..charE )
                  end
               end
            end
         end
         file:write( "},"..charE )
      end
      file:write( "}" )
      file:close()
   end
   
--// The Load Function
function table.load( sfile )
   local ftables,err = loadfile( sfile )
   if err then return _,err end
   local tables = ftables()
   for idx = 1,#tables do
      local tolinki = {}
      for i,v in pairs( tables[idx] ) do
         if type( v ) == "table" then
            tables[idx][i] = tables[v[1]]
         end
         if type( i ) == "table" and tables[i[1]] then
            table.insert( tolinki,{ i,tables[i[1]] } )
         end
      end
      -- link indices
      for _,v in ipairs( tolinki ) do
         tables[idx][v[2]],tables[idx][v[1]] =  tables[idx][v[1]],nil
      end
   end
   return tables[1]
end


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





