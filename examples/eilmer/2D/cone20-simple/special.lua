-- special.lua
-- Try out custom post-processing for Eilmer4
print("Begin custom script")
fsol = FlowSolution:new{jobName="cone20", dir=".", tindx=1, nBlocks=2}
print("fsol=", fsol)
nc = fsol:find_nearest_cell_centre{x=0.5, y=0.5}
print("closest cell indices=", nc["ib"], nc["i"], nc["j"], nc["k"])
print("for block 1, nic=", fsol:get_nic(1),
      " njc=", fsol:get_njc(1), 
      " nkc=", fsol:get_nkc(1))
vtx = fsol:get_vtx{ib=1, i=13, j=18}
print("vtx=", vtx)
varNames = fsol:get_var_names()
print("varNames=")
for i,v in pairs(varNames) do
   print("   ", i, v)
end
cellData = fsol:get_cell_data{ib=1, i=13, j=18}
print("cellData=")
for k,v in pairs(cellData) do
   print("   ", k, "=", v)
end
print("Finished custom script")
