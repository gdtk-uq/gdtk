-- Compare two solutions using a custom post-processing script
-- @author: Nick Gibbons

for blockfiles in io.popen("ls -1 ./flow/t0000 | wc -l"):lines() do
   nblocks = tonumber(blockfiles)
end

fsol0 = FlowSolution:new{jobName="swbli", dir=".", tindx=2, nBlocks=nblocks}
fsol6 = FlowSolution:new{jobName="swbli", dir=".", tindx=3, nBlocks=nblocks}

pdiff = 0.0

for ib=0,nblocks-1 do
   local ni = fsol0:get_nic(ib)
   local nj = fsol0:get_njc(ib)
   for ic=0,ni-1 do
      for jc=0,nj-1 do
         cellData0 = fsol0:get_cell_data{ib=ib, i=ic, j=jc}
         cellData6 = fsol6:get_cell_data{ib=ib, i=ic, j=jc}
         p0 = cellData0['p']
         p6 = cellData6['p']
         pdiff = pdiff + (p6-p0)*(p6-p0)
      end
   end
end

pdiff = math.sqrt(pdiff)
print("pdiff: ", pdiff)
