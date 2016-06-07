--Find the maximum divergence and its location as well as the RMS divergence in the Kelvin-Helmholtz MHD case. 

--Open file to write data to
local file = assert(io.open('plotdata.dat', 'wb'))
io.output(file)

--Number of timesteps to run through
nsteps = 100
for n = 0, nsteps do
	nb = 4
	fsol_clean = FlowSolution:new{jobName="KH-MHD-clean", dir="./KH-MHD-clean", tindx=n, nBlocks=nb}
	fsol_unclean = FlowSolution:new{jobName="KH-MHD-unclean", dir="./KH-MHD-unclean", tindx=n, nBlocks=nb}

	--Run through cells in each block
	first = true
	divB_RMS_clean = 0
	divB_RMS_unclean = 0
	NoCells = 0
	for ib = 0, nb-1 do
		--Cell domain is the same for clean and unclean
		local ni = fsol_clean:get_nic(ib)
		local nj = fsol_clean:get_njc(ib)
		for i = 0, ni-1 do
			for j = 0, nj-1 do
				cellData_clean = fsol_clean:get_cell_data{ib=ib, i=i, j=j}
				cellData_unclean = fsol_unclean:get_cell_data{ib=ib, i=i, j=j}
				--Find maximum divergence
				if (first) then
					divB_max_clean = cellData_clean["divB"]
					x_clean = cellData_clean["pos.x"]
					y_clean = cellData_clean["pos.y"]
				elseif (cellData_clean["divB"] > divB_max_clean) then
					divB_max_clean = cellData_clean["divB"]
					x_clean = cellData_clean["pos.x"]
					y_clean = cellData_clean["pos.y"]
				end
					
				if (first) then
					divB_max_unclean = cellData_unclean["divB"]
					x_unclean = cellData_unclean["pos.x"]
					y_unclean = cellData_unclean["pos.y"]
				elseif(cellData_unclean["divB"] > divB_max_unclean) then
					divB_max_unclean = cellData_unclean["divB"]
					x_unclean = cellData_unclean["pos.x"]
					y_unclean = cellData_unclean["pos.y"]
				end
				first = false

				--Add to divB_RMS
				divB_RMS_clean = divB_RMS_clean + cellData_clean["divB"]^2
				divB_RMS_unclean = divB_RMS_unclean + cellData_unclean["divB"]^2
			end
		end
		NoCells = NoCells + ni * nj
	end
	
	divB_RMS_clean = (divB_RMS_clean / NoCells)^0.5
	divB_RMS_unclean = (divB_RMS_clean / NoCells)^0.5
	--Write data to file
	io.write(n/nsteps," ", divB_max_clean," ", divB_RMS_clean," ", x_clean," ", y_clean," ", divB_max_unclean," ", divB_RMS_unclean," ", x_unclean," ", y_unclean, "\n")
end

