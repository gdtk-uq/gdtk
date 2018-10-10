--Extract the final pressure, density and velocity from the shock tube

-- Invoked with the run_validation.sh script

--open file to write data to
local file = assert(io.open('ShockTubeNumerical.dat', 'wrb+'))
io.output(file)

--Number of blocks in domain
nb = 2

--Number of recorded timesteps
tsteps = 20

--Take the last time step
fsol = FlowSolution:new{jobName="MHDShockTube", dir=".", tindx = tsteps, nBlocks = nb}

R = 287

gamma = 5/3

for ib = 0, nb-1 do
	local ni = fsol:get_nic(ib)
	for i = 0, ni-1 do
		celldata = fsol:get_cell_data{ib=ib, i=i, j=1}
		x = celldata["pos.x"]
		velx = celldata["vel.x"]
		p = celldata["p"]
		T = celldata["T"]
		rho = p / (T * R)
		vely = celldata["vel.y"]
		velz = celldata["vel.z"]
		Bx = celldata["B.x"]
		By = celldata["B.y"]
		Bz = celldata["B.z"]
		psi = math.atan(Bz/By)
		E = p / (gamma - 1) + 0.5 * rho * (velx^2 + vely^2 + velz^2) + 0.5 * (Bx^2 + By^2 + Bz^2)
		secondderiv = 0
		deriv = 0
		if i >= 1 and i <= ni-2 then
			celldata_m1 = fsol:get_cell_data{ib=ib, i=i-1, j=1}
			x_m1 = celldata_m1["pos.x"]
			velx_m1 = celldata_m1["vel.x"]
			p_m1 = celldata_m1["p"]
			T_m1 = celldata_m1["T"]
			rho_m1 = p_m1 / (T_m1 * R)
	
			celldata_p1 = fsol:get_cell_data{ib=ib, i=i+1, j=1}
			x_p1 = celldata_p1["pos.x"]
			velx_p1 = celldata_p1["vel.x"]
			p_p1 = celldata_p1["p"]
			T_p1 = celldata_p1["T"]
			rho_p1 = p_p1 / (T_p1 * R)

			secondderiv = (rho_p1 - 2 * rho + rho_m1) / (x - x_m1) ^ 2
			deriv = 0.5 * (rho_p1 - rho_m1) / (x - x_m1)
		end		
		io.write(x, " ", rho, " ", p, " ", E, " ", velx, " ", vely, " ", velz, " ", By, " ", Bz, " ", psi, "\n")
	end
end


