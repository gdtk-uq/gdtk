#!/usr/bin/env dgd-lua

-- Each of these should be specified as tables
N = {64, 128, 256}                      -- Number of cells per side
flux_calculator = {"ldfss2", "ausmdv", "roe"}  -- Which flux calculator
interpolation_order = {1, 2, 3}         -- Order of the interpolation

-- This is single-valued- blocking structure
NB = 16                                 -- Number of fluid blocks- must be divisible by 4
assert(NB % 4 == 0, "NB must be divisible by 4")
mpi = true                              -- Set up the MPI mapping

os.execute("prep-gas ideal-air.inp ideal-air.lua")
for _, FluxCalc in ipairs(flux_calculator) do
    for _, InterpOrder in ipairs(interpolation_order) do
        jobString = FluxCalc.."-InterpOrder-"..InterpOrder
        os.execute("mkdir -p "..jobString)

        -- Remove old data
        os.execute("rm -r "..jobString)

        -- A flag for the writing errors to file
        first = true
        for _, NCells in ipairs(N) do

            print("\n-----------------------------------------------\n")
            print("Running Isentropic Vortex with:\nFlux Calculator: "..FluxCalc.."\nInterpolation Order: "..InterpOrder.."\nCells: "..NCells.."x"..NCells)

            os.execute("mkdir -p "..jobString.."/"..NCells.."x"..NCells)

            -- Create the prep file to be read by IVA.lua
            f = assert(io.open("run-config.lua", "w+"))
            f:write(string.format("N = %d; flux_calculator = \"%s\"; interpolation_order = %d; NB = %d; mpi = %s; first = %s", NCells, FluxCalc, InterpOrder, NB, mpi, first))
            f:close()

            -- Do prep and run, and catch for failures
            assert(os.execute("e4shared --prep --job=IVA"), "Prep stage failed")
            if mpi then
                assert(os.execute("mpirun -np "..NB.." e4mpi --run --job=IVA"), "Run stage failed")
            else
                assert(os.execute("e4shared --run --job=IVA"), "Run stage failed")
            end

            -- Compute error metrics
            assert(os.execute("e4shared --custom-post --script-file=ComputeError.lua"), "Post stage failed")
            first = false
            
            os.execute("mv flow grid config "..jobString.."/"..NCells.."x"..NCells)
            os.execute("rm run-config.lua")
        end
        os.execute("mv ErrorMetrics.dat "..jobString)
    end
end




