-- post-processing script to compute the error metrics for the Isentropic vortex
-- Lachlan Whyborn 01/09/23

dofile("run-config.lua")

fsol_0 = FlowSolution:new{jobName = "IVA", dir = ".", nBlocks = NB, tindx = 0}
fsol_end = FlowSolution:new{jobName = "IVA", dir = ".", nBlocks = NB, tindx = 5}

-- This computes the L2 norm, as well as the dissipation/dispersion splitting
-- as described by Takacs in "A Two Step Scheme for the Advection Equation with
-- Minimized Dissipatio and Dispersion Errors", Monthly Weather Review vol. 113

L2 = 0; Exact_av = 0; Computed_av = 0; N_total = N * N
for ib = 0, NB-1 do
    nic = fsol_0:get_nic(ib); njc = fsol_0:get_njc(ib)
    for i = 0, nic-1 do
        for j = 0, njc-1 do
            ExactValue = fsol_0:get_cell_data{ib = ib, i = i, j = j}["rho"]
            ComputedValue = fsol_end:get_cell_data{ib = ib, i = i, j = j}["rho"]
            L2 = L2 + (ComputedValue - ExactValue)^2
            Exact_av = Exact_av + ExactValue
            Computed_av = Computed_av + ComputedValue
        end
    end
end

L2 = math.sqrt(L2 / N_total); Exact_av = Exact_av / N_total; Computed_av = Computed_av / N_total

-- Compute standard deviations
Exact_stdev = 0; Computed_stdev = 0
for ib = 0, NB-1 do
    nic = fsol_0:get_nic(ib); njc = fsol_0:get_njc(ib)
    for i = 0, nic-1 do
        for j = 0, njc-1 do
            ExactValue = fsol_0:get_cell_data{ib = ib, i = i, j = j}["rho"]
            ComputedValue = fsol_end:get_cell_data{ib = ib, i = i, j = j}["rho"]
            Exact_stdev = Exact_stdev + (Exact_av - ExactValue)^2
            Computed_stdev = Computed_stdev + (Computed_av - ComputedValue)^2
        end
    end
end

Exact_stdev = math.sqrt(Exact_stdev / N); Computed_stdev = math.sqrt(Computed_stdev / N)

-- Compute the correlation coefficient
CorrelationCoeff = 0
for ib = 0, NB-1 do
    nic = fsol_0:get_nic(ib); njc = fsol_0:get_njc(ib)
    for i = 0, nic-1 do
        for j = 0, njc-1 do
            ExactValue = fsol_0:get_cell_data{ib = ib, i = i, j = j}["rho"]
            ComputedValue = fsol_end:get_cell_data{ib = ib, i = i, j = j}["rho"]
            Exact_z = (ExactValue - Exact_av) / Exact_stdev
            Computed_z = (ComputedValue - Computed_av) / Computed_stdev
            CorrelationCoeff = CorrelationCoeff + Exact_z * Computed_z
        end
    end
end

CorrelationCoeff = CorrelationCoeff / N

DissipationError = math.sqrt((Exact_stdev - Computed_stdev)^2 + (Exact_av - Computed_av)^2)
DispersionError = math.sqrt(2 * (1 - CorrelationCoeff) * Exact_stdev * Computed_stdev)

f = assert(io.open("ErrorMetrics.dat", "a+"))
if first then
    f:write(string.format("# FluxCalc: %s InterpOrder: %d\n", flux_calculator, interpolation_order))
    f:write("# NCells\tL2\tDissipation\tDispersion\n")
end
f:write(string.format("%d\t%1.8e\t%1.8e\t%1.8e\n", N, L2, DissipationError, DispersionError))
f:close()
