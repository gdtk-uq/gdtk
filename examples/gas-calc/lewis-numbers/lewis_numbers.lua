-- A script to compute the effective Lewis numbers of high
-- temperature air using the binary diffusion routines.
--
-- Author: Nick Gibbons
-- Date: 2020-10-21
--
-- To run this calculation:
-- $ prep-gas air-5sp-1T.inp air-5sp-1T.lua
-- $ gas-calc lewis-numbers.lua
--

gasModelFile = 'air-5sp-1T.lua'
gmodel = GasModel:new{gasModelFile}

gs = GasState:new{gmodel}
gs.p = 1.0e5 -- Pa
gs.massf = {N2=0.70, O2=0.10, N=0.05, O=0.05, NO=0.1}
gs.T = 3000.0

gmodel:updateThermoFromPT(gs)
gmodel:updateTransCoeffs(gs)
D = gmodel:binary_diffusion_coefficients(gs)
molef = gmodel:massf2molef(gs)
cp = gmodel:Cp(gs)
print("Problem Description: ")
print(string.format("    T: %g", gs.T))
print(string.format("    p: %g", gs.p))
print(string.format("    cp: %g", cp))
print(string.format("    rho: %g", gs.rho))
print(string.format("    k: %g", gs.k))

-- Take the binary diffusion coefficients and compute a single 
-- average diffusion coefficient for each species.
-- Code from mass_diffusion.d
Dav = {}
for i,Di in pairs(D) do
    sum = 0.0
    for j,Dij in pairs(Di) do
        if (i~=j) then
            molefj = molef[gmodel:speciesName(j-1)]
            sum = sum + molefj/Dij
        end
    end
    molefi = molef[gmodel:speciesName(i-1)]
    Dav[i] = (1.0 - molefi)/sum
end

print("---------------------------------------")
print("Species Specific Lewis Numbers: ")
for i,Davi in pairs(Dav) do
    Lei = gs.rho*Davi*cp/gs.k
    print(string.format("    %s: %g", gmodel:speciesName(i-1), Lei))
end

print("Done.")
