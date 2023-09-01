Lachlan Whyborn, 01/09/23 (Friday afternoon)

The Advecting Isentropic Vortex, a verification case for the Euler equations.
The vortex described by the initial conditions should advect with the prescribed
background velocity across the domain indefinitely.

The boundaries are periodic, so the domain is in effect infinite (although care
should be taken to make sure the vortex does not self-interfere across the
boundaries). See "A Survey of the Isentropic Euler Vortex Problem using High-Order
Methods" by Spiegel et al. for a thorough description.

Run the problem by calling ./run-IVA.lua (you might have to sudo chmod +x the file
first). The preamble in this file prescribes the range of flux calculators, 
interpolation order and cells to run for. The actual prep file is IVA.lua.

After calling run-IVA.lua, a new directory will be created for each flux calculator/
interpolation order combination. Within those directories are the results for each
discretization. For example, "ausmdv" flux calculation with interpolation order = 2
creates a directory called "ausmdv-InterpOrder-2.

If you want to scale up to higher cell counts, the number of blocks is controlled
by NB and configured for running in MPI by setting mpi = true. Number of processes
is the same as the number of blocks i.e. NB.

ComputeError.lua computes the L2 error, as well as the dissipation and dispersion
errors computed at each discretization using the method of Takacs in "A Two-step 
Scheme for the Advection Equation with Minimized Dissipation and Dispersion Errors",
 and written to ErrorMetrics.dat. ErrorMetrics.dat is placed in the relevant flux
calculator/interpolation order directory, to be plotted using gnuplot or some such.
