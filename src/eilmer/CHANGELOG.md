# Changelog

This changelog is focused on changes related to our tagged releases.
The intent is to document all notable changes, with most recent listed first.

## [Unreleased]

- shift default interpolation type for StructuredGrids from "simple" to "linear"

## [4.1.0] - 2022-02-23
### Meta Change
- Repository hosting shifted from bitbucket to github.
  The repository is now located at: <https://github.com/gdtk-uq/gdtk>
  
### Changed
- Default install location `$HOME/dgdinst` --> `$HOME/gdtkinst`
- Specification of `cfl_schedule` in transient code is now of the form:

        config.cfl_schedule = {{0.0,0.5}, {0.5e-3,1.0}, {2.0e-3,10.0}}
	
  In each pair, the first number is time and the second number is CFL value.
- For specialist post-processing, class `BlockFlow` renamed to `FluidBlockLite`
  to better represent its intent. `FluidBlockLite` is used in post-processing
  to read, interact, manipulate fluid data at a block level.
  The "lite" refers to the fact that it is not a full-blown `FluidBlock` object
  as implemented in the code internals.
- `prep-gas` now selects the `CompositeGas` when requesting a thermally perfect gas model.
  This uses the the unified gas model machinery built in 2021, and begins the
  deprecation of the special-purpose implementations.
- Lua version embedded in code upgraded: 5.1.4 --> 5.4.3

### Fixed
- Clean up on MPI exit paths to avoid hanging processes when the gas
dynamics encounters adverse situations.
- Fixed bug in computation of transport properties and update of sound speed
  at a boundary interface for the special user-defined boundary interface effect.

### Added
- Control point form as option for algebraic grid generation on structured grids
  See technical note: <https://gdtk.uqcloud.net/docs/eilmer/technical-notes/control-point-grid-gen/>
- Conical inflow as option to supersonic inflow, constant flux and shock-fitting boundary conditions.
  See example in: `examples/eilmer/2D/nozzle-conical-ideal`
- `GridArray` class added to separate the role of grid division out from `FBArray`.
- `FBArray` can now be constructed from an array of grid objects.
  See example: `examples/eilmer/2D/shear-layer-periodic/psl-gridArray.lua`
- An option to compute mixture averaged diffusion coefficients from the binary
  diffusion coefficients of all interacting pairs in the mixture.
- A more robust calculation of wall distances for use in the Spalart-Allmaras turbulence model.
- `e4monitor` program now has option to run in foreground mode to start a simulation:

        $ e4monitor --job=cone20 --startup=1 --period=1 --command="mpirun -np 8 --oversubscribe e4mpi --run --job=cone20"

- `CubePatch` geometry object. `CubePatch` is like `SpherePatch` but produces a flat (square) surface with edge length 2a.
- Meta-data at Eilmer start-up now include more information: Revision-date, Build-date and Profiling.
- `RuledSurface` geometry object.

