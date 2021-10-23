Summary:
  Reflected shock tube analysis followed by supersonic nozzle expansion.
  There is a stand-alone program nenzf1d and a Python-callable library.

Authors:
  Peter J., Nick Gibbons and Rowan Gollan.
  School of Mechanical and Mining Engineering
  The University of Queensland

Versions:
  2020-09-26 Initial code built from Python prototype.
  2020-10-16 Ready for use with 1T thermally-perfect gas model, I believe.
  2020-12-09 Allow 2T gas models.
  2021-09-03 Refactor by NNG
  2021-10-23 Refactor Nick's refactor.

Notes:
  Starting from a known initial temperature, pressure and equilibrium composition,
  it is shock compressed by an incident shock, brought to rest by a reflected shock,
  and allowed to expand to an observed stagnation pressure in the nozzle-supply
  region of the shock tube.  Equilibrium thermochemistry is assumed for these states.
  The equilibrium gas is allowed to expand to slightly supersonic conditions near
  the start of the diverging (supersonic) nozzle expansion.
  With the change to a nonequilibrium thermochemical model for the nozzle expansion,
  there is a slight jump in sound-speed for the same nominal conditions, so it may
  be necessary to run the equilibrium expansion beyond Mach 1 to be sure that the
  expansion stepping with nonequilibrium thermochemistry takes the desired supersonic branch.
  We have found that the 2T gas models require a transition Mach number up to 1.2
  in order to expand supersonically.  1T thermochemical models seem to work well enough
  with transition Mach numbers much closer to 1.
