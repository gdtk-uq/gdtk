C     ideal_air_fortran_demo.f
C     Exercise the functions from a fortran program.
C     PJ, 2016-12-27
      program iaf_demo
      use ideal_air_fortran
      implicit none

      real(kind=8) p, T, rho, u, massf(1)
      p = 100.0e3
      T = 300.0
      print '(a, e10.3, a, e10.3)', 'pressure=', p, ' temperature=', T

      end program

