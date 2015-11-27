# make_profile.awk
# Set up an inflow profile for the inviscid vortex case
# PJ, 20-Feb-01, 14-Dec-06 write 1.0 for mass-fraction[0]
#
function pow( base, exponent ) {
   # print base, exponent
   return exp( exponent * log(base) )
}

BEGIN {
   Rgas  = 287      # J/kg.K
   g     = 1.4      # ratio of specific heats

   n     = 40
   r_i   = 1.0      # metres
   r_o   = 1.384
   dr    = (r_o - r_i) / n

   # Set flow properties ar the inner radius.
   p_i   = 100.0e3                  # kPa
   M_i   = 2.25
   rho_i = 1.0                      # kg/m**3
   T_i   = p_i / (Rgas * rho_i)     # K
   a_i   = sqrt( g * Rgas * T_i )   # m/s
   u_i   = M_i * a_i                # m/s
   # print p_i, M_i, rho_i, T_i, a_i, u_i

   # Generate the profile along the radial direction.
   print n > "profile.dat"
   for ( i = 1; i <= n; ++i ) {
      r   = r_i + dr * (i - 0.5)
      # print "i= ", i, "r=", r
      u   = u_i * r_i / r
      t1  = r_i / r
      t2  = 1.0 + 0.5 * (g - 1.0) * M_i * M_i * (1.0 - t1 * t1)
      rho = rho_i * pow( t2, 1.0/(g - 1.0) );
      p   = p_i * pow( rho/rho_i, g )
      T   = p / (rho * Rgas)
      # print p, u, 0.0, T, 1.0 > "profile.dat"
      print r/r_i, p/p_i, u/u_i, 0.0, T/T_i, 1.0 > "radial_profile_0.dat"
   } # end for
}
