# extract_radial.awk
# Extract the radial profile data from e4shared --post generated files.
BEGIN{
   r_i = 1.0; p_i = 100.0e3; u_i = 841.87; T_i = 348.43;
}

NR > 1 {
   x = $1; y = $2; p = $9; u = $6; v = $7; T = $20
   r = sqrt( x * x + y * y )
   speed = sqrt( u * u + v * v )
   print r/r_i, p/p_i, speed/u_i, 0.0, T/T_i
}
