/* scale_pruett.c */
#include <stdio.h>

int main() {
   double y, T, ux;
   double d_star, T_inf, ux_inf, R0;

   /* Normalising quantities. */
   R0 = 0.005;  /* cylinder radius in metres */
   ux_inf = 597.3; /* external velocity */
   T_inf = 222.0; /* External temperature in Kelvin */
   d_star = 4.224e-3; /* displacement thickness */

   while (!feof(stdin)) {
      scanf( "%lf %lf %lf", &y, &T, &ux );
      y  = R0 + y * d_star;
      T  = T_inf * T;
      ux = ux_inf * ux;
      printf( "%g %g %g\n", y, T, ux );
   } /* end while */

   return 0;
} /* end function main */
