-- Collision integral parameters fo Mars and Venus entries.
-- Sources:
--    Wright, Hwang and Schwenke (2007)
--    Recommended Collision Integrals for Transport Property Computations
--    Part 2: Mars and Venus Entries
--    AIAA Journal, 45(1), pp. 281--288
--
--    Wright, Bose, Palmer and Levin (2005)
--    Recommended Collision Integrals for Transport Property Computations
--    Part 1: Air species
--    AIAA Journal, 43(12), pp. 2558--2564
--
-- Data is fitted to Gupta-Yos expression by Dan Potter.

cis = {}
cis['CO2:CO2'] = {
   pi_Omega_11 = {A= -0.0115, B=  0.2886, C= -2.5580, D=  11.3378 },
   pi_Omega_22 = {A= -0.0148, B=  0.3731, C= -3.2655, D=  13.3972 }
}
cis['CO2:CO'] = {
   pi_Omega_11 = {A= -0.0078, B=  0.1780, C= -1.5361, D=  8.2160 },
   pi_Omega_22 = {A= -0.0119, B=  0.2770, C= -2.3226, D= 10.3647 }
}
cis['CO2:O2'] = {
   pi_Omega_11 = {A= -0.0068, B=  0.1637, C= -1.4935, D=  8.2050 },
   pi_Omega_22 = {A= -0.0110, B=  0.2678, C= -2.3236, D= 10.4795 }
}
cis['CO2:C'] = {
   pi_Omega_11 = {A= -0.0005, B= -0.0005, C= -0.1217, D=  4.4210 },
   pi_Omega_22 = {A= -0.0013, B=  0.0210, C= -0.2868, D=  4.9602 }
}
cis['CO2:O'] = {
   pi_Omega_11 = {A=  0.0002, B= -0.0134, C= -0.0206, D=  3.9632 },
   pi_Omega_22 = {A= -0.0009, B=  0.0120, C= -0.2145, D=  4.5663 }
}
cis['CO:CO'] = {
   pi_Omega_11 = {A= -0.0053, B=  0.0993, C= -0.8126, D=  6.0193},
   pi_Omega_22 = {A= -0.0081, B=  0.1664, C= -1.3274, D=  7.4106}
}
cis['CO:O2'] = {
   pi_Omega_11 = {A= -0.0031, B=  0.0602,  C= -0.5824, D=  5.5388 },
   pi_Omega_22 = {A= -0.0079, B=  0.1733,  C= -1.4406, D=  7.7641 }
}
cis['CO:C'] = {
   pi_Omega_11 = {A= -0.0019, B=  0.0294,  C= -0.3675, D=  4.9335 },
   pi_Omega_22 = {A= -0.0022, B=  0.0385,  C= -0.4262, D=  5.1878 }
}
cis['CO:O'] = {
   pi_Omega_11 = {A= -0.0013, B=  0.0152,  C= -0.2515, D=  4.4278 },
   pi_Omega_22 = {A= -0.0018, B=  0.0290,  C= -0.3565, D=  4.8211 }
}
cis['O2:O2'] = {
   pi_Omega_11 = {A= -0.0023, B=  0.0516,  C= -0.5785, D=  5.6041 },
   pi_Omega_22 = {A= -0.0089, B=  0.2066,  C= -1.7522, D=  8.6099 }
}
cis['O2:C'] = {
   pi_Omega_11 = {A= -0.0008, B=  0.0048,  C= -0.1689, D=  4.3108 },
   pi_Omega_22 = {A=  0.0002, B= -0.0159,  C= -0.0031, D=  4.0027 }
}
cis['O2:O'] = {
   pi_Omega_11 = {A= -0.0055, B=  0.1174,  C= -1.0770, D=  6.6923 },
   pi_Omega_22 = {A= -0.0049, B=  0.1020,  C= -0.9228, D=  6.3111 }
}
cis['C:C'] = {
   pi_Omega_11 = {A= -0.0138, B=  0.2950,  C= -2.3618, D= 10.0993 },
   pi_Omega_22 = {A= -0.0088, B=  0.1871,  C= -1.6061, D=  8.4436 }
}
cis['C:O'] = {
   pi_Omega_11 = {A= -0.0157, B=  0.3309,  C= -2.5544, D= 10.1976 },
   pi_Omega_22 = {A= -0.0125, B=  0.2665,  C= -2.1274, D=  9.3579 }
}
cis['O:O'] = {
   pi_Omega_11 = {A= -0.0055, B=  0.1245,  C= -1.2233, D=  7.2312 },
   pi_Omega_22 = {A= -0.0032, B=  0.0717,  C= -0.8058, D=  6.2443 }
}
