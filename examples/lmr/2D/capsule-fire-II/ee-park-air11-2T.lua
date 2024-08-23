-- VT Exchange Rates from Table 1 from Park, 1993
-- Some of these match the Millian White a and b values, but some don't
-- Ions are assumed to be the same as their neutral molecules

------------------ N2 ------------------
Mechanism{
   "(N2) ~~ (N,N+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=180.0, b = 0.0262}}
}

Mechanism{
   "(N2) ~~ (O,O+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=72.4, b = 0.0150}}
}

Mechanism{
   "(N2) ~~ (N2,N2+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=221.0, b = 0.0290}}
}

Mechanism{
   "(N2) ~~ (O2,O2+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=229.0, b = 0.0295}}
}

Mechanism{
   "(N2) ~~ (NO,NO+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=225.0, b = 0.0293}}
}

------------------ O2 ------------------
Mechanism{
   "(O2) ~~ (N,N+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=72.4, b = 0.015}}
}

Mechanism{
   "(O2) ~~ (O,O+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=47.7, b = 0.059}}
}

Mechanism{
   "(O2) ~~ (N2,N2+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=134.0, b = 0.0295}}
}

Mechanism{
   "(O2) ~~ (O2,O2+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=138.0, b = 0.0300}}
}

Mechanism{
   "(O2) ~~ (NO,NO+)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=136.0, b = 0.0298}}
}

------------------ NO ------------------

-- It's worth nothing that Gehre, 2012 objects to this treatment of the NO colliders being all the same
Mechanism{
   "(NO) ~~ (*heavy)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC2", submodel={"Millikan-White", a=49.5, b=0.042}}
}


-- ET Rates from Gnoffo 1989, since Park doesn't mention these and they're important for FireII
Mechanism{
    "(e-) ~~ (N2)",
    type = "E-T",
    exchange_cross_section = {type="GnoffoNeutral", a = 7.5e-20, b = 5.5e-24, c = -1.0e-28}
}

Mechanism{
    "(e-) ~~ (O2)",
    type = "E-T",
    exchange_cross_section = {type="GnoffoNeutral", a = 2.0e-20, b = 6.0e-24, c = 0.0e-00}
}

Mechanism{
    "(e-) ~~ (N)",
    type = "E-T",
    exchange_cross_section = {type="GnoffoNeutral", a = 5.0e-20, b = 0.0e+00, c = 0.0e+00}
}

Mechanism{
    "(e-) ~~ (O)",
    type = "E-T",
    exchange_cross_section = {type="GnoffoNeutral", a = 1.2e-20, b = 1.7e-24, c = -2.0e-29}
}
    
Mechanism{
    "(e-) ~~ (NO)",
    type = "E-T",
    exchange_cross_section = {type="GnoffoNeutral", a = 1.0e-19, b = 0.0e+00, c = 0.0e+00}
}

Mechanism{
    "(e-) ~~ (*ions)",
    type = "E-T",
    exchange_cross_section = {type="Coulomb"}
}
