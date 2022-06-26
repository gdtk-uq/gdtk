-- Set WITH_ELECTRONS appropriately so that electron energy exchange mechanisms 
-- are enabled if required.
WITH_ELECTRONS = false

Mechanism{
   "(*molcs) ~~ (*heavy)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"ParkHTC", submodel={"Millikan-White"}}
}

if WITH_ELECTRONS then

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

end
