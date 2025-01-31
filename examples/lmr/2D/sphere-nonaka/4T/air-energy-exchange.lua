-- Kinetics to mimic that used in Candler's DPLR code for nozzle flows
-- author: Robert Watt

Mechanism {
    "(*molcs) ~~ (*all)",
    type = "V-T",
    rate = "Landau-Teller",
    relaxation_time = {"Millikan-White"},
}

ChemistryCouplingMechanism{
    reaction_labels = {"r1","r2","r3","r4"},
    rate = "Marrone-Treanor",
    coupling_model = {model="ImpartialDissociation"}
}

-- remove this for 2T simulations
Mechanism{
   "(*molcs) ~~ (*molcs)",
   type = "V-V",
   rate = "Thivet-SSH",
   relaxation_time = {"Schwartz-Slawsky-Herzfeld"}
}
