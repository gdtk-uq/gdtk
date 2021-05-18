Mechanism{
   "(*molcs) ~~ (*all)",
   type = "V-T",
   rate = "Landau-Teller",
   relaxation_time = {"Millikan-White"}
}

ChemistryCouplingMechanism{
    reaction_labels = {"r1","r2","r3","r4"},
    rate = "Marrone-Treanor",
    coupling_model = {model="ImpartialDissociation"}
}
