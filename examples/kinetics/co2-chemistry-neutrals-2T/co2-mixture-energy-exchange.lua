-- Author: Rowan J. Gollan
-- Date: 2021-06-12
--
-- Reference:
--   Park, Howe, Jaffe, and Candler (1994)
--   Review of Chemical-Kinetic Problems of Future NASA Missions,
--   II: Mars Entries
--   Journal of Thermophysics and Heat Transfer, 8(1), pp. 9--22
--
-- From Table 1

Mechanism{
   'CO ~~ (O, C)',
   type = 'V-T',
   rate = 'Landau-Teller',
   relaxation_time = {"Millikan-White", a=47.7, b=0.05}
}

Mechanism{
   'CO ~~ (CO, CO2)',
   type = 'V-T',
   rate = 'Landau-Teller',
   relaxation_time = {"Millikan-White"},
   note = "Ref states to use Eq.(3) which implementation does by default."
}

Mechanism{
   'CO2 ~~ (O, C, CO)',
   type = 'V-T',
   rate = 'Landau-Teller',
   relaxation_time = {"Millikan-White"},
   note = "Ref states to use Eq.(3) which implementation does by default."
}

Mechanism{
   'CO2 ~~ CO2',
   type = 'V-T',
   rate = 'Landau-Teller',
   relaxation_time = {"Millikan-White", a=36.5, b=-0.0193}
}

Mechanism{
   'O2 ~~ (*all)',
   type = 'V-T',
   rate = 'Landau-Teller',
   relaxation_time = {"Millikan-White"}
}

ChemistryCouplingMechanism{
    reaction_labels = {'r1','r2','r3','r4','r5','r6','r7'},
    rate = "Marrone-Treanor",
    coupling_model = {model="ImpartialDissociation"}
}
