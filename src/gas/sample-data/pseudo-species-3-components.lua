-- Author: Rowan J. Gollan
-- Date: 2018-08-06
--
-- Demo file for defining pseudo-species components in a 
-- pseudo-species gas model.

number_pseudo_species = 3

pseudo_species = {}

pseudo_species[0] = {
   name = "N2_first_level",
   type = "single_state",
   energy = 0.0
}

pseudo_species[1] = {
   name = "N2_second_level",
   type = "single_state",
   energy = 3.6
}

pseudo_species[2] = {
   name = "N2_third_level",
   type = "single_state",
   energy = 4.2
}

