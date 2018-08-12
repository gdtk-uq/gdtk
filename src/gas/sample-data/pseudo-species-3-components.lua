-- Author: Rowan J. Gollan
-- Date: 2018-08-06
--
-- Demo file for defining pseudo-species components in a 
-- pseudo-species gas model.

number_pseudo_species = 3

pseudo_species = {}

pseudo_species[0] = {
   name = 'N_e2p3Minus4S',
   type = "single_state",
   M = 0.0140067,
   DOF_base_mode = 3,
   energy = 0.0,
}

pseudo_species[1] = {
   name = 'N2_eX1SIGGPlus_v0',
   type = "single_state",
   M = 0.0280134,
   DOF_base_mode = 5,
   energy = -9.754, -- eV
}

pseudo_species[2] = {
   name = 'N2_eX1SIGGPlus_v1',
   type = "single_state",
   M = 0.0280134,
   DOF_base_mode = 5,
   energy = -9.465,
}

