lumped_species = {}

lumped_species[0] = {
   name = 'N',
   type = 'atom',
   M = 0.0140067,
   DOF_base_mode = 3
}

lumped_species[1] = {
   name = 'N2',
   type = 'molecule',
   M = 0.0280134,
   DOF_base_mode = 5
}

pseudo_species = {}

pseudo_species[0] = {
   name = 'N_e2p3Minus4S',
   lumped_species_idx = 0,
   energy = 0.0
}

pseudo_species[1] = {
   name = 'N2_eX1SIGGPlus_v0',
   lumped_species_idx = 1,
   energy = 0.0
}

pseudo_species[2] = {
   name = 'N2_eX1SIGGPlus_v1',
   lumped_species_idx = 1,
   energy = 0.0
}

pseudo_species[3] = {
   name = 'N2_eX1SIGGPlus_v2',
   lumped_species_idx = 1,
   energy = 0.0
}


