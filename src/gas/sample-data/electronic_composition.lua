-- Commenting nonsense at the top
-- some more
-- 
-- why not have a blank line
-- electronic gas model

-- First 8 are NI (1-46), second 8 are OI (1-40)
number_electronic_species=21
s1 = 0.0
T1 = 300.0
p1 = 1e5

electronic_species = {}

electronic_species[0] = {
    name = 'NI 1',
    level = 1,
    individual_range = "1 - 1",
    M = 0.0140067,
    group_energy = 0,
    group_degeneracy = 4,
    dof = 3,
}

electronic_species[1] = {
    name = 'NI 2',
    level = 2,
    individual_range = "2 - 2",
    M = 0.0140067,
    group_energy = 2.384,
    group_degeneracy = 10,
    dof = 3,
}

electronic_species[2] = {
    name = 'NI 3',
    level = 3,
    individual_range = "3 - 3",
    M = 0.0140067,
    group_energy = 3.576,
    group_degeneracy = 6,
    dof = 3,
}

electronic_species[3] = {
    name = 'NI 4',
    level = 4,
    individual_range = "4 - 6",
    M = 0.0140067,
    group_energy = 10.641,
    group_degeneracy = 30,
    dof = 3,
}

electronic_species[4] = {
    name = 'NI 5',
    level = 5,
    individual_range = "7 - 13",
    M = 0.0140067,
    group_energy = 11.9508,
    group_degeneracy = 64,
    dof = 3,
}

electronic_species[5] = {
    name = 'NI 6',
    level = 6,
    individual_range = "14 - 21",
    M = 0.0140067,
    group_energy = 12.9848,
    group_degeneracy = 110,
    dof = 3,
}

electronic_species[6] = {
    name = 'NI 7',
    level = 7,
    individual_range = "22 - 27",
    M = 0.0140067,
    group_energy = 13.342,
    group_degeneracy = 64,
    dof = 3,
}

electronic_species[7] = {
    name = 'NI 8',
    level = 8,
    individual_range = "28-46",
    M = 0.0140067,
    group_energy = 13.9876,
    group_degeneracy = 922,
    dof = 3,
}

electronic_species[8] = {
    name = 'NII',
    level = 1,
    individual_range = "all",
    M = 0.0140067,
    group_energy =  14.53413,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[9] = {
    name = 'OI 1',
    level = 1,
    individual_range = "1 - 1",
    M = 0.0159994,
    group_energy = 0,
    group_degeneracy = 9,
    dof = 3,
}

electronic_species[10] = {
    name = 'OI 2',
    level = 2,
    individual_range = "2 - 2",
    M = 0.0159994,
    group_energy = 1.97,
    group_degeneracy = 5,
    dof = 3,
}

electronic_species[11] = {
    name = 'OI 3',
    level = 3,
    individual_range = "3 - 3",
    M = 0.0159994,
    group_energy = 4.19,
    group_degeneracy = 1,
    dof = 3,
}

electronic_species[12] = {
    name = 'OI 4',
    level = 4,
    individual_range = "4 - 6",
    M = 0.0159994,
    group_energy = 10.2345,
    group_degeneracy = 23,
    dof = 3,
}

electronic_species[13] = {
    name = 'OI 5',
    level = 5,
    individual_range = "7 - 13",
    M = 0.0159994,
    group_energy = 12.0181,
    group_degeneracy = 81,
    dof = 3,
}

electronic_species[14]  ={
    name = 'OI 6',
    level = 6,
    individual_range = "14 - 21",
    M = 0.0159994,
    group_energy = 12.7522,
    group_degeneracy = 139,
    dof = 3,
}

electronic_species[15] = {
    name = 'OI 7',
    level = 7,
    individual_range = "22 - 27",
    M = 0.0159994,
    group_energy = 13.0729,
    group_degeneracy = 128,
    dof = 3,
}

electronic_species[16] = {
    name = 'OI 8',
    level = 8,
    individual_range = "28 - 40",
    M = 0.0159994,
    group_energy = 13.3392,
    group_degeneracy = 428,
    dof = 3,
}

electronic_species[17] = {
    name = 'OII',
    level = 1,
    individual_range = "all",
    M = 0.0159994,
    group_energy = 13.618055,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[18] = {
    name = 'free electron',
    level = 0,
    individual_range = "none",
    M = 0.0000005485799,
    group_energy = 0,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[19] = {
    name = "N2",
    level = 0,
    individual_range = "all",
    M = 0.028013,
    group_energy = 0.0,
    group_degeneracy = 0.0,
    dof = 5,
}

electronic_species[20] = {
    name = "O2",
    level = 0,
    individual_range = "all",
    M = 0.0319988,
    group_energy = 0.0,
    group_degeneracy = 0.0,
    dof = 5,
}