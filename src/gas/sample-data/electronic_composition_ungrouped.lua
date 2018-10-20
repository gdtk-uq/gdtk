-- Commenting nonsense at the top
-- some more
-- 
-- why not have a blank line
-- electronic gas model

-- First 8 are NI (1-46), second 8 are OI (1-40)
model = "ElectronicallySpecificGas"
number_electronic_species=91
grouped = "false"
s1 = 0.0
T1 = 300.0
p1 = 1e5

electronic_species = {}

electronic_species[0] = {
	name = 'NI - 1',
	individual_range = 1 - 1,
	M = 0.0140067,
	group_energy = 0.0,
	group_degeneracy = 4.0,
	dof = 3,
}

electronic_species[1] = {
	name = 'NI - 2',
	individual_range = 2 - 2,
	M = 0.0140067,
	group_energy = 2.384,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[2] = {
	name = 'NI - 3',
	individual_range = 3 - 3,
	M = 0.0140067,
	group_energy = 3.576,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[3] = {
	name = 'NI - 4',
	individual_range = 4 - 4,
	M = 0.0140067,
	group_energy = 10.332,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[4] = {
	name = 'NI - 5',
	individual_range = 5 - 5,
	M = 0.0140067,
	group_energy = 10.687,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[5] = {
	name = 'NI - 6',
	individual_range = 6 - 6,
	M = 0.0140067,
	group_energy = 10.927,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[6] = {
	name = 'NI - 7',
	individual_range = 7 - 7,
	M = 0.0140067,
	group_energy = 11.603,
	group_degeneracy = 2.0,
	dof = 3,
}

electronic_species[7] = {
	name = 'NI - 8',
	individual_range = 8 - 8,
	M = 0.0140067,
	group_energy = 11.759,
	group_degeneracy = 20.0,
	dof = 3,
}

electronic_species[8] = {
	name = 'NI - 9',
	individual_range = 9 - 9,
	M = 0.0140067,
	group_energy = 11.842,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[9] = {
	name = 'NI - 10',
	individual_range = 10 - 10,
	M = 0.0140067,
	group_energy = 11.996,
	group_degeneracy = 4.0,
	dof = 3,
}

electronic_species[10] = {
	name = 'NI - 11',
	individual_range = 11 - 11,
	M = 0.0140067,
	group_energy = 12.006,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[11] = {
	name = 'NI - 12',
	individual_range = 12 - 12,
	M = 0.0140067,
	group_energy = 12.125,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[12] = {
	name = 'NI - 13',
	individual_range = 13 - 13,
	M = 0.0140067,
	group_energy = 12.357,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[13] = {
	name = 'NI - 14',
	individual_range = 14 - 14,
	M = 0.0140067,
	group_energy = 12.856,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[14] = {
	name = 'NI - 15',
	individual_range = 15 - 15,
	M = 0.0140067,
	group_energy = 12.919,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[15] = {
	name = 'NI - 16',
	individual_range = 16 - 16,
	M = 0.0140067,
	group_energy = 12.972,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[16] = {
	name = 'NI - 17',
	individual_range = 17 - 17,
	M = 0.0140067,
	group_energy = 12.984,
	group_degeneracy = 28.0,
	dof = 3,
}

electronic_species[17] = {
	name = 'NI - 18',
	individual_range = 18 - 18,
	M = 0.0140067,
	group_energy = 13.0,
	group_degeneracy = 26.0,
	dof = 3,
}

electronic_species[18] = {
	name = 'NI - 19',
	individual_range = 19 - 19,
	M = 0.0140067,
	group_energy = 13.02,
	group_degeneracy = 20.0,
	dof = 3,
}

electronic_species[19] = {
	name = 'NI - 20',
	individual_range = 20 - 20,
	M = 0.0140067,
	group_energy = 13.035,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[20] = {
	name = 'NI - 21',
	individual_range = 21 - 21,
	M = 0.0140067,
	group_energy = 13.202,
	group_degeneracy = 2.0,
	dof = 3,
}

electronic_species[21] = {
	name = 'NI - 22',
	individual_range = 22 - 22,
	M = 0.0140067,
	group_energy = 13.245,
	group_degeneracy = 20.0,
	dof = 3,
}

electronic_species[22] = {
	name = 'NI - 23',
	individual_range = 23 - 23,
	M = 0.0140067,
	group_energy = 13.268,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[23] = {
	name = 'NI - 24',
	individual_range = 24 - 24,
	M = 0.0140067,
	group_energy = 13.294,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[24] = {
	name = 'NI - 25',
	individual_range = 25 - 25,
	M = 0.0140067,
	group_energy = 13.322,
	group_degeneracy = 4.0,
	dof = 3,
}

electronic_species[25] = {
	name = 'NI - 26',
	individual_range = 26 - 26,
	M = 0.0140067,
	group_energy = 13.343,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[26] = {
	name = 'NI - 27',
	individual_range = 27 - 27,
	M = 0.0140067,
	group_energy = 13.624,
	group_degeneracy = 12.0,
	dof = 3,
}

electronic_species[27] = {
	name = 'NI - 28',
	individual_range = 28 - 28,
	M = 0.0140067,
	group_energy = 13.648,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[28] = {
	name = 'NI - 29',
	individual_range = 29 - 29,
	M = 0.0140067,
	group_energy = 13.679,
	group_degeneracy = 90.0,
	dof = 3,
}

electronic_species[29] = {
	name = 'NI - 30',
	individual_range = 30 - 30,
	M = 0.0140067,
	group_energy = 13.693,
	group_degeneracy = 126.0,
	dof = 3,
}

electronic_species[30] = {
	name = 'NI - 31',
	individual_range = 31 - 31,
	M = 0.0140067,
	group_energy = 13.717,
	group_degeneracy = 24.0,
	dof = 3,
}

electronic_species[31] = {
	name = 'NI - 32',
	individual_range = 32 - 32,
	M = 0.0140067,
	group_energy = 13.77,
	group_degeneracy = 2.0,
	dof = 3,
}

electronic_species[32] = {
	name = 'NI - 33',
	individual_range = 33 - 33,
	M = 0.0140067,
	group_energy = 13.792,
	group_degeneracy = 38.0,
	dof = 3,
}

electronic_species[33] = {
	name = 'NI - 34',
	individual_range = 34 - 34,
	M = 0.0140067,
	group_energy = 13.824,
	group_degeneracy = 4.0,
	dof = 3,
}

electronic_species[34] = {
	name = 'NI - 35',
	individual_range = 35 - 35,
	M = 0.0140067,
	group_energy = 13.872,
	group_degeneracy = 10.0,
	dof = 3,
}

electronic_species[35] = {
	name = 'NI - 36',
	individual_range = 36 - 36,
	M = 0.0140067,
	group_energy = 13.925,
	group_degeneracy = 6.0,
	dof = 3,
}

electronic_species[36] = {
	name = 'NI - 37',
	individual_range = 37 - 37,
	M = 0.0140067,
	group_energy = 13.969,
	group_degeneracy = 18.0,
	dof = 3,
}

electronic_species[37] = {
	name = 'NI - 38',
	individual_range = 38 - 38,
	M = 0.0140067,
	group_energy = 13.988,
	group_degeneracy = 60.0,
	dof = 3,
}

electronic_species[38] = {
	name = 'NI - 39',
	individual_range = 39 - 39,
	M = 0.0140067,
	group_energy = 13.999,
	group_degeneracy = 126.0,
	dof = 3,
}

electronic_species[39] = {
	name = 'NI - 40',
	individual_range = 40 - 40,
	M = 0.0140067,
	group_energy = 14.054,
	group_degeneracy = 32.0,
	dof = 3,
}

electronic_species[40] = {
	name = 'NI - 41',
	individual_range = 41 - 41,
	M = 0.0140067,
	group_energy = 14.149,
	group_degeneracy = 18.0,
	dof = 3,
}

electronic_species[41] = {
	name = 'NI - 42',
	individual_range = 42 - 42,
	M = 0.0140067,
	group_energy = 14.16,
	group_degeneracy = 90.0,
	dof = 3,
}

electronic_species[42] = {
	name = 'NI - 43',
	individual_range = 43 - 43,
	M = 0.0140067,
	group_energy = 14.164,
	group_degeneracy = 126.0,
	dof = 3,
}

electronic_species[43] = {
	name = 'NI - 44',
	individual_range = 44 - 44,
	M = 0.0140067,
	group_energy = 14.202,
	group_degeneracy = 20.0,
	dof = 3,
}

electronic_species[44] = {
	name = 'NI - 45',
	individual_range = 45 - 45,
	M = 0.0140067,
	group_energy = 14.26,
	group_degeneracy = 108.0,
	dof = 3,
}

electronic_species[45] = {
	name = 'NI - 46',
	individual_range = 46 - 46,
	M = 0.0140067,
	group_energy = 14.316,
	group_degeneracy = 18.0,
	dof = 3,
}

electronic_species[46] = {
    name = 'NII',
    level = 1,
    individual_range = "all",
    M = 0.0140067,
    group_energy =  14.53413,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[47] = {
	name = 'OI - 1',
	individual_range = 1 - 1,
	M = 0.0159994,
	group_energy = 0.0,
	group_degeneracy = 9.0,
	dof = 3,
}

electronic_species[48] = {
	name = 'OI - 2',
	individual_range = 2 - 2,
	M = 0.0159994,
	group_energy = 1.97,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[49] = {
	name = 'OI - 3',
	individual_range = 3 - 3,
	M = 0.0159994,
	group_energy = 4.19,
	group_degeneracy = 1.0,
	dof = 3,
}

electronic_species[50] = {
	name = 'OI - 4',
	individual_range = 4 - 4,
	M = 0.0159994,
	group_energy = 9.146,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[51] = {
	name = 'OI - 5',
	individual_range = 5 - 5,
	M = 0.0159994,
	group_energy = 9.521,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[52] = {
	name = 'OI - 6',
	individual_range = 6 - 6,
	M = 0.0159994,
	group_energy = 10.74,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[53] = {
	name = 'OI - 7',
	individual_range = 7 - 7,
	M = 0.0159994,
	group_energy = 10.99,
	group_degeneracy = 9.0,
	dof = 3,
}

electronic_species[54] = {
	name = 'OI - 8',
	individual_range = 8 - 8,
	M = 0.0159994,
	group_energy = 11.838,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[55] = {
	name = 'OI - 9',
	individual_range = 9 - 9,
	M = 0.0159994,
	group_energy = 11.93,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[56] = {
	name = 'OI - 10',
	individual_range = 10 - 10,
	M = 0.0159994,
	group_energy = 12.09,
	group_degeneracy = 25.0,
	dof = 3,
}

electronic_species[57] = {
	name = 'OI - 11',
	individual_range = 11 - 11,
	M = 0.0159994,
	group_energy = 12.1,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[58] = {
	name = 'OI - 12',
	individual_range = 12 - 12,
	M = 0.0159994,
	group_energy = 12.3,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[59] = {
	name = 'OI - 13',
	individual_range = 13 - 13,
	M = 0.0159994,
	group_energy = 12.37,
	group_degeneracy = 9.0,
	dof = 3,
}

electronic_species[60] = {
	name = 'OI - 14',
	individual_range = 14 - 14,
	M = 0.0159994,
	group_energy = 12.55,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[61] = {
	name = 'OI - 15',
	individual_range = 15 - 15,
	M = 0.0159994,
	group_energy = 12.67,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[62] = {
	name = 'OI - 16',
	individual_range = 16 - 16,
	M = 0.0159994,
	group_energy = 12.71,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[63] = {
	name = 'OI - 17',
	individual_range = 17 - 17,
	M = 0.0159994,
	group_energy = 12.74,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[64] = {
	name = 'OI - 18',
	individual_range = 18 - 18,
	M = 0.0159994,
	group_energy = 12.76,
	group_degeneracy = 25.0,
	dof = 3,
}

electronic_species[65] = {
	name = 'OI - 19',
	individual_range = 19 - 19,
	M = 0.0159994,
	group_energy = 12.77,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[66] = {
	name = 'OI - 20',
	individual_range = 20 - 20,
	M = 0.0159994,
	group_energy = 12.78,
	group_degeneracy = 56.0,
	dof = 3,
}

electronic_species[67] = {
	name = 'OI - 21',
	individual_range = 21 - 21,
	M = 0.0159994,
	group_energy = 12.86,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[68] = {
	name = 'OI - 22',
	individual_range = 22 - 22,
	M = 0.0159994,
	group_energy = 12.89,
	group_degeneracy = 9.0,
	dof = 3,
}

electronic_species[69] = {
	name = 'OI - 23',
	individual_range = 23 - 23,
	M = 0.0159994,
	group_energy = 13.03,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[70] = {
	name = 'OI - 24',
	individual_range = 24 - 24,
	M = 0.0159994,
	group_energy = 13.05,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[71] = {
	name = 'OI - 25',
	individual_range = 25 - 25,
	M = 0.0159994,
	group_energy = 13.08,
	group_degeneracy = 40.0,
	dof = 3,
}

electronic_species[72] = {
	name = 'OI - 26',
	individual_range = 26 - 26,
	M = 0.0159994,
	group_energy = 13.087,
	group_degeneracy = 56.0,
	dof = 3,
}

electronic_species[73] = {
	name = 'OI - 27',
	individual_range = 27 - 27,
	M = 0.0159994,
	group_energy = 13.13,
	group_degeneracy = 15.0,
	dof = 3,
}

electronic_species[74] = {
	name = 'OI - 28',
	individual_range = 28 - 28,
	M = 0.0159994,
	group_energy = 13.14,
	group_degeneracy = 9.0,
	dof = 3,
}

electronic_species[75] = {
	name = 'OI - 29',
	individual_range = 29 - 29,
	M = 0.0159994,
	group_energy = 13.22,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[76] = {
	name = 'OI - 30',
	individual_range = 30 - 30,
	M = 0.0159994,
	group_energy = 13.23,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[77] = {
	name = 'OI - 31',
	individual_range = 31 - 31,
	M = 0.0159994,
	group_energy = 13.25,
	group_degeneracy = 168.0,
	dof = 3,
}

electronic_species[78] = {
	name = 'OI - 32',
	individual_range = 32 - 32,
	M = 0.0159994,
	group_energy = 13.33,
	group_degeneracy = 5.0,
	dof = 3,
}

electronic_species[79] = {
	name = 'OI - 33',
	individual_range = 33 - 33,
	M = 0.0159994,
	group_energy = 13.34,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[80] = {
	name = 'OI - 34',
	individual_range = 34 - 34,
	M = 0.0159994,
	group_energy = 13.353,
	group_degeneracy = 96.0,
	dof = 3,
}

electronic_species[81] = {
	name = 'OI - 35',
	individual_range = 35 - 35,
	M = 0.0159994,
	group_energy = 13.412,
	group_degeneracy = 8.0,
	dof = 3,
}

electronic_species[82] = {
	name = 'OI - 36',
	individual_range = 36 - 36,
	M = 0.0159994,
	group_energy = 13.418,
	group_degeneracy = 40.0,
	dof = 3,
}

electronic_species[83] = {
	name = 'OI - 37',
	individual_range = 37 - 37,
	M = 0.0159994,
	group_energy = 13.459,
	group_degeneracy = 8.0,
	dof = 3,
}

electronic_species[84] = {
	name = 'OI - 38',
	individual_range = 38 - 38,
	M = 0.0159994,
	group_energy = 13.464,
	group_degeneracy = 40.0,
	dof = 3,
}

electronic_species[85] = {
	name = 'OI - 39',
	individual_range = 39 - 39,
	M = 0.0159994,
	group_energy = 13.493,
	group_degeneracy = 3.0,
	dof = 3,
}

electronic_species[86] = {
	name = 'OI - 40',
	individual_range = 40 - 40,
	M = 0.0159994,
	group_energy = 13.496,
	group_degeneracy = 40.0,
	dof = 3,
}

electronic_species[87] = {
    name = 'OII',
    level = 1,
    individual_range = "all",
    M = 0.0159994,
    group_energy = 13.618055,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[88] = {
    name = 'e',
    level = 0,
    individual_range = "none",
    M = 0.0000005485799,
    group_energy = 0,
    group_degeneracy = 0,
    dof = 3,
}

electronic_species[89] = {
    name = "N2",
    level = 0,
    individual_range = "all",
    M = 0.028013,
    group_energy = 0.0,
    group_degeneracy = 0.0,
    dof = 5,
}

electronic_species[90] = {
    name = "O2",
    level = 0,
    individual_range = "all",
    M = 0.0319988,
    group_energy = 0.0,
    group_degeneracy = 0.0,
    dof = 5,
}