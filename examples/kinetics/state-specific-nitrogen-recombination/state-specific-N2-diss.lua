number_of_reactions = 2

test = "dissociation-by-atom"

if test == "dissociation-by-atom" then

reaction = {}
reaction[0]={
	equation = 'N2_eX1SIGGPlus_v0 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05},
	brc={model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00},
	label=r0,
	type="dissociation-by-atom",
	molecule_idx = 1,
	atom_idx = 0
}

reaction[1]={
	equation = 'N2_eX1SIGGPlus_v1 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05},
	brc={model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00},
	label=r1,
	type="dissociation-by-atom",
        molecule_idx = 2,
        atom_idx = 0
}

elseif test == "collision-2B-3B" then

reaction = {}
reaction[0]={
	equation = 'N2_eX1SIGGPlus_v0 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05},
	brc={model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00},
	label=r0,
	type="collision-2B-3B",
	reactants_idx = {1, 0},
	products_idx = {0, 0, 0},
}

reaction[1]={
	equation = 'N2_eX1SIGGPlus_v1 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05},
	brc={model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00},
	label=r1,
	type="collision-2B-3B",
	reactants_idx = {2, 0},
	products_idx = {0, 0, 0},
}

end


