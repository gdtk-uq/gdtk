number_of_reactions = 2

reaction = {}
reaction[0]={
	equation = 'N2_eX1SIGGPlus_v0 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=4.15e22,  n=-1.5, C=1.131e5},
	brc={model='Arrhenius', A=2.32e21,  n=-1.5, C=0.},
	label=r0,
	type="dissociation-by-atom",
	molecule_idx = 1,
	atom_idx = 0
}

reaction[1]={
	equation = 'N2_eX1SIGGPlus_v1 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=4.15e22,  n=-1.5, C=1.131e5},
	brc={model='Arrhenius', A=2.32e21,  n=-1.5, C=0.},
	label=r1,
	type="dissociation-by-atom",
        molecule_idx = 2,
        atom_idx = 0
}


