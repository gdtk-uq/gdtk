reaction = {}
reaction[0]={
	equation = 'N2_eX1SIGGPlus_v0 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=6.13762e+17,  n=-0.674, C=113226.258, rctIndex=-1},
	brc={model='Arrhenius', A=2.28731e+18,  n=-1.171, C=27.496, rctIndex=-1},
	label=r0,
	type="dissociation-by-atom",
	molecule_idx = 1,
	atom_idx = 0
}

reaction[1]={
	equation = 'N2_eX1SIGGPlus_v1 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=1.02881e+18,  n=-0.730, C=109871.545, rctIndex=-1},
	brc={model='Arrhenius', A=3.86812e+18,  n=-1.228, C=23.563, rctIndex=-1},
	label=r1,
	type="dissociation-by-atom",
        molecule_idx = 2,
        atom_idx = 0
}

reaction[2]={
	equation = 'N2_eX1SIGGPlus_v2 + N_e2p3Minus4S <=> N_e2p3Minus4S + N_e2p3Minus4S + N_e2p3Minus4S',
	frc={model='Arrhenius', A=4.92300e+17,  n=-0.653, C=106565.747, rctIndex=-1},
	brc={model='Arrhenius', A=1.86716e+18,  n=-1.151, C=27.155, rctIndex=-1},
	label=r2,
	type="dissociation-by-atom",
        molecule_idx = 3,
        atom_idx = 0
}
