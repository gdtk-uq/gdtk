ncellsCoarsest = 8
nLevels = 6
gridLevels = [None]*nLevels
for i in range(nLevels):
    ncells = 2**(nLevels-i-1)*ncellsCoarsest
    gridLevels[i] = { 'dx': 1.0/ncells,
                      'ncells': ncells }
