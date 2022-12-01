"""
Python code for interacting with reaction mechanisms. To use this file,
first edit the gmodel_filename and rmechanism_filename to match your
simulation and then open the post-processed results in Paraview.

Select "Filters" -> "Alphabetical" -> "Programmable Filter", and paste this
code in the Script box.

(Ubuntu users may need to $ sudo apt install python3-paraview)

You will also need to build the gas libary using
$ cd gdtk/src/gas
$ make install

And update PYTHONPATH in your .bashrc:
export PYTHONPATH=${PYTHONPATH}:${DGD}/lib

@author: Nick Gibbons
"""

from gdtk.gas import GasModel, GasState, ReactionMechanism 

def get_array(name):
    for inp in inputs:
        if name in inp.CellData.keys():
            arr = inp.CellData[name]
            return arr
    else:
        raise KeyError('{} not found'.format(name))

def post_array(arr,name):
    print("Saving", name, "...")
    output.CellData.append(arr, name)
    return

gmodel_filename = "air-5sp-2T.gas"
rmechanism_filename =  "air-5sp-6r-2T.chem"

gmodel = GasModel(gmodel_filename)
rmech = ReactionMechanism(gmodel, rmechanism_filename)

nr = rmech.n_reactions
print("number of reactions is: ", nr)

# Note you may have to comment out the T_modes lines if you are working with a single T gas
T = get_array('T')
T_modes = get_array('T_modes[0]')
p = get_array('p')

massf = {}
for i,name in enumerate(gmodel.species_names):
	massf[name] = get_array('massf[{}]-{}'.format(i,name))

ncells = len(T)
reaction_tickrates = [zeros(ncells) for n in range(nr)]

state = GasState(gmodel)

for i in range(ncells):
    state.T = T[i]
    state.p = p[i]
    state.T_modes = [T_modes[i]]
    state.massf = {k:v[i] for k,v in massf.items()}
    state.update_thermo_from_pT()
    rates = rmech.reaction_tickrates(state)
    for j in range(nr): reaction_tickrates[j][i] = rates[j]

for j in range(nr):
    post_array(reaction_tickrates[j], 'rate_{}'.format(str(j).zfill(2)))

