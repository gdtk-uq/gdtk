# t1d_classes.py
# Cell and Face classes for constructing a quasi-one-dimensional
# simulation of gas flow in a duct.
#
# PJ, 2020-02-12
#
from eilmer.gas import GasModel, GasState

class Cell(object):
    __slots__ = ['x', 'vol', 'gas', 'vel', 'mass', 'momentum', 'energy']

    def __init__(self, gmodel, x=0.0, vol=1.0, p=1.0e5, T=300.0, vel=0.0):
        self.x = x
        self.vol = vol
        self.gas = GasState(gmodel)
        self.gas.p = p
        self.gas.T = T
        self.gas.update_thermo_from_pT()
        self.vel = vel
        self.encode_conserved()

    def encode_conserved(self):
        # Conserved quantities are per unit volume.
        self.mass = self.gas.rho
        self.momentum = self.gas.rho * self.vel
        self.energy = self.gas.rho * (self.gas.u + 0.5*self.vel*self.vel)
        return

    def decode_conserved(self):
        self.gas.rho = self.mass
        self.vel = self.momentum / self.gas.rho
        self.gas.u = (self.energy / self.gas.rho) - 0.5*self.vel*self.vel
        self.gas.update_thermo_from_rhou()
        return

    def __repr__(self):
        text = "Cell({}, {}, {}, {}, {})"
        return text.format(self.x, self.vol, self.gas.p, self.gas.T, self.vel)

class Face(object):
    __slots__ = ['x', 'area', 'gas', 'vel']

    def __init__(self, gmodel, x=0.0, area=1.0, p=1.0e5, T=300.0, vel=0.0):
        self.x = x
        self.area = area
        self.gas = GasState(gmodel)
        self.gas.p = p
        self.gas.T = T
        self.gas.update_thermo_from_pT()
        self.vel = vel

    def __repr__(self):
        text = "Face({}, {}, {}, {}, {})"
        return text.format(self.x, self.area, self.gas.p, self.gas.T, self.vel)
