# initial-conditions.py
# 20th February 2018
# Jamie Border <jamie.border@uq.net.au>
# Solving for flow conditions for a desired shock wave Mach number.

class FlowState(object):
    def __init__(self, T1, p1, Ms):
        self._T1 = T1 	# K
        self._p1 = p1	# Pa
        self._Ms = Ms 	# -
        self._y  = 1.4	# -
        self._R  = 287 	# J / (kg.K)
        
        self.calc_states()
        self.return_states()

    def calc_states(self):
        y = self._y
        R = self._R
        Ms = self._Ms

        us = Ms * (y * R * self._T1)**(0.5)
        p2 = self._p1 * ((2 * y * Ms**2 - y + 1) / (y + 1))
        T2 = self._T1 * (((2 * y * Ms**2 - y + 1)
            * ((y - 1) * Ms**2 + 2))
            / ((y + 1)**2 * Ms**2))
        Mb = ((1 + (y - 1) / 2 * Ms**2) / (y * Ms**2 - (y - 1) / 2))**(0.5)
        u2 = us - Mb * (y * R * T2)**(0.5)
        M2 = u2 / ((y * R * T2)**(0.5))
        
        self._u1 = 0.0
        
        self._u2 = u2
        self._p2 = p2
        self._T2 = T2
        self._M2 = M2
        
    def return_states(self):

        print("Flow States: \n " +
                " ------------------------------\n" +
                "                s              \n" +
                "       (2)      s     (1)      \n" +
                "                s              \n" + 
                " -------------------           \n" +
                "                   |           \n" +
                "                   |           \n" +
                "                   ------------ \n")
        print("Conditions in (1):\n \tp1 = {0:3f} Pa\n\tT1 = {1:3f} K\n\tu1 = {2:3f} m/s".format(
				self._p1, self._T1, self._u1))
        print("Conditions in (2):\n\tp2 = {0:3f} Pa\n\tT2 = {1:3f} K\n\tu2 = {2:3f} m/s".format(
				self._p2, self._T2, self._u2))

if __name__ == "__main__":
    FS = FlowState(298.0, 1000.0, 5.09)
