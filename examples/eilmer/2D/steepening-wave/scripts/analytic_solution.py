"""
Get the exact solution for the steepening wave problem.

"Fluid Mechanics", L. D. Landau and E. M. Lifshitz
Chapter 10, Heinemann, 1987

@authors: Nick Gibbons and Lachlan Whyborn
"""

import numpy as np

def Bisection(f, a, b, tol):
    if f(a) * f(b) > 0:
        print("Guesses do not bracket the root or there are an even number of roots in the range")

    steps = 0

    while np.abs(a - b) > tol:
        c = (a + b) / 2.
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c

        steps = steps + 1

        if steps > 1000:
            print("Solution not converging")
            break

    return (a + b) / 2

class SteepeningWave(object):
    def __init__(self, x, rho_inf, p_inf, M_inf, Runi = 8.31446261815324, mMass = 0.02896000, gamma = 1.4):
        Rgas = Runi / mMass
        T_inf = p_inf / (Rgas * rho_inf);

        c0 = np.sqrt(gamma * Rgas * T_inf)
        u_inf = M_inf * c0
        
        self.x = x; self.rho_inf = rho_inf; self.p_inf = p_inf;
        self.Runi = Runi; self.mMass=mMass; self.gamma=gamma;
        self.Rgas = Rgas; self.T_inf = T_inf;
        self.c0 = c0; self.u_inf = u_inf;
        return

    def shock_formation_time(self):
        return np.abs(2.0/(self.u_inf * np.pi * (self.gamma + 1)))

    def SWP_solve(self, u_inf, x, t, gamma, c0):

        def f(u):
            return u / u_inf - np.sin(np.pi * (x - t * (c0 + 0.5 * (gamma + 1) * u)))
        
        def df(u):
            return 1. / u_inf - 0.5 * (gamma + 1) * t * np.pi * np.cos(np.pi * (x - t * (c0 + 0.5 * (gamma + 1) * u)))

        return Bisection(f, u_inf, -u_inf, 1e-12)

    def gen_profiles(self, t):
        """
        This function will compute the steepening wave profile at time t.
        """

        gamma = self.gamma
        n = len(self.x)
        u = np.empty(n)
        p = np.empty(n)
        rho = np.empty(n)

        for i in range(n):
            u[i] = self.SWP_solve(self.u_inf, self.x[i], t, gamma, self.c0)
            term = (1.0 + (gamma - 1.0) * u[i] / (2.0 * self.c0))
            p[i] = self.p_inf * term ** (2.0 * gamma / (gamma - 1.0))
            rho[i] = self.rho_inf * term ** (2.0 / (gamma - 1.0))
            
        return [p, rho, u]

if __name__=='__main__':
    x = np.linspace(0.0, 2.0)
    M_inf = -1.0
    p_inf = 1e5
    rho_inf = 1.0
    wave = SteepeningWave(x, rho_inf, p_inf, M_inf)
    t = 0.5*wave.shock_formation_time()
    p, rho, u = wave.gen_profiles(t)

    from pylab import plot,show
    plot(x, u)
    show()

