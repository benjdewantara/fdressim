# this file is just a demonstration to show us that
# for a slightly compressible fluid
# initial pressure distribution with respect to z does not vary much
# so it is safe to assume that with respect to z, we can ignore gravity effect

from scipy import integrate
import numpy as np



#rhoGDeltaZDimMultiplier = 0.3048/(144*9.80665)
rhoGDeltaZDimMultiplier = 1
class Fluid(object):
    def __init__(self, refRho, refPres, compress, mu):
        '''
        A Fluid object will represent a fluid residing in a reservoir.
        The properties of a fluid we are interested in are:
        density (self.rho, in lbm/ft^3), and
        viscosity (self.mu, in cP).
        The behavior of density with respect to pressure is explained using
        compressibility (compressibility() function, in psi^-1)
        
        '''
        self.refRho = refRho
        self.refPres = refPres
        self.compress = compress
        self.mu = mu
    
    def getRho(self, pres):
        return self.refRho*np.exp(self.compress*(pres - self.refPres))
    





f = Fluid(refRho=62.428, refPres=14.7, compress=3.5*1e-6, mu=10)


gDeltaZ = 0.3048/(144*9.80665) * 32.1740485 * (75*75)

integrateLHS = integrate.quad(f.getRho, 6000, 6000.5)


print("integrateLHS[0]=%s gDeltaZ=%s" %(integrateLHS[0], gDeltaZ))














