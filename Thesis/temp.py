from scipy import special
import numpy as np

eiTerms = np.loadtxt('temp.txt')

eiTerms = -1*special.expi(eiTerms)

np.savetxt('temp-solution.txt', eiTerms)

