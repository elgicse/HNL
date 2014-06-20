
from ROOT import *

def ComputeBRandLT(U2=6e-7):

    
    alpha_LT = 3*1e-6*6e-7
    alpha_BR = 3.6*1e-8/(6e-7)

    LT = alpha_LT/U2
    BR = alpha_BR*U2


    return {'LifeTime':LT, 'BR':BR}


def N(eff, U2=6e-7, fb=50., factor=1.0):

    N = (U2/6e-7)**2*eff*factor*20.*fb

    return N

#the neutrino lifetime is 3*10^(-6) s, not 3*10^(-7) s
#the branching of D->N+mu is 3*10^(-7), not 10^(-7)







