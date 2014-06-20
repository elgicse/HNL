#decay = 'Ds->KX,X->PiPi'
decay = 'Ds->muN,N->muPi'


MeV = 1000.
'''Ds = 1.9685*MeV
D = 1.86962*MeV

massN = 1.6*MeV
massMu = 0.10565*MeV
massPi = 0.13957*MeV
massK = 0.493677*MeV
massX = 0.300*MeV
'''
decays = ['Ds->KX,X->PiPi','Ds->muN,N->muPi','D->piX,X->PiPi']

masses = {'N':1.6*MeV,
          'mu':0.10565*MeV,
          'pi':0.13957*MeV,
          'K':0.493677*MeV,
          'X':0.300*MeV,
          'Ds':1.9685*MeV,
          'D':1.86962*MeV}

mGluino = 0.


if not decay in decays:
    print 'Decay %s NOT found.'%decay
    print 'The possibilities are: %s'%decays
    assert(False)
elif decay == "Ds->muN,N->muPi":
    mother = "Ds"
    kid1 = "N"
    kid2 = "muTag"
    gkid1 = "pi"
    gkid2 = "mu"
    outName = "DsNmu,NmuPi"
    lifetime = 3e-6 #lifetimeN

elif decay == "Ds->KX,X->PiPi":
    mother = "Ds"
    kid1 = "X"
    kid2 = "K"
    gkid1 = "pi1"
    gkid2 = "pi2"
    outName = "DsXK,XPiPi" 
    lifetime = getLifetimeX(mGluino, F, mX)

elif decay == "Ds->piX,X->PiPi":
    mother = "D"
    kid1 = "N"
    kid2 = "piTag"
    gkid1 = "pi1"
    gkid2 = "pi2"
    outname = "DXpi,XPiPi"
    lifetime = getLifetimeX(mGluino, F, mX)

name2particle = {'pi':'pi','pi1':'pi', 'pi2':'pi','piTag':'pi',
                 'mu':'mu','muTag':'mu',
                 'K':'K',
                 'X':'X',
                 'N':'N',
                 'D':'D',
                 'Ds':'Ds'}


initialP = masses[name2particle[mother]]
# masses of the first dacay final state (e.g.: Ds -> mu N)
m1 = masses[name2particle[kid1]]
m2 = masses[name2particle[kid2]]
# masses of the second decay final state (e.g: N -> mu pi)
m3 = masses[name2particle[gkid1]]
m4 = masses[name2particle[gkid2]]
