

from ROOT import *


#lifetime = 0.0018
#lifetime = 1.3e-5
#lifetime = 3.e-6


#filename = "Neutrinos_acceptance_LT_%s.root"%(lifetime)
#filename = "NeutrinosLHC_acceptance_LT_3e-06.root"


filename = "D2Nmu_M1.6_LT3e-06.root"

f = TFile(filename)
t = f.Get("tree_acceptance")

nentries = t.GetEntries()





##############################################################################
def Get4Vector(name, t):
    X = t.__getattr__(name+"_Px")
    Y = t.__getattr__(name+"_Py")
    Z = t.__getattr__(name+"_Pz")
    E = t.__getattr__(name+"_E")    
    lorentz = TLorentzVector(Z,X,Y,E)
    
    return lorentz

def Get3Vector(name, t, scale=10):
    X = t.__getattr__(name+"x")/scale
    Y = t.__getattr__(name+"y")/scale
    Z = t.__getattr__(name+"z")/scale  
    vect = TVector3(Z,X,Y)
    return vect


def GetXY(Z, Origin, direction):

    Lambda =  (Z-Origin[0])/direction[0]
    X = Origin[1]+Lambda*direction[1]
    Y = Origin[2]+Lambda*direction[2]

    return {'X':X, 'Y':Y}
    

def IsInAcceptance(R0, Impact):

    X = Impact['X']/100.
    Y = Impact['Y']/100.
    R = TMath.sqrt(X**2+Y**2)
    if R<R0:
        return True
    return False


############ Check if the decay is in the acceptance of the first element
def InTheFirstElement(pmu, ppi, Origin):

    #### I first check that the decay is inside the first vacuum tank
    if Origin[0]/100>60. and Origin[0]/100.<100.:        
        if TMath.sqrt((Origin[1]/100.)**2+(Origin[2]/100.)**2)<2.5:
            ### Now I check that the decay product go through the magnet
            #return True
            ##### The end of the maget is located at 95+5 m (probably) 
            impact_mu = GetXY(10000., Origin, pmu)
            impact_pi = GetXY(10000., Origin, ppi)

            Rmu = TMath.sqrt((impact_mu["X"]/100.)**2+(impact_mu["Y"]/100.)**2)
            Rpi = TMath.sqrt((impact_pi["X"]/100.)**2+(impact_pi["Y"]/100.)**2)

            if Rmu<2.5 and Rpi<2.5:
                return True

    return False
    
############ Check if the decay is in the acceptance of the first element
def InTheSecondElement(pmu, ppi, Origin):

    #### I first check that the decay is inside the first vacuum tank
    if Origin[0]/100>110. and Origin[0]/100.<140:
        if TMath.sqrt((Origin[1]/100.)**2+(Origin[2]/100.)**2)<2.5:
            ### Now I check that the decay product go through the magnet
            ##### The end of the maget is located at 130+5 m (probably) 
            impact_mu = GetXY(14000., Origin, pmu)
            impact_pi = GetXY(14000., Origin, ppi)

            Rmu = TMath.sqrt((impact_mu["X"]/100.)**2+(impact_mu["Y"]/100.)**2)
            Rpi = TMath.sqrt((impact_pi["X"]/100.)**2+(impact_pi["Y"]/100.)**2)

            if Rmu<2.5 and Rpi<2.5:
                return True

    return False


##############################################################################################################

mass_mu = 0.1057
mass_pi = 0.1396

##### dump+muon filter 0-60m

########## The first element
##### vacuum tank 40m starts 60 ends at 100
##### first spectrometer 10m (100m, 110m)

########## The first element
##### vacuum tank 40m starts 110 ends at 150
##### second spectrometer 10m (140m, 160m)


firstCount = 0.
secondCount = 0.

for i in xrange(nentries-1):

    if (i%10000)==0:
        print "Event  ",i


    t.GetEntry(i)

    if t.SVz>-1000.:
        ##### Saving the momentum of the particles in the final state
        pPi = Get4Vector(name="Pi", t=t)
        pmu = Get4Vector(name="Mu", t=t)
        pmu0 = Get4Vector(name="Mu0", t=t)
        pN = Get4Vector(name="N", t=t)
        #pN = pmu+pPi
        pD = pN+pmu0
        ##### Getting the Origin of the particles
        OriginN = Get3Vector(name="SV", t=t, scale=1.)
        DecayN = Get3Vector(name="NDecay_", t=t, scale=1.)


        ###### now I count how many decay within the first element
        firstEl = InTheFirstElement(pmu, pPi, DecayN)
        if firstEl:
            firstCount+=1.

        ###### now I count how many decay within the first element
        secondEl = InTheSecondElement(pmu, pPi, DecayN)
        if secondEl:
            secondCount+=1.


t.GetEntry(nentries-1)
Norm = t.Norm

print "Acc First Element  %s"%(float(firstCount)/float(Norm))
print "based on %s events"%(firstCount)
print "Acc Second Element   %s"%(float(secondCount)/float(Norm))
print "based on %s events"%(secondCount)


