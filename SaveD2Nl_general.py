
from ROOT import * 
from array import array

from options import *
##########################################################################



###### Add a branch variable to the tree
def __AddVar__(self, var, name, res):
    has_branche = res.has_key(name)
            
    if has_branche:
        res[name][0] = float(var)
    else:
        comp = array('f', [0])
        self.Branch(name, comp, name+'/F')
        comp[0]= float(var)
        aux = {name: comp}
        res.update(aux)
    return res

TTree.AddVar = __AddVar__

######## Copy all the branches of the tree
def CopyTree(tree, newtree, res_tuple):

    branches = tree.GetListOfBranches()
    var_list = []
    for b in branches:
        var_list.append(b.GetName())

    for n in var_list:
        res_tuple = newtree.AddVar(tree.__getattr__(n), n, res_tuple)

##########################################################################


#######################################################################################################################
def GetNDecayPoint(Origin, Momentum, lifetime):

    c = 3*1e10
    v = Momentum.Beta()*c  
    DT =gRandom.Exp(lifetime)*Momentum.Gamma()
    DL = DT*v

    Direction = Momentum.Vect().Unit()
    EndVertex = TVector3(Momentum.Vect().Unit()[0]*DL, Momentum.Vect().Unit()[1]*DL, Momentum.Vect().Unit()[2]*DL) 
    EndVertex = Origin+EndVertex

    return EndVertex


def Get3Vector(name, t):
    X = t.__getattr__(name+"x")/10.
    Y = t.__getattr__(name+"y")/10.
    Z = t.__getattr__(name+"z")/10.    
    vect = TVector3(Z,X,Y)

    return vect

    
def Get4Vector(name, t):
    X = t.__getattr__(name+"_Px")
    Y = t.__getattr__(name+"_Py")
    Z = t.__getattr__(name+"_Pz")
    E = t.__getattr__(name+"_E")    
    lorentz = TLorentzVector(Z,X,Y,E)
    
    return lorentz
#######################################################################################################################

geoManager = TGeoManager.Import("CNGS_geometry.root")

chain = TChain("newTree", "newTree")
finput = "%s_M%s.root"%(outName,m1/1000.)
chain.AddFile(finput)

#chain.AddFile("D2Nmu_445_301_600.root")
#chain.AddFile("D2Nmu_445_601_900.root")
######## 


#D2Nmu_445_601_900.root


write2file = True


nentries = chain.GetEntries()

####U2 = 6e-7: 3e-06, 'BR': 4.5e-9
#lifetime =  3e-6
####U2 = 1e-9: 0.0018, 'BR': 7.5e-12
#lifetime =  0.0018


###### the lifetime is between 3e-6 (for U2=6e-7) and 0.002 (for U2=1e-9)
#lifetime = 3e-5

###### Create a sample with Hans Lifetime
#lifetime = 1.3e-5


count = 0.
Norm = 0.
Ncount = 0.
geoDet = 0.
Multiply = 20

tree_acceptance = TTree("tree_acceptance", "tree_acceptance")
res = {}






for i in xrange(nentries):


    if (i%10000)==0:
        print "Event number ",i

    chain.GetEntry(i)
    ### I define the 3D Point that I will use
    OriginN = Get3Vector(name="SV", t=chain)
    #print OriginN[0], OriginN[1], OriginN[2]
    pKid1 = Get4Vector(name=kid1, t=chain)
    pGKid1 = Get4Vector(name=gkid1, t=chain)
    pGKid2 = Get4Vector(name=gkid2, t=chain)
    pKid2 = Get4Vector(name=kid2, t=chain)



    for j in xrange(Multiply):

        Norm+=1.
        EndVertex = GetNDecayPoint(OriginN, pKid1, lifetime)


        
    ####Look if the decay of the neutrino is in the tunnel 
        if EndVertex[0]/100.<200.:

            if write2file:
                CopyTree(chain, tree_acceptance, res)
                tree_acceptance.AddVar(EndVertex[0], "NDecay_z", res)
                tree_acceptance.AddVar(EndVertex[1], "NDecay_x", res)
                tree_acceptance.AddVar(EndVertex[2], "NDecay_y", res)
                #tree_acceptance.AddVar(Norm, "Norm", res)
                tree_acceptance.SetDirectory(0)
                tree_acceptance.Fill()

tree_acceptance.AddVar(Norm, "Norm", res)
tree_acceptance.Fill()

print "prob %s"%(count/float(Norm))
print "acc %s"%(geoDet/float(nentries))

if write2file:
    filename = finput.replace(".root", "_LT%s.root"%(lifetime))
    #filename = "NeutrinosLHC_acceptance_LT_%s.root"%(lifetime)
    facc = TFile(filename, "RECREATE")
    tree_acceptance.Write()
    facc.Save()
    facc.Close()







