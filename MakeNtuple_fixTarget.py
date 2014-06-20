
from ROOT import * 
from array import array

# where the decay options are set
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


###### Opening the ntuple of production
massN = m1 #1.6
fname = "CharmFixTarget.root"
#fname = "CharmLHC_CM.root"
fprod = TFile(fname)
tprod = fprod.Get("newTree;1")
#####  Opening the ntuple for decay 
fname = "%s_CM_M%.1f.root"%(outName,m1/1000.) 

fdec = TFile("%s"%(fname))
tdec = fdec.Get("tuple;1")


##### Creating the new ntuple

tree = TTree("newTree", "newTree")
res = {}


nentries = tprod.GetEntries()

start = 1
end = 100 #tdec.GetEntries()



for j in xrange(start, end):  
    if (j%10)==0:
        print "Doing ",j
    tdec.GetEntry(j)
    for i in xrange(nentries):

        #if (i%1000)==0:
        #    print "Processing event  %s "%(i)
        tprod.GetEntry(i)
        ##### if I am outside the range I restart to loop from the first event

    
        pMother =  TLorentzVector(tprod.CharmPx, tprod.CharmPy, tprod.CharmPz, tprod.CharmE)
        BoostD = pMother.BoostVector()
        SV =  TLorentzVector(tprod.SVx, tprod.SVy, tprod.SVz, tprod.DecayTime)
        ## kid1
        pKid1 = TLorentzVector(tdec.__getattr__("%s_Px"%kid1)/1000., tdec.__getattr__("%s_Py"%kid1)/1000., tdec.__getattr__("%s_Pz"%kid1)/1000., tdec.__getattr__("%s_E"%kid1)/1000.)
        pKid1.Boost(BoostD)
        ## gkid1
        pGKid1 = TLorentzVector(tdec.__getattr__("%s_Px"%gkid1)/1000., tdec.__getattr__("%s_Py"%gkid1)/1000., tdec.__getattr__("%s_Pz"%gkid1)/1000., tdec.__getattr__("%s_E"%gkid1)/1000.)
        pGKid1.Boost(BoostD)
        ## gkid2
        pGKid2 = TLorentzVector(tdec.__getattr__("%s_Px"%gkid2)/1000., tdec.__getattr__("%s_Py"%gkid2)/1000., tdec.__getattr__("%s_Pz"%gkid2)/1000., tdec.__getattr__("%s_E"%gkid2)/1000.)
        pGKid2.Boost(BoostD)
        ## kid2
        pKid2 = TLorentzVector(tdec.__getattr__("%s_Px"%kid2)/1000., tdec.__getattr__("%s_Py"%kid2)/1000., tdec.__getattr__("%s_Pz"%kid2)/1000., tdec.__getattr__("%s_E"%kid2)/1000.)
        pKid2.Boost(BoostD)
        
    ##### Getting The decay vertex
        #if tprod.CharmPID==431:# or tprod.CharmPID==421 or tprod.CharmPID==411:
        if tprod.CharmPID==431 or tprod.CharmPID==411:
            tree.AddVar(pMother.X(), "%s_Px"%mother, res)
            tree.AddVar(pMother.Y(), "%s_Py"%mother, res)
            tree.AddVar(pMother.Z(), "%s_Pz"%mother, res)
            tree.AddVar(pMother.E(), "%s_E"%mother, res)
            
            tree.AddVar(pKid1.X(), "%s_Px"%kid1, res)
            tree.AddVar(pKid1.Y(), "%s_Py"%kid1, res)
            tree.AddVar(pKid1.Z(), "%s_Pz"%kid1, res)
            tree.AddVar(pKid1.E(), "%s_E"%kid1, res)    

            tree.AddVar(pGKid1.X(), "%s_Px"%gkid1, res)
            tree.AddVar(pGKid1.Y(), "%s_Py"%gkid1, res)
            tree.AddVar(pGKid1.Z(), "%s_Pz"%gkid1, res)
            tree.AddVar(pGKid1.E(), "%s_E"%gkid1, res)
    
            tree.AddVar(pGKid2.X(), "%s_Px"%gkid2, res)
            tree.AddVar(pGKid2.Y(), "%s_Py"%gkid2, res)
            tree.AddVar(pGKid2.Z(), "%s_Pz"%gkid2, res)
            tree.AddVar(pGKid2.E(), "%s_E"%gkid2, res)

            tree.AddVar(pKid2.X(), "%s_Px"%kid2, res)
            tree.AddVar(pKid2.Y(), "%s_Py"%kid2, res)
            tree.AddVar(pKid2.Z(), "%s_Pz"%kid2, res)
            tree.AddVar(pKid2.E(), "%s_E"%kid2, res)
    
            tree.AddVar(SV.X(), "SVx", res)
            tree.AddVar(SV.Y(), "SVy", res)
            tree.AddVar(SV.Z(), "SVz", res)
            tree.AddVar(SV.T(), "SVt", res)
            
            tree.SetDirectory(0)
            tree.Fill()

###########################################################################
###########################################################################

#foutput = "D2Nmu_LHC_%s_%s_%s.root"%(n, start, end)
foutput = "%s_M%s.root"%(outName,m1/1000.)
fnew = TFile(foutput, "RECREATE")
tree.Print()
tree.Write()
fnew.Save()
fnew.Close()



