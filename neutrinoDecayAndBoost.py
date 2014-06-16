#import ROOT as r
import numpy as np
import gc
from options import *

class myBranches():
    """ Contains all the silly lines of code to fill
    5 silly particles into a silly tree """
    def __init__(self):
        # Variables for the trees
        self.mum_Px = np.zeros(1, dtype=float)
        self.mum_Py = np.zeros(1, dtype=float)
        self.mum_Pz = np.zeros(1, dtype=float)
        self.mum_E = np.zeros(1, dtype=float)
        self.mum_P = np.zeros(1, dtype=float)
        self.mum_Theta = np.zeros(1, dtype=float)
        self.kid1_Px = np.zeros(1, dtype=float)
        self.kid1_Py = np.zeros(1, dtype=float)
        self.kid1_Pz = np.zeros(1, dtype=float)
        self.kid1_E = np.zeros(1, dtype=float)
        self.kid1_P = np.zeros(1, dtype=float)
        self.kid1_Theta = np.zeros(1, dtype=float)
        self.kid2_Px = np.zeros(1, dtype=float)
        self.kid2_Py = np.zeros(1, dtype=float)
        self.kid2_Pz = np.zeros(1, dtype=float)
        self.kid2_E = np.zeros(1, dtype=float)
        self.kid2_P = np.zeros(1, dtype=float)
        self.kid2_Theta = np.zeros(1, dtype=float)
        self.gKid1_Px = np.zeros(1, dtype=float)
        self.gKid1_Py = np.zeros(1, dtype=float)
        self.gKid1_Pz = np.zeros(1, dtype=float)
        self.gKid1_E = np.zeros(1, dtype=float)
        self.gKid1_P = np.zeros(1, dtype=float)
        self.gKid1_Theta = np.zeros(1, dtype=float)
        self.gKid2_Px = np.zeros(1, dtype=float)
        self.gKid2_Py = np.zeros(1, dtype=float)
        self.gKid2_Pz = np.zeros(1, dtype=float)
        self.gKid2_E = np.zeros(1, dtype=float)
        self.gKid2_P = np.zeros(1, dtype=float)
        self.gKid2_Theta = np.zeros(1, dtype=float)
    def fillNumbers(self, pMum, pKid1, pKid2, pGKid1, pGKid2):
        self.mum_Px[0]      = pMum.Px()
        self.mum_Py[0]      = pMum.Py()
        self.mum_Pz[0]      = pMum.Pz()
        self.mum_E[0]       = pMum.E()
        self.mum_P[0]       = pMum.P()
        self.mum_Theta[0]   = pMum.Theta()
        self.kid1_Px[0]     = pKid1.Px()
        self.kid1_Py[0]     = pKid1.Py()
        self.kid1_Pz[0]     = pKid1.Pz()
        self.kid1_E[0]      = pKid1.E()
        self.kid1_P[0]      = pKid1.P()
        self.kid1_Theta[0]  = pKid1.Theta()
        self.kid2_Px[0]     = pKid2.Px()
        self.kid2_Py[0]     = pKid2.Py()
        self.kid2_Pz[0]     = pKid2.Pz()
        self.kid2_E[0]      = pKid2.E()
        self.kid2_P[0]      = pKid2.P()
        self.kid2_Theta[0]  = pKid2.Theta()
        self.gKid1_Px[0]    = pGKid1.Px()
        self.gKid1_Py[0]    = pGKid1.Py()
        self.gKid1_Pz[0]    = pGKid1.Pz()
        self.gKid1_E[0]     = pGKid1.E()
        self.gKid1_P[0]     = pGKid1.P()
        self.gKid1_Theta[0] = pGKid1.Theta()
        self.gKid2_Px[0]    = pGKid2.Px()
        self.gKid2_Py[0]    = pGKid2.Py()
        self.gKid2_Pz[0]    = pGKid2.Pz()
        self.gKid2_E[0]     = pGKid2.E()
        self.gKid2_P[0]     = pGKid2.P()
        self.gKid2_Theta[0] = pGKid2.Theta()

    def branchMyTree(self, f, tree):
        # Self-explaining
        tree.Branch("mum_Px",      self.mum_Px,     "mum_Px/D")
        tree.Branch("mum_Py",      self.mum_Py,     "mum_Py/D")
        tree.Branch("mum_Pz",      self.mum_Pz,     "mum_Pz/D")
        tree.Branch("mum_E",       self.mum_E ,     "mum_E/D")
        tree.Branch("mum_P",       self.mum_P ,     "mum_P/D")
        tree.Branch("mum_Theta",   self.mum_Theta , "mum_Theta/D")
        tree.Branch("kid1_Px",     self.kid1_Px,     "kid1_Px/D")
        tree.Branch("kid1_Py",     self.kid1_Py,     "kid1_Py/D")
        tree.Branch("kid1_Pz",     self.kid1_Pz,     "kid1_Pz/D")
        tree.Branch("kid1_E",      self.kid1_E ,     "kid1_E/D")
        tree.Branch("kid1_P",      self.kid1_P ,     "kid1_P/D")
        tree.Branch("kid1_Theta",  self.kid1_Theta , "kid1_Theta/D")
        tree.Branch("kid2_Px",     self.kid2_Px,     "kid2_Px/D")
        tree.Branch("kid2_Py",     self.kid2_Py,     "kid2_Py/D")
        tree.Branch("kid2_Pz",     self.kid2_Pz,     "kid2_Pz/D")
        tree.Branch("kid2_E",      self.kid2_E ,     "kid2_E/D")
        tree.Branch("kid2_P",      self.kid2_P ,     "kid2_P/D")
        tree.Branch("kid2_Theta",  self.kid2_Theta , "kid2_Theta/D")
        tree.Branch("gKid1_Px",    self.gKid1_Px,     "gKid1_Px/D")
        tree.Branch("gKid1_Py",    self.gKid1_Py,     "gKid1_Py/D")
        tree.Branch("gKid1_Pz",    self.gKid1_Pz,     "gKid1_Pz/D")
        tree.Branch("gKid1_E",     self.gKid1_E ,     "gKid1_E/D")
        tree.Branch("gKid1_P",     self.gKid1_P ,     "gKid1_P/D")
        tree.Branch("gKid1_Theta", self.gKid1_Theta , "gKid1_Theta/D")
        tree.Branch("gKid2_Px",    self.gKid2_Px,     "gKid2_Px/D")
        tree.Branch("gKid2_Py",    self.gKid2_Py,     "gKid2_Py/D")
        tree.Branch("gKid2_Pz",    self.gKid2_Pz,     "gKid2_Pz/D")
        tree.Branch("gKid2_E",     self.gKid2_E ,     "gKid2_E/D")
        tree.Branch("gKid2_P",     self.gKid2_P ,     "gKid2_P/D")
        tree.Branch("gKid2_Theta", self.gKid2_Theta , "gKid2_Theta/D")
        return f, tree


def makeRFNtuple(u, mass, process = 'Ds -> mu N, N -> mu pi', nEvents = 100):
    """ Let a charmed meson decay into sterile neutrino in its rest frame,
    then let the sterile neutrino decay in the charmed meson rest frame. """
    # Setup physics
    pp = physicsParameters()
    pp.setNMass(mass)
    pp.setNCoupling(u)

    # Set process
    ev = myNEvent(pp, process)
    treeName = "%s_%s" %(''.join(ev.production.particles), ''.join(ev.decay.particles))

    # Prepare the tree
    #rootFileName = "out/NTuples/restFrame_U%s_m%s.root"%(u,mass)
    f = r.TFile(ev.rootFileName,"update")
    # Check if this process was already simulated
    fileContent = [key.GetName() for key in r.gDirectory.GetListOfKeys()]
    if treeName in fileContent:
        print "makeNtuple(%s, %s, %s, %s): TTree %s already present in %s, SKIPPING!"%(u,mass,process,nEvents,treeName,ev.rootFileName)
        return ev, f, f.Get(treeName)
    tree = r.TTree(treeName,treeName)
    # Branch the tree
    b = myBranches()
    b.branchMyTree(f, tree)

    # Populate the tree
    for i in xrange(nEvents):
        if not i%10000:
            print "Making rest frame decay %s..."%i
            gc.collect()

        pKid1, pKid2 = ev.production.makeDecay()
        ev.decay.setPMother(pKid2)
        pGKid1, pGKid2 = ev.decay.makeDecay()

        b.fillNumbers(ev.production.pMother, pKid1, pKid2, pGKid1, pGKid2)
        tree.Fill()

    tree.Write()
    gc.collect()
    #f.Close()
    #return ev, b, f, tree
    return ev, f, tree

def boostEvents(u, mass, nCharm = 0, process = 'Ds -> mu N, N -> mu pi', nEvents = 100):
    """ Read charmed mesons four-vectors from a file,
    and boost the products of the selected decay chain """
    # Take rest frame decays
    #ev, b, f, RFtree = makeRFNtuple(u, mass, process, nEvents)
    ev, f, RFtree = makeRFNtuple(u, mass, process, nEvents)
    btName = "boosted_%s_%s" %(''.join(ev.production.particles), ''.join(ev.decay.particles))
    fileContent = [key.GetName() for key in r.gDirectory.GetListOfKeys()]
    if btName in fileContent:
        #rootFileName = "out/NTuples/restFrame_U%s_m%s.root"%(u,mass)
        print "boostEvents(%s, %s, %s, %s, %s): TTree %s already present in %s, SKIPPING!"%(u,mass,nCharm,process,nEvents,btName,ev.rootFileName)
        return ev#, f, f.Get(btName)
    nRF = RFtree.GetEntriesFast()
    # Take charmed mesons
    charmFile = r.TFile(ev.pp.charmSourceFile, "read")
    charmTree = charmFile.Get(ev.pp.charmTreeName)
    if not nCharm:
        nCharm = charmTree.GetEntriesFast()
    print "Preparing to store %s events!"%(nRF*nCharm)
    # Crate a new tree with boosted particles
    bTree = r.TTree(btName, btName)
    # Branch the tree
    bb = myBranches()
    bb.branchMyTree(f, bTree)
    pVtx_X = np.zeros(1, dtype=float)
    pVtx_Y = np.zeros(1, dtype=float)
    pVtx_Z = np.zeros(1, dtype=float)
    pVtx_T = np.zeros(1, dtype=float)
    bTree.Branch("prodVtx_X",    pVtx_X,     "prodVtx_X/D")
    bTree.Branch("prodVtx_Y",    pVtx_Y,     "prodVtx_Y/D")
    bTree.Branch("prodVtx_Z",    pVtx_Z,     "prodVtx_Z/D")
    bTree.Branch("prodVtx_T",    pVtx_T,     "prodVtx_T/D")    

    # Populate the tree
    nBoosted = 0
    for i in xrange(nCharm):#nCharm
        charmTree.GetEntry(i)
        pCharmed = r.TLorentzVector(charmTree.CharmPx, charmTree.CharmPy, charmTree.CharmPz, charmTree.CharmE)
        prodVtx = r.TLorentzVector(charmTree.SVx, charmTree.SVy, charmTree.SVz, charmTree.DecayTime)
        pBoost = pCharmed.BoostVector()
        for j in xrange(nRF):#nrf
            nBoosted += 1
            if not nBoosted%10000:
                print "Boosting decay chain %s..."%nBoosted
                gc.collect()
            RFtree.GetEntry(j)
            pKid1 = r.TLorentzVector(RFtree.kid1_Px, RFtree.kid1_Py, RFtree.kid1_Pz, RFtree.kid1_E)
            pKid2 = r.TLorentzVector(RFtree.kid2_Px, RFtree.kid2_Py, RFtree.kid2_Pz, RFtree.kid2_E)
            pGKid1 = r.TLorentzVector(RFtree.gKid1_Px, RFtree.gKid1_Py, RFtree.gKid1_Pz, RFtree.gKid1_E)
            pGKid2 = r.TLorentzVector(RFtree.gKid2_Px, RFtree.gKid2_Py, RFtree.gKid2_Pz, RFtree.gKid2_E)
            pKid1.Boost(pBoost)
            pKid2.Boost(pBoost)
            pGKid1.Boost(pBoost)
            pGKid2.Boost(pBoost)

            bb.fillNumbers(pCharmed, pKid1, pKid2, pGKid1, pGKid2)
            pVtx_X[0] = prodVtx.X()
            pVtx_Y[0] = prodVtx.Y()
            pVtx_Z[0] = prodVtx.Z()
            pVtx_T[0] = prodVtx.T()
            bTree.SetDirectory(0)
            bTree.Fill()

    f.cd()
    bTree.Write()
    gc.collect()
    f.Close()
    charmFile.Close()
    return ev#, f, bTree