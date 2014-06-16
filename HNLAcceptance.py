from __future__ import division
#import ROOT as r
import numpy as np
import gc
from neutrinoDecayAndBoost import *

def HNLAcceptance(u, mass, nMultiply = 10, nCharm = 100, process = 'Ds -> mu N, N -> mu pi', nEvents = 100):
    #ev, f, bTree = boostEvents(u, mass, nCharm, process, nEvents)
    ev = boostEvents(u, mass, nCharm, process, nEvents)
    ev.pp.computeNLifetime()
    ep = experimentParams(ev.pp, 'SHiP')
    f = r.TFile(ev.rootFileName,"read")
    btName = "boosted_%s_%s" %(''.join(ev.production.particles), ''.join(ev.decay.particles))
    bTree = f.Get(btName)
    nBoostedHNLs = bTree.GetEntries()#Fast()
    Norm = 0.
    inVol1 = 0.
    inVol2 = 0.
    for i in xrange(nBoostedHNLs):
        bTree.GetEntry(i)
        pKid2 = r.TLorentzVector(bTree.kid2_Px, bTree.kid2_Py, bTree.kid2_Pz, bTree.kid2_E) # pN
        pGKid1 = r.TLorentzVector(bTree.gKid1_Px, bTree.gKid1_Py, bTree.gKid1_Pz, bTree.gKid1_E)
        pGKid2 = r.TLorentzVector(bTree.gKid2_Px, bTree.gKid2_Py, bTree.gKid2_Pz, bTree.gKid2_E)
        NProdVtx = r.TVector3(bTree.prodVtx_X, bTree.prodVtx_Y, bTree.prodVtx_Z)
        for j in xrange(nMultiply):
            Norm = Norm + 1.
            if not Norm%100000:
                print "Checking acceptance for event %s of %s..."%(int(Norm), nBoostedHNLs*nMultiply)
                gc.collect()
            NDecayVtx = ev.pp.drawNDecayVtx(pKid2, NProdVtx)
            inVol1 += int( ep.inAcceptance(NDecayVtx, pGKid1, pGKid2, 1) )
            inVol2 += int( ep.inAcceptance(NDecayVtx, pGKid1, pGKid2, 2) )
    f.Close()
    return inVol1, inVol2, Norm


