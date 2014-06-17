from __future__ import division
#import ROOT as r
import numpy as np
import gc
#from neutrinoDecayAndBoost import *
from saveHNLPdfs import *

def makeDecayInVolume(hh, decayString, volume, statistic=100):
    pN = r.TLorentzVector()
    vec = r.TVector3()
    inVol = 0
    if volume == 1:
        sourceHist = hh.weightedProdHistVol1
    elif volume == 2:
        sourceHist = hh.weightedProdHistVol2
    else:
        print 'makeDecayInVolume ERROR: please select volume 1 or 2'
        return 0
    for i in xrange(statistic):
        momN, thetaN = r.Double(), r.Double()
        sourceHist.GetRandom2(momN, thetaN)
        vec.SetMagThetaPhi(momN, thetaN, 0.)
        pN.SetVect(vec)
        pN.SetE(hh.pp.energy(momN, hh.pp.MN))
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        NDecayVtx = hh.ep.makeVtxInVolume(pN, hh.pp.c*hh.pp.NLifetime, volume)
        inVol += int( ep.inAcceptance(NDecayVtx, pGKid1, pGKid2, volume) )
    return inVol/statistic, statistic    





if __name__ == '__main__':
    pp = physicsParameters()
    pp.setNMass(1.)
    pp.setNCoupling([0.25e-08, 1.e-08, 0.5e-08])
    ep = experimentParams(pp, 'SHIP')
    hh = HistoHandler(pp, ep)
    rawPDF = hh.makeProductionPDF()
    accv1, accv2 = hh.scaleProductionPDF([0.25e-08, 1.e-08, 0.5e-08])
    print accv1, accv2
    fracV1, tot = makeDecayInVolume(hh, 'N -> pi mu', 1)
    print fracV1, tot
    hh.weightedPDFoutfile.Close()
    hh.prodPDFoutfile.Close()
    hh.charmFile.Close()









"""
def HNLAcceptance(u, mass, nMultiply = 10, nCharm = 100, process = 'Ds -> mu N, N -> mu pi', nEvents = 100):
    gc.collect()
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
            #print NProdVtx.Z(), pKid2.Pz(), NDecayVtx.Z(), ev.pp.NLifetime
            inVol1 += int( ep.inAcceptance(NDecayVtx, pGKid1, pGKid2, 1) )
            inVol2 += int( ep.inAcceptance(NDecayVtx, pGKid1, pGKid2, 2) )
    f.Close()
    return inVol1, inVol2, Norm


"""