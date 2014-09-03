from __future__ import division
import numpy as np
import gc
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
        kids = HNLDecayChain(hh, decayString, pN)
        NDecayVtx = hh.ep.makeVtxInVolume(pN, hh.pp.c*hh.pp.NLifetime, volume)
        inVol += int( hh.ep.inAcceptance(NDecayVtx, kids, volume) )
    return inVol/statistic, statistic    

def HNLDecayChain(hh, decayString, pN):
    """ Returns the 4-momenta of the final detectable daughters """
    kids = []
    if decayString == 'N -> nu nu nu':
        kids = []
    elif decayString in hh.pp.decays[:7]: # 3-body, return the first two
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2, pGKid3 = hh.ev.decay.makeDecay()
        kids = [pGKid1, pGKid2]
    elif decayString == 'N -> rho nu':
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        hh.ev.decay.readString('rho -> pi pi')
        hh.ev.decay.setPMother(pGKid1)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        pPi1, pPi2 = r.TLorentzVector(pGKid1), r.TLorentzVector(pGKid2)
        kids = [pPi1, pPi2]
    elif decayString in hh.pp.decays[11:13]: # N -> rho l
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay() # rho+lepton
        pLepton = r.TLorentzVector(pGKid2)
        hh.ev.decay.readString('rho -> pi pi0')
        hh.ev.decay.setPMother(pGKid1)
        pGKid1, pGKid3 = hh.ev.decay.makeDecay() # pi+pi0
        pPi = r.TLorentzVector(pGKid1)
        hh.ev.decay.readString('pi0 -> gamma gamma')
        hh.ev.decay.setPMother(pGKid3)
        pGKid3, pGKid4 = hh.ev.decay.makeDecay() # 2 gamma
        pG1, pG2 = r.TLorentzVector(pGKid3), r.TLorentzVector(pGKid4)
        kids = [pLepton, pPi, pG1, pG2] # lepton+pi+2gamma
    elif decayString == 'N -> pi0 nu': # DIFFICILE DA MISURARE, RICHIEDE STUDIO BG
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        hh.ev.decay.readString('pi0 -> gamma gamma')
        hh.ev.decay.setPMother(pGKid1)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        pG1, pG2 = r.TLorentzVector(pGKid1), r.TLorentzVector(pGKid2)
        kids = [pG1, pG2]
    else: #two-body, charged
        hh.ev.decay.readString(decayString)
        hh.ev.decay.setPMother(pN)
        pGKid1, pGKid2 = hh.ev.decay.makeDecay()
        kids = [pGKid1, pGKid2]   
    return kids


def computeNEvents(model, mass, coupling, root_dir_path='/home/elena/Desktop/PhD-Work/0-SHIP/HiddenParticlesSensitivity/SterileNeutrinos'):
    """ Choose model 1, 2, 3, 4 or 5 """
    pp = physicsParameters(root_dir_path)
    pp.setNMass(mass)                   
    # Check kinematics
    sources = []
    if pp.MN < max(pp.m_max['charm']):
        sources.append('charm')
    if pp.MN < max(pp.m_max['beauty']):
        sources.append('beauty')
    else:
        return 0.
    pp.setNCoupling([coupling*f for f in pp.models[model-1]])
    ep = experimentParams(pp, 'SHIP')

    BRprod = []
    NEv = 0.
    for source in sources:
        hh = HistoHandler(pp, ep, source, model)
        if hh.makeProductionPDF():
            accv1, accv2 = hh.scaleProductionPDF(pp.U2)
            decList = hh.pp.HNLAllowedDecays()
            weight1 = 0.
            weight2 = 0.
            for dec in decList:
                if decList[dec] == 'yes' and dec != 'N -> pi0 nu' and dec != 'N -> nu nu nu':
                    if accv1:
                        acc1, tot1 = makeDecayInVolume(hh, dec, 1)
                    else:
                        acc1 = 0.
                    if accv2:
                        acc2, tot2 = makeDecayInVolume(hh, dec, 2)
                    else:
                        acc2 = 0.
                    weight1 += hh.pp.findBranchingRatio(dec)*acc1
                    weight2 += hh.pp.findBranchingRatio(dec)*acc2
            BRprod_s = 0.
            U2tot = sum(pp.U2)
            for flavour in [1,2,3]:
                BRprod_s += productionBR(pp, source, flavour) * pp.U2[flavour-1]/U2tot
            NEv += (accv1*weight1 + accv2*weight2)*2.*BRprod_s*hh.ep.protonFlux
            BRprod.append(BRprod_s)
        hh.sourceFile.Close()

    outFilePath = root_dir_path + "/out/TextData/sensitivityScan-HNL-model%s.txt"%(model)
    with open(outFilePath,"a") as ofile:
        try:
            ofile.write('%s \t %s \t %s \t %s\n'%(mass, coupling, sum(BRprod), NEv))
        except KeyboardInterrupt:
            pass
    return NEv





"""
if __name__ == '__main__':
    pp = physicsParameters()
    pp.setNMass(1.)
    pp.setNCoupling([0.25e-08, 1.e-08, 0.5e-08])
    ep = experimentParams(pp, 'SHIP')
    hh = HistoHandler(pp, ep)
    rawPDF = hh.makeProductionPDF()
    accv1, accv2 = hh.scaleProductionPDF([0.25e-08, 1.e-08, 0.5e-08])
    print accv1, accv2
    decList = HNLAllowedDecays(hh.pp)
    #print decList
    tot = 0.
    weight1 = 0.
    weight2 = 0.
    for dec in decList:
        if decList[dec] == 'yes' and dec != 'N -> pi0 nu' and dec != 'N -> nu nu nu':
            acc1, tot1 = makeDecayInVolume(hh, dec, 1)
            acc2, tot2 = makeDecayInVolume(hh, dec, 2)
            print dec + '\t', acc1, acc2, '\t BR: ', hh.pp.findBranchingRatio(dec)
            tot += hh.pp.findBranchingRatio(dec)
            weight1 += hh.pp.findBranchingRatio(dec)*acc1
            weight2 += hh.pp.findBranchingRatio(dec)*acc2
    #fracV1, tot = makeDecayInVolume(hh, 'N -> mu mu nu', 1)
    print tot, weight1, weight2
    NEv = (accv1*weight1 + accv2*weight2)*2.*hh.pp.Xcc*hh.pp.computeNProdBR(2)*hh.ep.protonFlux
    #print fracV1, tot
    print 'NEv: ', NEv
    hh.weightedPDFoutfile.Close()
    hh.prodPDFoutfile.Close()
    hh.charmFile.Close()

"""







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