from __future__ import division
from options import *

class HistoHandler():
    """ Manage PDFs for sterile neutrino production and decay """
    def __init__(self, pp, ep, process = 'Ds -> mu N, N -> mu pi', binsp = 90, binstheta = 80):
        # Setup input
        self.pp = pp
        self.ep = ep
        self.charmFile = r.TFile(self.pp.charmSourceFile, 'read')
        self.charmTree = self.charmFile.Get(self.pp.charmTreeName)
        # Setup physics
        self.ev = myNEvent(self.pp, process)
        self.binstheta = binstheta
        self.binsp = binsp
    def makeProductionPDF(self):
        """ To run once for every different HNL mass """
        # Make mass histogram
        self.prodHistName = "%s"%(''.join(self.ev.production.particles))
        self.tmax = self.ep.v1ThetaMax
        self.tmin = 0.
        self.thetaStep = (self.tmax-self.tmin)/self.binstheta
        self.pMax = math.sqrt(self.ep.protonEnergy**2. - self.pp.MN**2.)
        self.pMin = self.pp.MN
        self.pStep = (self.pMax-self.pMin)/self.binsp
        self.prodHist = r.TH2F("prodPDF_m%s"%(pp.MN) ,"prodPDF_m%s"%(pp.MN),
            self.binsp,self.pMin-0.5*self.pStep,self.pMax-0.5*self.pStep,
            self.binstheta,self.tmin-0.5*self.thetaStep,self.tmax+0.5*self.thetaStep)
        self.prodHist.SetTitle("PDF for N production (m_{N}=%s GeV)"%(pp.MN))
        self.prodHist.GetXaxis().SetTitle("P_{N} [GeV]")
        self.prodHist.GetYaxis().SetTitle("#theta_{N} [rad]")
        # Fill mass histogram
        for charm in self.charmTree:
            pCharm = r.TLorentzVector(charm.CharmPx, charm.CharmPy, charm.CharmPz, charm.CharmE)
            self.ev.production.setPMother(pCharm)
            pKid1, pKid2 = self.ev.production.makeDecay()
            self.prodHist.Fill(pKid2.P(), pKid2.Theta())
        # Make it a PDF
        histInt = self.prodHist.Integral("width")
        self.prodHist.Scale(1./histInt)
        self.prodPDFfilename = 'out/NTuples/prodPDF_m%s.root'%(self.pp.MN)
        self.prodPDFoutfile = r.TFile(self.prodPDFfilename,'recreate')
        self.prodHist.Write("",5)
        #self.prodPDFoutfile.Close()
        #self.charmFile.Close()
    def scaleProductionPDF(self, couplings):
        # Now take the couplings and scale the PDF
        self.pp.setNCoupling(couplings)
        self.pp.computeNLifetime()
        ct = self.pp.c*self.pp.NLifetime
        # Prepare the weighted histograms
        self.couplingString = '_'.join([str(ui) for ui in self.pp.U2])
        self.weightedProdHistVol1 = r.TH2F("weightedPDF_vol1_m%s_couplings%s"%(pp.MN, self.couplingString),
            "weightedPDF_vol1_m%s_couplings%s"%(pp.MN, self.couplingString),
            self.binsp,self.pMin-0.5*self.pStep,self.pMax-0.5*self.pStep,
            self.binstheta,self.tmin-0.5*self.thetaStep,self.tmax+0.5*self.thetaStep)
        self.weightedProdHistVol2 = r.TH2F("weightedPDF_vol2_m%s_couplings%s"%(pp.MN, self.couplingString),
            "weightedPDF_vol2_m%s_couplings%s"%(pp.MN, self.couplingString),
            self.binsp,self.pMin-0.5*self.pStep,self.pMax-0.5*self.pStep,
            self.binstheta,self.tmin-0.5*self.thetaStep,self.tmax+0.5*self.thetaStep)
        fourMom = r.TLorentzVector()
        vec = r.TVector3()
        index = 0
        self.accVol1 = 0.
        self.accVol2 = 0.
        # Fill the weighted histogram
        for p in xrange(self.binsp):
            for th in xrange(self.binstheta):
                index += 1
                weight = self.prodHist.GetBinContent(p, th)
                if weight > 0.:
                    mom = self.prodHist.GetXaxis().GetBinCenter(p)
                    angle = self.prodHist.GetYaxis().GetBinCenter(th)
                    binWeight = weight*self.pStep*self.thetaStep
                    vec.SetMagThetaPhi(mom, angle, 0.)
                    fourMom.SetE(pp.energy(mom, pp.MN))
                    fourMom.SetVect(vec)
                    probVtx1 = ep.probVtxInVolume(fourMom, ct, 1)
                    probVtx2 = ep.probVtxInVolume(fourMom, ct, 2)
                    #px = fourMom.Px()
                    #pz = fourMom.Pz()
                    accGeo1 = ep.GeometricAcceptance(fourMom, 1)
                    accGeo2 = ep.GeometricAcceptance(fourMom, 2)
                    acc1 = (binWeight * accGeo1 * probVtx1)
                    acc2 = (binWeight * accGeo2 * probVtx2)
                    self.accVol1 += acc1
                    self.accVol2 += acc2
                    self.weightedProdHistVol1.Fill(mom,angle,acc1)
                    self.weightedProdHistVol2.Fill(mom,angle,acc2)
        # Save the PDF
        self.outFileName = 'out/NTuples/m%s_couplings%s.root'%(self.pp.MN, self.couplingString)
        outfile = r.TFile(self.outFileName,'update')
        self.weightedProdHistVol1.Write("",5)
        self.weightedProdHistVol2.Write("",5)
        outfile.Close()
        self.prodPDFoutfile.Close()
        self.charmFile.Close()
        return self.accVol1, self.accVol2


if __name__ == '__main__':
    pp = physicsParameters()
    pp.setNMass(1.)
    pp.setNCoupling([0.25e-08, 1.e-08, 0.5e-08])
    ep = experimentParams(pp, 'SHIP')
    hh = HistoHandler(pp, ep)
    rawPDF = hh.makeProductionPDF()
    accv1, accv2 = hh.scaleProductionPDF([0.25e-08, 1.e-08, 0.5e-08])
    print accv1, accv2
