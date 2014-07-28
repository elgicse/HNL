

#from ROOT import *
import ROOT as r
from array import array

#######################################################
###### Here I compute the form factors
def fplus(q2):
    val = ((2.007)**2/((2.007)**2-q2))*0.726
    return val

def fzero(q2):
    val = ((1.870)**2/((1.870)**2-q2))*0.726
    return val

def fminus(q2, mH, mh):
    val = (fzero(q2)-fplus(q2))*(mH**2-mh**2)/q2
    return val
#######################################################

#### Function tat
def PhaseSpace3Body(q2, EN, mH, mh, mN, ml):

    fp = fplus(q2)
    fm = fminus(q2, mH, mh)
    PhaseSpace = (fm**2)*(q2*(mN**2+ml**2)-(mN**2-ml**2)**2)+2.*fp*fm*((mN**2)*(2*mH**2-2*mh**2-4*EN*mH-ml**2+mN**2+q2)+(ml**2)*(4.*EN*mH+ml**2-mN**2-q2))
    PhaseSpace += (fp**2)*((4.*EN*mH+ml**2-mN**2-q2)*(2.*mH**2-2.*mh**2-4.*EN*mH-ml**2+mN**2+q2)-(2.*mH**2+2.*mh**2-q2)*(q2-mN**2-ml**2))

    return PhaseSpace



def Integrate3Body(mH, mh, mN, ml, nToys=1000):

    #### Setting up the parameters for generation
    W = r.TLorentzVector(0., 0., 0., mH)
    masses = array('d', [mh, mN, ml])
    event = r.TGenPhaseSpace()
    event.SetDecay(W, 3, masses)

    
    Nq2 = 20 
    NEN = 20
    ENMax = (mH-mh)
    ENMax =  r.TMath.Sqrt(((ENMax**2-ml**2-mN**2)/2.)**2+mN**2)
    hist = r.TH2F("hist", "", Nq2, (ml+mN)**2, (mH-mh)**2, NEN, mN, ENMax)

    ###### Integral
    Integral = 0.

    #### For loop in order to integrate
    for i in xrange(nToys):
        event.Generate()
        #### Getting momentum of the daughters
        ph = event.GetDecay(0)
        pN    = event.GetDecay(1)
        pl    = event.GetDecay(2)
        #### Computing the parameters to compute the phase space
        q = pl+pN
        q2 = q.M2()
        EN = pN.E()

        iBin = hist.Fill(q2, EN)
        
        if hist.GetBinContent(iBin)==1.:
            val = PhaseSpace3Body(q2, EN, mH, mh, mN, ml)
            Integral+=val*hist.GetXaxis().GetBinWidth(1)*hist.GetYaxis().GetBinWidth(1)

    hist.Delete()
    return {'Int':Integral, "PHSP":hist}






def ComputeBR3BodyApp(U2, mass, tauH, mH, Vckm, fH, ml):

    GF2 = 0.608498183794
    const = (GF2*tauH*U2*Vckm**2)
    x = mass/mH
    y = ml/mH
    Diagram = (tauH*(GF2)*(fH**2)*mH)/(8.*r.TMath.Pi())*(Vckm**2)*(U2)
    PhaseSpace = (1-x**2+2*y**2+(y/x)*(1-y**2))*r.TMath.Sqrt((1+x**2-y**2)**2-4.*x**2)*mass**2
    #PhaseSpace = ml**2*(1.-(ml/mH)**2)**2

    #print "********************************************"
    print "Phase Space  ",PhaseSpace
    #print "********************************************"

    BR = Diagram*PhaseSpace

    return BR


def ComputeBR3Body(U2, mass, tauH, mH, mh, Vckm, ml):
    GF2 = 0.608498183794
    const = (GF2*tauH*U2*Vckm**2)
    Diagram = const/(64.*r.TMath.Pi()**3*mH**2)
    a = Integrate3Body(mH, mh, mass, ml, nToys=5000)
    PhaseSpace = a['Int']
    BR=Diagram*PhaseSpace
    print "Diagram ", Diagram
    print "PHSP ", PhaseSpace

    return BR



##### Parameters
ml = 0.1057
mDs = 1.9685
mD = 1.8696
mD0 = 1.865
fDs = 0.2801
fD = 0.2226

#### ctau for the Ds
tauDs = 149.9
tauD = 311.8
tauD0 = 122.9
mass = 0.5#1.
U2= 6e-7
Vcs = 0.97345
Vcd = 0.2252
Vcb = 1.0
mK = 0.497
mKst = 0.892
mpi = 0.140
mB = 5270.5
tauB = 492.



def phsp2body(mH, mN, ml):
    p1 = 1. - mN**2./mH**2. + 2.*(ml**2./mH**2.) + (ml**2./mH**2.)*(1. - (ml**2./mH**2.))
    p2 = ( 1. + (mN**2./mH**2.) - (ml**2./mH**2.) )**2. - 4.*(mN**2./mH**2.)
    p3 = r.TMath.Sqrt( p2 )
    return p1*p3


BR = ComputeBR3Body(U2, mass, tauD, mD, mK, Vcs, ml)
print "The BR(D) is ",BR

BR = ComputeBR3Body(U2, mass, tauD0, mD0, mK, Vcs, ml)
print "The BR(D0) is ",BR


BR = ComputeBR3Body(U2, mass, tauD0, mD0, mpi, Vcd, ml)
print "The BR(D0) is ",BR

#mass = 2.0
#U2 = 3e-7
#BR = ComputeBR3Body(U2, mass, tauB, mB, mD, Vcb, ml)
#print "The BR(B) is ",BR


mass = 0.5#1.
U2= 6e-7
a = Integrate3Body(mD, mK, mass, ml, nToys=1000)

#wk = (Vcd**2./Vcs**2.) * (1./mass**2.) * (1./(0.28**2.)) * (1./mDs**3.) * (1./(8.*r.TMath.Pi()**2.)) * a['int']
wk = (1./mass**2.) * (1./fD**2.) * (1./mDs**2.) * (1./mD) * (1./(8.*r.TMath.Pi()**2.)) * a['Int'] / phsp2body(mDs, mass, ml)

print "(D -> K mu N) / (Ds -> mu N) = ", wk

ml = 511.e-6

a1 = Integrate3Body(mD, mK, mass, ml, nToys=1000)

wk = (1./mass**2.) * (1./fD**2.) * (1./mDs**2.) * (1./mD) * (1./(8.*r.TMath.Pi()**2.)) * a1['Int'] / phsp2body(mDs, mass, ml)

print "(D -> K e N) / (Ds -> e N) = ", wk
