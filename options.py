from __future__ import division
import sys
import math
import ROOT as r
import numpy as np
from array import array
from scipy import interpolate
import tauToN


def roundToN(x, n=2):
    if x:
        result = round(x, -int(math.floor(math.log10(math.fabs(x)))) + (n - 1))
    else:
        result = 0.
    return result


class physicsParameters():
    """ Set the physics parameters here """
    def __init__(self):
        self.charmSourceFile = 'CharmFixTarget.root'
        self.beautySourceFile = 'BeautyFixTarget.root'
        self.charmTreeName = 'newTree'
        # Contents of the source file
        self.nD = 5635
        self.nD0 = 11465
        self.nDs = 1553
        self.nB0 = 1027
        self.nB = 1088
        self.nBs = 173
        self.nTotCharm = self.nD + self.nD0 + self.nDs
        #############################
        self.MeV = 1. #everything in GeV for now
        self.modelsSq = [(52.,1.,1.), (1.,16.,3.8), (0.061,1.,4.3)]
        models = self.modelsSq#[tuple([math.sqrt(i) for i in model]) for model in self.modelsSq]
        self.models = [tuple([c/models[l][l] for c in models[l]]) for l in [0,1,2]]
        self.factors = [1. + 1./26.,
                        1. + 1./16. + 3.8/16.,
                        1. + 0.061/4.3 + 1./4.3]
        self.processes = ['Ds -> mu N,N -> pi mu', 'Ds -> mu N,N -> nu nu nu']
        self.decays = ['N -> nu nu nu',
                        'N -> e e nu',
                        'N -> e mu nu',
                        'N -> mu mu nu',
                        'N -> tau tau nu',
                        'N -> e tau nu',
                        'N -> mu tau nu',
                        'N -> pi0 nu',
                        'N -> pi e',
                        'N -> pi mu',
                        'N -> rho nu',
                        'N -> rho e',
                        'N -> rho mu',
                        'rho -> pi pi',
                        'rho -> pi pi0',
                        'pi0 -> gamma gamma',
                        'Ds -> mu N',
                        'Ds -> nu tau',
                        'Ds -> e N',
                        'B -> e N', 'B -> mu N',
                        'tau -> mu N nu', 'tau -> e N nu',
                        'D -> K mu N', 'D -> K e N',
                        'B -> D0 e N', 'B -> D0 mu N',
                        'B0 -> D e N', 'B0 -> D mu N',
                        'Bs -> Ds e N', 'Bs -> Ds mu N']
        self.masses = {'mu':0.10565*self.MeV,
                    'e':0.000511*self.MeV,
                    'tau':1.7768*self.MeV,
                    'pi':0.13957*self.MeV,
                    'pi0':0.13498*self.MeV,
                    'eta':0.54785*self.MeV,
                    'omega':0.78265*self.MeV,
                    'eta1':0.95778*self.MeV,
                    'K':0.493677*self.MeV,
                    'Ds':1.9685*self.MeV,
                    'D':1.86962*self.MeV,
                    'D0':1.865,
                    'Dst':2.007,
                    'p':0.938*self.MeV,
                    'rho':0.775*self.MeV,
                    'nu':0.*self.MeV,
                    'gamma':0.*self.MeV,
                    'up':2.3e-3,
                    'down':4.8e-3,
                    'strange':95.e-3,
                    'charm':1.28,
                    'bottom':4.18,
                    'top':173.,
                    'W':80.39,
                    'Z':91.19,
                    'B':5.279,
                    'Bs':5.367,
                    'B0':5.280,
                    'Bst':5.325,
                    'Bsst':5.415}
        self.name2particle = {'pi':'pi','pi1':'pi', 'pi2':'pi','piTag':'pi',
                            'pi0':'pi0',
                            'mu':'mu','muTag':'mu',
                            'e':'e','eTag':'e',
                            'tau':'tau','tauTag':'tau',
                            'K':'K',
                            'eta':'eta',
                            'omega':'omega',
                            'eta1':'eta1',
                            'X':'X',
                            'N':'N',
                            'D':'D',
                            'Ds':'Ds',
                            'Dst':'Dst',
                            'p':'p',
                            'rho':'rho',
                            'nu':'nu',
                            'gamma':'gamma',
                            'up':'up',
                            'down':'down',
                            'strange':'strange',
                            'charm':'charm',
                            'top':'top',
                            'bottom':'bottom',
                            'W':'W',
                            'Z':'Z', 'Z0':'Z',
                            'B':'B',
                            'Bs':'Bs',
                            'B0':'B0'}
        self.particle2id = {'D':411,
                            'D0':421, # NOTA: ATTENZIONE QUESTO SAREBBE D0
                            'Ds':431}
        self.alphaQED = 1./137.
        self.heV = 6.58211928*pow(10.,-16)
        self.hGeV = self.heV * pow(10.,-9)
        self.c = 3. * pow(10.,8)
        self.GF = 1.166379e-05 # Fermi's constant (GeV^-2)
        self.s2thetaw = 0.23126 # square sine of the Weinberg angle
        self.fpi = 0.1307 # from http://pdg.lbl.gov/2006/reviews/decaycons_s808.pdf
        self.fD = 0.2226
        self.fDs = 0.2801
        self.fB = 0.190 # from gorbunov, shaposhnikov
        self.fBs = 0.230 # from gorbunov, shaposhnikov
        self.tauB = 1.64e-12
        self.tauB0 = 1.52e-12
        self.tauBs = 1.52e-12
        self.tauD = 1.e-12 #sec
        self.tauDs = 0.5e-12
        self.tauTau = 2.91e-13 #sec
        self.decayConstant = {'pi':0.1307, #GeV
                            'pi0':0.130, #GeV
                            'rho':0.102, #GeV
                            'eta':1.2*0.130, #GeV
                            'eta1':-0.45*0.130} #GeV^2
        self.CKM = CKMmatrix()
        self.CKMelemSq = {'pi':self.CKM.Vud**2.,
                        'rho':self.CKM.Vud**2.,
                        'K':self.CKM.Vus**2.,
                        # Quarks:
                        # 0=u, 1=d, 2=s, 3=c, 4=b, 5=t
                        (0,0):1., (1,1):1., (2,2):1., (3,3):1., (4,4):1., (5,5):1., #flavour with same flavour
                        (0,3):0., (3,0):0., (0,5):0., (5,0):0., (3,5):0., (5,3):0., #up-type with up-type
                        (1,2):0., (2,1):0., (1,4):0., (4,1):0., (2,4):0., (4,2):0., #down-type with down-type
                        (0,1):self.CKM.Vud**2., (1,0):self.CKM.Vud**2.,
                        (0,2):self.CKM.Vus**2., (2,0):self.CKM.Vus**2.,
                        (0,4):self.CKM.Vub**2., (4,0):self.CKM.Vub**2.,
                        (3,1):self.CKM.Vcd**2., (1,3):self.CKM.Vcd**2.,
                        (3,2):self.CKM.Vcs**2., (2,3):self.CKM.Vcs**2.,
                        (3,4):self.CKM.Vcb**2., (4,3):self.CKM.Vcb**2.,
                        (5,1):self.CKM.Vtd**2., (1,5):self.CKM.Vtd**2.,
                        (5,2):self.CKM.Vts**2., (2,5):self.CKM.Vts**2.,
                        (5,4):self.CKM.Vtb**2., (4,5):self.CKM.Vtb**2.}
        self.Xcc = 0.45e-03
        self.Xbb = 4.52e-8
        self.BRDsToTau = 0.0543
        #self.w2body = {}
        self.w3body = {}
        #self.lifetimeFun = interpNLifetime('NLifetime.dat')

    def fplus(self, family, q2):
        if family == 'D':
            Mv = self.masses['Dst']
            f0 = 0.727
        elif family == 'B':
            Mv = self.masses['Bst']
            f0 = 0.4
        elif family == 'Bs':
            Mv = self.masses['Bsst']
            f0 = 0.4
        else:
            print 'fplus: unknown meson family!'
            sys.exit(-1)
        val = (Mv**2/(Mv**2-q2))*f0
        return val
    def fzero(self, family, q2):
        if family == 'D':
            Ms = self.masses['D']
            f0 = 0.727
        elif family == 'B':
            Ms = self.masses['B']
            f0 = 0.4
        elif family == 'Bs':
            Ms = self.masses['Bs']
            f0 = 0.4
        else:
            print 'fzero: unknown meson family!'
            sys.exit(-1)
        val = (Ms**2/(Ms**2-q2))*f0
        return val
    def fminus(self, family, q2, mH, mh):
        val = (self.fzero(family, q2)-self.fplus(family, q2))*(mH**2-mh**2)/q2
        return val
    def PhaseSpace3Body(self, family, q2, EN, mH, mh, mN, ml):
        fp = self.fplus(family, q2)
        fm = self.fminus(family, q2, mH, mh)
        PhaseSpace = (fm**2)*(q2*(mN**2+ml**2)-(mN**2-ml**2)**2)+2.*fp*fm*((mN**2)*(2*mH**2-2*mh**2-4*EN*mH-ml**2+mN**2+q2)+(ml**2)*(4.*EN*mH+ml**2-mN**2-q2))
        PhaseSpace += (fp**2)*((4.*EN*mH+ml**2-mN**2-q2)*(2.*mH**2-2.*mh**2-4.*EN*mH-ml**2+mN**2+q2)-(2.*mH**2+2.*mh**2-q2)*(q2-mN**2-ml**2))
        return PhaseSpace
    def phsp2body(self, mH, mN, ml):
        if mN >= (mH-ml):
            return 0.
        p1 = 1. - mN**2./mH**2. + 2.*(ml**2./mH**2.) + (ml**2./mN**2.)*(1. - (ml**2./mH**2.))
        p2 = ( 1. + (mN**2./mH**2.) - (ml**2./mH**2.) )**2. - 4.*(mN**2./mH**2.)
        p3 = math.sqrt( p2 )
        return p1*p3
    def Integrate3Body(self, family, mH, mh, mN, ml, nToys=1000):
        #### Setting up the parameters for generation
        W = r.TLorentzVector(0., 0., 0., mH)
        masses = array('d', [mh, mN, ml])
        event = r.TGenPhaseSpace()
        event.SetDecay(W, 3, masses)
        Nq2 = 20 
        NEN = 20
        ENMax = (mH-mh)
        ENMax = r.TMath.Sqrt(((ENMax**2-ml**2-mN**2)/2.)**2+mN**2) # converte massa max ad energia max
        hist = r.TH2F("hist", "", Nq2, (ml+mN)**2, (mH-mh)**2, NEN, mN, ENMax)
        ###### Integral
        Integral = 0.
        #### For loop in order to integrate
        for i in xrange(nToys):
            event.Generate()
            #### Getting momentum of the daughters
            ph = event.GetDecay(0)
            pN = event.GetDecay(1)
            pl = event.GetDecay(2)
            #### Computing the parameters to compute the phase space
            q = pl+pN
            q2 = q.M2()
            EN = pN.E()
            iBin = hist.Fill(q2, EN)
            if hist.GetBinContent(iBin)==1.:
                val = self.PhaseSpace3Body(family, q2, EN, mH, mh, mN, ml)
                Integral+=val*hist.GetXaxis().GetBinWidth(1)*hist.GetYaxis().GetBinWidth(1)
        hist.Delete()    
        return Integral
    def computeProductionWeights(self, lepton):
        #if self.MN > self.masses[self.name2particle['Ds']] - self.masses[self.name2particle[lepton]]:
        #    self.w2body[lepton] = 0.
        #else:
        if self.MN > self.masses[self.name2particle['D']] - self.masses[self.name2particle[lepton]] - self.masses[self.name2particle['K']]:
            self.w3body[lepton] = 0.
            #self.w2body[lepton] = 1.
        else:
            wk = (1./self.MN**2.) * (1./self.fDs**2.) * (1./self.masses[self.name2particle['D']]**2.) * (1./self.masses[self.name2particle['Ds']]) * (1./(8.*math.pi**2.))
            wk *= (self.tauD/self.tauDs) #* (self.CKM.Vcd**2./self.CKM.Vcs**2.) #no, e' lo stesso Vcs!!!
            wk *= self.Integrate3Body('D',self.masses[self.name2particle['D']], self.masses[self.name2particle['K']], self.MN, self.masses[self.name2particle[lepton]], nToys=300)
            wk = wk / self.phsp2body(self.masses[self.name2particle['Ds']], self.MN, self.masses[self.name2particle[lepton]])
            if wk > 1.:
                wk = 1.
            self.w3body[lepton] = wk
            #self.w2body[lepton] = 1. - wk


    def setNCoupling(self, couplings):
        self.U2 = couplings
        self.U = [math.sqrt(ui) for ui in self.U2]
    def computeU2tot(self):
        self.U2tot = sum(self.U2)
        return self.U2tot
    def setNMass(self, mass):
        self.MN = mass
        self.NM = mass
        if 'N' in self.masses:
            self.masses['N'] = self.MN
        else:
            self.masses.update({'N':self.MN})
    def computeNProdBR(self, alpha): # only for 2-body production!!
        #alpha_BR = 3.6*1e-8/(6e-7)
        #self.BR = alpha_BR*self.U2[alpha]
        if alpha<2: #prod from Ue, Umu
            l = ['e','mu','tau']
            #alpha_BR = (self.tauD / self.hGeV) * (self.GF**2. * self.fD**2. * self.masses['D']**2. * self.CKM.Vcd**2.) / (8. * math.pi)
            alpha_BR = (self.tauDs / self.hGeV) * (self.GF**2. * self.fDs**2. * self.masses['Ds'] * self.CKM.Vcs**2.) / (8. * math.pi)
            ps = self.phsp2body(self.masses[self.name2particle['Ds']], self.MN, self.masses[self.name2particle[l[alpha]]])
            self.BR = ps * alpha_BR * self.U2[alpha] * self.MN**2. * 2. #Majorana neutrinos --> x2!
        elif alpha==2: #prod from Utau
            self.BR = tauToN.brTauToNuEllN(self,'e') + tauToN.brTauToNuEllN(self,'mu') + tauToN.brTauToPiN(self)
        return self.BR
    #def computeNLifetime(self):
    #    alpha_LT = 3*1e-6*6e-7
    #    self.LT = alpha_LT/self.U2
    #    return self.LT
    def computeNLifetime(self):
        #self.NLifetime = self.lifetimeFun(self.NM, self.U)
        # Da slides Ruchaivski workshop:
        #self.NLifetime = 0.3 * pow(self.MN, -4.) # seconds 
        #self.NLifetime = self.hGeV*86.*(math.pi**3.) / ( (self.U2)*(self.GF**2.)*(self.MN**5.) ) # seconds 
        self.NLifetime = self.hGeV / self.NDecayWidth()
        return self.NLifetime
    def drawNDecayVtx(self, momentum, originVtx):
        speed = momentum.Beta()*self.c
        DeltaT = momentum.Gamma()*r.gRandom.Exp(self.NLifetime)
        DeltaL = speed*DeltaT
        direction = momentum.Vect().Unit()
        endVtx = r.TVector3(direction.X()*DeltaL, direction.Y()*DeltaL, direction.Z()*DeltaL)
        decayVtx = originVtx + endVtx
        return decayVtx
    def energy(self, p, m):
        return math.sqrt(p*p + m*m)
    def momentum(self, E, m):
        return math.sqrt(E*E - m*m)
    def Width_H0_nu(self, H, alpha):
        if self.MN < (self.masses[self.name2particle[H]]):
            return 0.
        width = (math.fabs(self.U2[alpha-1])/(32.*math.pi))*(self.GF**2.)*(self.decayConstant[self.name2particle[H]]**2.)
        par = (self.MN**3.)*((1 - ((self.masses[self.name2particle[H]]**2.)/(self.MN**2.)))**2.)
        if self.name2particle[H] == 'rho':
            par = par*2./(self.masses[self.name2particle[H]]**2.)
            par = par*(1 + 2.*((self.masses[self.name2particle[H]]**2.)/(self.MN**2.)))
        width = width*par
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    def Width_H_l(self, H, alpha):
        l = [None,'e','mu','tau']
        if self.MN < (self.masses[self.name2particle[H]] + self.masses[self.name2particle[l[alpha]]]):
            return 0.
        width = (math.fabs(self.U2[alpha-1])/(16.*math.pi))*(self.GF**2.)*(self.decayConstant[self.name2particle[H]]**2.)
        width = width*(self.MN**3.)*self.CKMelemSq[self.name2particle[H]]
        par = ( ((1 - ((self.masses[self.name2particle[l[alpha]]]**2.)/(self.MN**2.)))**2.) 
                - ( (self.masses[self.name2particle[H]]**2.)/(self.MN**2.)
                * (1 + ((self.masses[self.name2particle[l[alpha]]]**2.)/(self.MN**2.))) ) )
        if self.name2particle[H] == 'rho':
            par = ( ((1 - ((self.masses[self.name2particle[l[alpha]]]**2.)/(self.MN**2.)))**2.) 
                    + ( (self.masses[self.name2particle[H]]**2.)/(self.MN**2.)
                    * (1 + (((self.masses[self.name2particle[l[alpha]]]**2. - 2.*self.masses[self.name2particle[H]]**2.))/(self.MN**2.))) ) )
            par = par*2./(self.masses[self.name2particle[H]]**2.)
            
        width = width*par
        rad = math.sqrt( ( 1-((self.masses[self.name2particle[H]]-self.masses[self.name2particle[l[alpha]]])**2.)/(self.MN**2.) )
                * ( ( 1-((self.masses[self.name2particle[H]]+self.masses[self.name2particle[l[alpha]]])**2.)/(self.MN**2.) ) ) )
        width = width*rad
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    def Width_3nu(self):
        width = (self.GF**2.)*(self.MN**5.)*sum(self.U2)/(192.*(math.pi**3.))
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    def Width_l1_l2_nu(self, alpha, beta, gamma): # alpha, beta for the two leptons, gamma for the neutrino
        l = [None,'e','mu','tau','up','down','strange','charm','bottom','top']
        if self.MN < (self.masses[self.name2particle[l[alpha]]] + self.masses[self.name2particle[l[beta]]]):
            return 0.

        if (alpha > 3) and (beta == alpha): # N -> nu qq or N -> nu ll (mainly Z0)
            width = (self.GF**2.)*(self.MN**5.)*self.U2[gamma-1]/(192.*(math.pi**3.))

        elif ((alpha in [4, 7, 9]) and (beta in [5, 6, 8])) or ((alpha in [5, 6, 8]) and (beta in [4, 7, 9])): # N -> l W, W -> u dbar, index gamma stands for massive lepton
            if gamma == 3: # tau too massive for this parametrisation
                return 0.
            if self.MN < (self.masses[self.name2particle[l[alpha]]] + self.masses[self.name2particle[l[beta]]] + self.masses[self.name2particle[l[gamma]]]):
                return 0.
            width = (self.GF**2.)*(self.MN**5.)*self.U2[gamma-1]/(192.*(math.pi**3.))*self.CKMelemSq[(alpha-4, beta-4)]

        elif (alpha <= 3) and (beta <= 3):
            width = (self.GF**2.)*(self.MN**5.)*self.U2[alpha-1]/(192.*(math.pi**3.))

        else:
            return 0.

        if alpha != beta:
            xl = max([self.masses[self.name2particle[l[alpha]]] , self.masses[self.name2particle[l[beta]]]])/self.MN
            width = width*(1. - 8.*xl**2. + 8.*xl**6. - (12.*xl**4.)*math.log(xl**2.))
        else:
            xl = self.masses[self.name2particle[l[alpha]]] / self.MN
            lo = 0.
            logContent = (1. - 3.*xl**2. - (1.-xl**2.)*math.sqrt(1. - 4.*xl**2.) ) / ( (xl**2.)*(1 + math.sqrt(1. - 4.*xl**2.)) )
            if logContent > 0:
                lo = math.log( logContent )
            c1 = 0.25*(1. - 4.*self.s2thetaw + 8.*self.s2thetaw**2.)
            c2 = 0.5*self.s2thetaw*(2.*self.s2thetaw -1.)
            c3 = 0.25*(1. + 4.*self.s2thetaw + 8.*self.s2thetaw**2.)
            c4 = 0.5*self.s2thetaw*(2.*self.s2thetaw +1.)
            d = (alpha == gamma)
            width = width*( (c1*(1-d)+c3*d)*( (1.-14.*xl**2. -2.*xl**4. -12.*xl**6.)*math.sqrt(1-4.*xl**2) +12.*xl**4. *(-1.+xl**4.)*lo )
                            + 4.*(c2*(1-d)+c4*d)*( xl**2. *(2.+10.*xl**2. -12.*xl**4.) * math.sqrt(1.-4.*xl**2) + 6.*xl**4. *(1.-2.*xl**2+2.*xl**4)*lo ) )
        #if alpha>3 and beta>3: # N -> q q nu
        #    width = width*self.CKMelemSq[(alpha-4, beta-4)]
        width = 2.*width # Majorana case (charge conjugate channels)
        return width

    def NDecayWidth(self):
        if self.MN < 1.:
            totalWidth = ( self.Width_3nu()
                    + sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])
                    + sum([self.Width_l1_l2_nu(a,b,g) for a in [1,2,3] for b in [1,2,3] for g in [1,2,3]]) )
        elif self.MN > 2.:
            totalWidth = ( self.Width_3nu() + sum([self.Width_l1_l2_nu(a,b,g) for a in range(1,10) for b in range(1,10) for g in [1,2,3]]) )
        else:
            m1, m2 = 1., 2.
            tempMass = self.MN
            self.MN = m1
            w1 = ( self.Width_3nu()
                    + sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])
                    + sum([self.Width_l1_l2_nu(a,b,g) for a in [1,2,3] for b in [1,2,3] for g in [1,2,3]]) )
            self.MN = m2
            w2 = ( self.Width_3nu() + sum([self.Width_l1_l2_nu(a,b,g) for a in range(1,10) for b in range(1,10) for g in [1,2,3]]) )
            self.MN = tempMass
            totalWidth = w1 + (self.MN - m1)*(w2 - w1)/(m2 - m1)
        #if self.MN < 10.:
        #    totalWidth = ( self.Width_3nu()
        #            + sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])
        #            + sum([self.Width_l1_l2_nu(a,b,g) for a in [1,2,3] for b in [1,2,3] for g in [1,2,3]]) )
        #else:
        #    totalWidth = ( self.Width_3nu() + sum([self.Width_l1_l2_nu(a,b,g) for a in range(1,10) for b in range(1,10) for g in [1,2,3]]) )
        return totalWidth
    def findBranchingRatio(self, decayString):
        br = 0.
        totalWidth = self.NDecayWidth()
        if decayString == 'N -> pi e':
            br = self.Width_H_l('pi',1) / totalWidth
        elif decayString == 'N -> pi0 nu' or decayString == 'N -> pi nu':
            br = sum([self.Width_H0_nu('pi0',l) for l in [1,2,3]]) / totalWidth
        elif decayString == 'N -> pi mu':
            br = self.Width_H_l('pi',2) / totalWidth
        elif decayString == 'N -> rho nu' or decayString == 'N -> rho0 nu':
            br = sum([self.Width_H0_nu('rho',l) for l in [1,2,3]]) / totalWidth
        elif decayString == 'N -> rho e':
            br = self.Width_H_l('rho',1) / totalWidth
        elif decayString == 'N -> rho mu':
            br = self.Width_H_l('rho',2) / totalWidth
        elif decayString == 'N -> e e nu':
            br = sum([self.Width_l1_l2_nu(1,1,l) for l in [1,2,3]]) / totalWidth
        elif decayString == 'N -> mu mu nu':
            br = sum([self.Width_l1_l2_nu(2,2,l) for l in [1,2,3]]) / totalWidth
        elif decayString == 'N -> tau tau nu':
            br = sum([self.Width_l1_l2_nu(3,3,l) for l in [1,2,3]]) / totalWidth
        elif decayString == 'N -> e mu nu':
            br = (sum([self.Width_l1_l2_nu(1,2,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(2,1,l) for l in [1,2,3]])) / totalWidth
        elif decayString == 'N -> e tau nu':
            br = (sum([self.Width_l1_l2_nu(1,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,1,l) for l in [1,2,3]])) / totalWidth
        elif decayString == 'N -> mu tau nu':
            br = (sum([self.Width_l1_l2_nu(2,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,2,l) for l in [1,2,3]])) / totalWidth
        elif decayString == 'N -> nu nu nu' or decayString == 'N -> 3nu':
            br = self.Width_3nu() / totalWidth
        elif decayString == 'N -> hadrons':
            if self.MN < 1.:
                br = sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])/totalWidth
            elif self.MN > 2.:
                br = sum([self.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/totalWidth
            else:
                m1, m2 = 1., 2.
                tempMass = self.MN
                self.MN = m1
                br1 = sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])/totalWidth
                self.MN = m2
                br2 = sum([self.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/totalWidth
                self.MN = tempMass
                br = br1 + (self.MN - m1)*(br2 - br1)/(m2 - m1)
        elif decayString == 'N -> charged hadrons':
            if self.MN < 1.:
                br = sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) for l in [1,2,3]])/totalWidth
            elif self.MN > 2.:
                br = sum([self.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/totalWidth
            else:
                m1, m2 = 1., 2.
                tempMass = self.MN
                self.MN = m1
                br1 = sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) for l in [1,2,3]])/totalWidth
                self.MN = m2
                br2 = sum([self.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/totalWidth
                self.MN = tempMass
                br = br1 + (self.MN - m1)*(br2 - br1)/(m2 - m1)
            #if self.MN < 10.:
            #    br = sum([self.Width_H_l('pi',l) + self.Width_H_l('rho',l) + self.Width_H0_nu('pi0',l) + self.Width_H0_nu('rho',l) + self.Width_H0_nu('eta',l) + self.Width_H0_nu('eta1',l) for l in [1,2,3]])/totalWidth
            #else:
            #    br = sum([self.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/totalWidth
        elif decayString == 'N -> TLEP visible':
            # eenu, mumunu, emunu
            br0 = ( (sum([self.Width_l1_l2_nu(1,1,l) for l in [1,2,3]]) 
                   + sum([self.Width_l1_l2_nu(2,2,l) for l in [1,2,3]])
                   + (sum([self.Width_l1_l2_nu(1,2,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(2,1,l) for l in [1,2,3]]))
                   )/ totalWidth )
            br = br0
            if self.MN > 2.:
                br1 = ( (sum([self.Width_l1_l2_nu(a,b,g) for a in [4,7,9] for b in [5,6,8] for g in [1,2,3]]) +
                        sum([self.Width_l1_l2_nu(a,b,g) for a in [5,6,8] for b in [4,7,9] for g in [1,2,3]]) +
                        sum([self.Width_l1_l2_nu(3,3,l) for l in [1,2,3]]) + #tau tau nu
                        (sum([self.Width_l1_l2_nu(1,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,1,l) for l in [1,2,3]])) + #e tau nu
                        (sum([self.Width_l1_l2_nu(2,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,2,l) for l in [1,2,3]])) #mu tau nu
                        )/ totalWidth )
                br = br0 + br1
            if 1. <= self.MN <= 2.:
                m1, m2 = 1., 2.
                tempMass = self.MN
                self.MN = m2
                br2 = ( (sum([self.Width_l1_l2_nu(a,b,g) for a in [4,7,9] for b in [5,6,8] for g in [1,2,3]]) +
                        sum([self.Width_l1_l2_nu(a,b,g) for a in [5,6,8] for b in [4,7,9] for g in [1,2,3]]) +
                        sum([self.Width_l1_l2_nu(3,3,l) for l in [1,2,3]]) + #tau tau nu
                        (sum([self.Width_l1_l2_nu(1,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,1,l) for l in [1,2,3]])) + #e tau nu
                        (sum([self.Width_l1_l2_nu(2,3,l) for l in [1,2,3]]) + sum([self.Width_l1_l2_nu(3,2,l) for l in [1,2,3]])) #mu tau nu
                        )/ totalWidth ) + br0
                self.MN = tempMass
                br = br0 + (self.MN - m1)*(br2 - br0)/(m2 - m1)
        else:
            print 'findBranchingRatio ERROR: unknown decay %s'%decayString
            sys.exit(-1)
        return br
    def HNLAllowedDecays(self):
        m = self.MN
        allowedDecays = {'N -> nu nu nu':'yes'}
        if m > 2.*self.masses[self.name2particle['e']]:
            allowedDecays.update({'N -> e e nu':'yes'})
            if m > self.masses[self.name2particle['e']] + self.masses[self.name2particle['mu']]:
                allowedDecays.update({'N -> e mu nu':'yes'})
            if m > self.masses[self.name2particle['pi0']]:
                allowedDecays.update({'N -> pi0 nu':'yes'})
            if m > self.masses[self.name2particle['pi']] + self.masses[self.name2particle['e']]:
                allowedDecays.update({'N -> pi e':'yes'})
                if m > 2.*self.masses[self.name2particle['mu']]:
                    allowedDecays.update({'N -> mu mu nu':'yes'})
                    if m > self.masses[self.name2particle['pi']] + self.masses[self.name2particle['mu']]:
                        allowedDecays.update({'N -> pi mu':'yes'})
                        if m > self.masses[self.name2particle['rho']]:
                            allowedDecays.update({'N -> rho nu':'yes'})
                            if m > self.masses[self.name2particle['rho']] + self.masses[self.name2particle['e']]:
                                allowedDecays.update({'N -> rho e':'yes'})
                                if m > self.masses[self.name2particle['rho']] + self.masses[self.name2particle['mu']]:
                                    allowedDecays.update({'N -> rho mu':'yes'})
                                    if m > self.masses[self.name2particle['e']] + self.masses[self.name2particle['tau']]:
                                        allowedDecays.update({'N -> e tau nu':'yes'})
                                        if m > self.masses[self.name2particle['mu']] + self.masses[self.name2particle['tau']]:
                                            allowedDecays.update({'N -> mu tau nu':'yes'})
                                            if m > 2.*self.masses[self.name2particle['tau']]:
                                                allowedDecays.update({'N -> tau tau nu':'yes'})
        for decay in self.decays:
            if decay not in allowedDecays and decay.startswith('N'):
                allowedDecays.update({decay:'no'})
        return allowedDecays



def expon(x, par):
    return math.exp( (-1.)*float(par[0])*float(x[0]) )


class experimentParams():
    """ The experimental settings (beam, geometry...)
        Requires an instance pp of physicsParameters() """
    def __init__(self, pp, name):
        self.LTfun = r.TF1("lifetime",expon,0., 160.,2)
        if name == "SHiP" or name == "ship" or name == "SHIP":
            self.protonEnergy = 400.*pp.MeV # GeV/c
            self.protonMomentum = pp.momentum(self.protonEnergy, pp.masses['p'])
            self.protonFlux = 2.*pow(10.,20.)
            self.firstVolume = [60., 100., 2.5] # start, end, radius (m)
            self.secondVolume = [110., 150., 2.5] # start, end, radius (m)
            self.v1ThetaMax = self.firstVolume[2]/self.firstVolume[0]
            self.v2ThetaMax = self.secondVolume[2]/self.secondVolume[0]
        elif name == "TLEP" or name == "tlep" or name == "Tlep":
            self.nZ = 1.e12
            self.nW = 1.e11
            self.efficiency = 1.
            self.minSVdistance = 1.e-3 # 1 mm (maybe increase to 5 mm)
            self.Rmin = 1.e-3
            self.Rmax = 1. # 1 m (maybe increase?)
            #self.NGamma = pp.masses[pp.name2particle['Z']]/2.
            self.BRZnunu = 0.08
            self.BRWlnu = 0.108
        else:
            print "experimentParams ERROR: please provide the name of the experiment!"
    def makeVtxInVolume(self, momentum, ct, volume):
        if volume == 1:
            vol = self.firstVolume
        elif volume == 2:
            vol = self.secondVolume
        else:
            print "ERROR: select decay volume 1 or 2"
            return 0
        gamma = momentum.Gamma()
        Direction = momentum.Vect().Unit()
        costheta = math.fabs(momentum.CosTheta())
        dxdz = momentum.Px()/momentum.Pz()
        dydz = momentum.Py()/momentum.Pz()
        Origin = r.TVector3( vol[0]*dxdz, vol[0]*dydz, vol[0] )
        self.LTfun.SetParameter(0, 1./(gamma*ct))
        maxlength = 40.*costheta
        DL = self.LTfun.GetRandom(0., maxlength)
        EndVertex = r.TVector3(Direction[0]*DL, Direction[1]*DL, Direction[2]*DL)
        EndVertex = Origin + EndVertex
        return EndVertex
    def TLEPacceptance(self,pp,boson):
        if (boson is not 'Z') and (boson is not 'W'):
            print 'TLEPacceptance: please define the source of neutrinos (Z or W)!'
            sys.exit(-1)
        lt = pp.computeNLifetime()
        mb = pp.masses[pp.name2particle[boson]]
        gamma = (mb/(2.*pp.MN)) + (pp.MN/(2.*mb)) #this neglects ml in W->lN
        cgammatau = lt*gamma*pp.c
        esp1 = (-1.)*self.Rmin/cgammatau
        esp2 = (-1.)*self.Rmax/cgammatau
        #np.seterr(all='raise')
        try:
            #result = np.nan_to_num(np.fabs( np.exp(esp1) - np.exp(esp2) ))
            result = np.nan_to_num(np.exp(esp1) - np.exp(esp2))
        except (ValueError, FloatingPointError):#, RuntimeWarning):
            result = 0.
        return result
    def inAcceptance(self, vtx, particleList, volume):
        if not particleList:
            return False
        if volume == 1:
            vol = self.firstVolume
            tMax = self.v1ThetaMax
        elif volume == 2:
            vol = self.secondVolume
            tMax = self.v2ThetaMax
        else:
            print "experimentParams.inAcceptance ERROR: please select volume 1 or 2!"
        detectable = []
        if (vtx.Z() > vol[0]) and (vtx.Z() < vol[1]): # longitudinal vtx acc
            if (vtx.X()**2. + vtx.Y()**2.) < vol[2]: # transverse vtx acc
                for particle in particleList:
                    tx = particle.Px() / particle.Pz()
                    ty = particle.Py() / particle.Pz()
                    endPos1 = r.TVector3()
                    endPos1.SetZ(vol[1])
                    endPos1.SetX( vtx.X() + tx*(endPos1.Z() - vtx.Z()) )
                    endPos1.SetY( vtx.Y() + ty*(endPos1.Z() - vtx.Z()) )
                    #print "mass: ", math.sqrt(math.fabs(particle.E()**2. - particle.P()**2.)), " endPos1 ", endPos1.X(), endPos1.Y(), endPos1.Z(), ((endPos1.X()**2. + endPos1.Y()**2.) < vol[2]**2.)
                    if (endPos1.X()**2. + endPos1.Y()**2.) < vol[2]**2.:
                        detectable.append(True)
                    else:
                        detectable.append(False)
                return bool(np.product(detectable))
                    ## Check if child 2 is detectable:
                    #tx2 = mom2.Px() / mom2.Pz()
                    #ty2 = mom2.Py() / mom2.Pz()
                    #endPos2 = r.TVector3()
                    #endPos2.SetZ(vol[1])
                    #endPos2.SetX( vtx.X() + tx2*(endPos2.Z() - vtx.Z()) )
                    #endPos2.SetY( vtx.Y() + ty2*(endPos2.Z() - vtx.Z()) )
                    #if (endPos2.X()**2. + endPos2.Y()**2.) < vol[2]:
                    #    return True
        # Otherwise
        return False
    def probVtxInVolume(self, momentum, ct, volume):
        if volume == 1:
            vol = self.firstVolume
        elif volume == 2:
            vol = self.secondVolume
        else:
            print "ERROR: select decay volume 1 or 2"
            return 0
        gamma = momentum.Gamma()
        costheta = np.fabs(momentum.CosTheta())
        tantheta = np.tan(momentum.Theta())
        start = vol[0]
        end = vol[1]
        rad = vol[2]
        stop = min(rad/tantheta, end)
        if stop < start:
            return 0.
        esp1 = (-1.) * (start/costheta) / (gamma*ct)
        esp2 = (-1.) * (stop/costheta) / (gamma*ct)
        #np.seterr(all='raise')
        try:
            result = np.nan_to_num(np.fabs( np.exp(esp1) - np.exp(esp2) ))
        except (ValueError, FloatingPointError):#, RuntimeWarning):
            result = 0.
        return result
    def GeometricAcceptance(self, momentum, volume):
        if volume == 1:
            vol = self.firstVolume
        elif volume == 2:
            vol = self.secondVolume
        else:
            print "ERROR: select decay volume 1 or 2"
            return 0
        #if (math.fabs((px/pz)*vol[1])<vol[2]) and (pz>0):
        #    return True
        if momentum.Theta() < vol[2]/vol[0] and momentum.Pz() > 0.:
            return True
        return False

        

class decayNBody():
    """ General 2-bodies decay representation.
        Requires an instance pp of physicsParameters(). """
    def __init__(self, pp):
        self.pp = pp
        self.mother = None
        self.kid1 = None
        self.kid2 = None
        self.particles = [None]
        self.lifetime = None
        self.gen = r.TGenPhaseSpace()
    def readString(self, decayName):
        if decayName not in self.pp.decays:# + [process.split(',')[0] for process in self.pp.processes]:
            print 'decayNBody::readString ERROR: decay %s not found in database!'%decayName
            sys.exit(-1)
        self.name = decayName
        self.particles = [self.pp.name2particle[p] for p in self.name.replace('->',' ').split()]
        self.mother = self.particles[0]
        self.motherMass = self.pp.masses[self.mother]
        self.pMother = r.TLorentzVector(0., 0., 0., self.motherMass)
        self.kid1 = self.particles[1]
        self.kid2 = self.particles[2]
        self.childrenMasses = array('d', [self.pp.masses[self.kid1], self.pp.masses[self.kid2]])
        self.nbody = 2
        if self.name in (self.pp.decays[:7]+self.pp.decays[21:]):
            self.nbody = 3
            self.kid3 = self.particles[3]
            self.childrenMasses = array('d', [self.pp.masses[self.kid1], self.pp.masses[self.kid2], self.pp.masses[self.kid3]])
        if sum(self.childrenMasses) > self.motherMass:
            print 'decayNBody::readString ERROR: children too heavy! %s'%decayName
            sys.exit(-1)
    def setPMother(self, pMother):
        self.pMother = pMother
    def setLifetime(self, lifetime):
        self.lifetime = lifetime
    def makeDecay(self):
        """ if no setPMother method is called, this method makes a RF decay """
        if not self.childrenMasses:
            print "decayNBody ERROR: I have no particles!"
            sys.exit(-1)
        self.gen.SetDecay(self.pMother, self.nbody, self.childrenMasses)
        self.gen.Generate()
        self.pKid1, self.pKid2 = self.gen.GetDecay(0), self.gen.GetDecay(1)
        #if self.name in self.pp.decays[:4]:
        if self.nbody == 3:
            self.pKid3 = self.gen.GetDecay(2)
            return self.pKid1, self.pKid2, self.pKid3
        return self.pKid1, self.pKid2
    def boostChildren(self, pBoost):
        """ requires a TLorentzVector::BoostVector pBoost """
        self.pKid1.Boost(pBoost)
        self.pKid2.Boost(pBoost)

class myNEvent():
    """ Sterile neutrino production and decay.
        Requires an instance pp of physicsParameters(). """
    def __init__(self, pp):#, evName):
        if not 'N' in pp.masses:
            print "physicsParameters ERROR: sterile neutrino mass undefined!"
            sys.exit(-1)
        if not pp.U:
            print "physicsParameters ERROR: sterile neutrino coupling undefined!"
            sys.exit(-1)
        self.pp = pp #physics
        #self.evName = evName
        #if not evName in self.pp.processes:
        #    print "myNEvent ERROR: decay chain not found in db!"
        #    sys.exit(-1)
        #self.prodEvtName, self.NdecayName = self.evName.split(',')
        self.production = decayNBody(self.pp)
        #self.production.readString(self.prodEvtName)
        self.decay = decayNBody(self.pp)
        #self.decay.readString(self.NdecayName)
        #self.rootFileName = "out/NTuples/%s_U%s_m%s.root"%(self.evName,'_'.join([str(ui) for ui in self.pp.U2]),self.pp.MN)






class CKMmatrix():
    # From http://pdg.lbl.gov/2013/reviews/rpp2013-rev-ckm-matrix.pdf
    def __init__(self):
        self.Vud = 0.9742
        self.Vus = 0.2252
        self.Vub = 4.2e-03
        self.Vcd = 0.23
        self.Vcs = 1.
        self.Vcb = 4.1e-02
        self.Vtd = 8.4e-03
        self.Vts = 4.3e-02
        self.Vtb = 0.89
















def interpNLifetime(infile):
    """ DA SISTEMARE """
    NM, Umu, points, lifetime = [], [], [], []
    with open(infile,'r') as f:
        for line in f:
            line = line.split()
            try:
                val_NM = float(line[0])
                val_Umu = float(line[1])
                val_lifetime = float(line[2])
                NM.append(roundToN(np.log10(val_NM),4))
                Umu.append(roundToN(np.log10(val_Umu),4))
                points.append((val_NM, val_Umu))
                lifetime.append(np.log10(val_lifetime))
            except:
                continue
    print len(NM), len(Umu), len(lifetime)
    #NM, Umu, lifetime = zip(* sorted( zip(NM, Umu, lifetime), key=lambda x:(x[0], -x[1]) ) )
    NM, Umu, lifetime = zip(* sorted( zip(NM, Umu, lifetime), key=lambda x:x[1] ) )
    NM, Umu, lifetime = list(NM), list(Umu), list(lifetime)
    NM, Umu, lifetime = zip(* sorted( zip(NM, Umu, lifetime), key=lambda x:x[0] ) )
    NM, Umu, lifetime = list(NM), list(Umu), list(lifetime)
    print NM[0], NM[-1]
    print Umu[0], Umu[-1]
    print lifetime[0], lifetime[-1]
    #gridx, gridy = np.mgrid[NM[0]:NM[-1]:(np.fabs(NM[-1]-NM[0])/len(NM)), Umu[0]:Umu[-1]:(np.fabs(Umu[-1]-Umu[0])/len(Umu))]
    #print gridx, gridy
    #fun = interpolate.interp2d(NM, Umu, lifetime, kind='cubic')
    #fun = interpolate.SmoothBivariateSpline(NM, Umu, lifetime)
    #print Umu
    griddedM = list(set(NM[:]))
    print griddedM#, len(griddedM)
    griddedUmu = list(set(Umu[:]))
    print griddedUmu#, len(griddedUmu)
    griddedLt = np.array(lifetime)
    griddedLt = griddedLt.reshape(len(griddedM), len(griddedUmu))
    print griddedLt.shape
#    print len(set(NM)), len(set(Umu))
    #grid = np.zeros(len(set(NM)), len(set(Umu)), 'f')
#    grid[griddedM, griddedUmu] = griddedLt
    fun = interpolate.RectBivariateSpline(griddedM, griddedUmu, griddedLt)
    #griddedLF = interpolate.griddata((NM, Umu), lifetime, (gridx, gridy), method='cubic')
    #print type(griddedLF), len(griddedLF), len(gridx), len(gridy)
    ##fun = interpolate.interp2d(gridx, gridy, griddedLF)
    #fun = interpolate.RectBivariateSpline(gridx, gridy, griddedLF)
    return fun
