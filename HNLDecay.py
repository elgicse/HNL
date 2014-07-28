from array import array
from options import *

"""Show branching fractions"""
#usage: makeBRPlot([(0.25*0.25)*1.e-8,1.e-8,0.25e-8])
cSaver = []
#r.gROOT.ProcessLine(".x lhcbstyle.C")

def BAU(m,pp):
    return 2.e-6*(1/m**2.)*((1.-(m**2./pp.masses[pp.name2particle['W']]**2.))**2.)

#fBau = r.TF1("Bau", "3.*(1e-8)*(1/x)", 1., 89.)

def makeTLEPacceptance(model, mass, logU2min, logU2max):
    pp = physicsParameters()
    pp.setNMass(mass)
    ep = experimentParams(pp, 'TLEP')
    if pp.MN > pp.masses[pp.name2particle['Z']]:
        print 'Too heavy for a Z factory!'
        return 0.
    model = model - 1
    if model not in xrange(3):
        print 'Please select model 1, 2 or 3'
        return 0.
    mz = pp.masses[pp.name2particle['Z']]
    gamma = (mz/(2.*pp.MN)) + (pp.MN/(2.*mz))
    coups = np.logspace(logU2min,logU2max,300).tolist()
    accsV1, accsV2, lt, esp1, esp2 = [], [], [], [], []
    for coupling in coups:
        if model == 0:
            couplings = [coupling, pp.models[model][1]*coupling, pp.models[model][2]*coupling]
        elif model == 1:
            couplings = [pp.models[model][0]*coupling, coupling, pp.models[model][2]*coupling]
        elif model == 2:
            couplings = [pp.models[model][0]*coupling, pp.models[model][1]*coupling, coupling]
        pp.setNCoupling(couplings)
        lt.append(pp.computeNLifetime())
        accsV1.append(ep.TLEPacceptance(pp))
        ep.Rmin = 1.e-4
        ep.Rmax = 5.
        accsV2.append(ep.TLEPacceptance(pp))
        ep.Rmin = 1.e-3
        ep.Rmax = 1.
        try:
            esp1.append(np.exp(-1.*1.e-3/(gamma*pp.c*pp.computeNLifetime())) - np.exp(-1./(gamma*pp.c*pp.computeNLifetime())))
        except FloatingPointError:
            esp1.append(0.0)
        try:
            esp2.append(np.exp(-1.*1.e-4/(gamma*pp.c*pp.computeNLifetime())) - np.exp(-5./(gamma*pp.c*pp.computeNLifetime())))
        except FloatingPointError:
            esp2.append(0.0)
    accgrV1 = r.TGraph(len(coups),array('f',coups),array('f',accsV1))
    accgrV2 = r.TGraph(len(coups),array('f',coups),array('f',accsV2))
    c1 = r.TCanvas()
    cSaver.append(c1)
    accgrV1.SetLineWidth(3)
    accgrV1.SetLineColor(r.kBlue)
    accgrV1.SetMarkerColor(r.kBlue)
    accgrV2.SetLineWidth(3)
    accgrV2.SetLineColor(r.kRed)
    accgrV2.SetMarkerColor(r.kRed)
    accgr = r.TMultiGraph()
    accgr.Add(accgrV1)
    accgr.Add(accgrV2)
    c1.SetLogy()
    accgr.Draw('alp')
    accgr.SetTitle('TLEP acceptance for HNL of mass M_{N}=%s GeV'%mass)
    accgr.GetXaxis().SetTitle('U_{%s}^{2}'%(model+1))
    accgr.GetXaxis().SetLimits(10.**logU2min, 10.**logU2max)
    accgr.GetYaxis().SetTitle('P(vtx in detector)')
    c1.SetLogx()
    c1.Update()

    mgesp = r.TMultiGraph()
    gresp1 = r.TGraph(len(coups),array('f',coups),array('f',esp1))
    gresp1.SetLineColor(r.kBlue)
    gresp1.SetMarkerColor(r.kBlue)
    gresp2 = r.TGraph(len(coups),array('f',coups),array('f',esp2))
    gresp2.SetLineColor(r.kRed)
    gresp2.SetMarkerColor(r.kRed)
    mgesp.Add(gresp1)
    mgesp.Add(gresp2)
    c3 = r.TCanvas()
    c3.SetLogy()
    cSaver.append(c3)
    mgesp.Draw('alp')
    c3.SetLogx()

    ltgr = r.TGraph(len(lt),array('f',coups),array('f',lt))
    c2 = r.TCanvas()
    cSaver.append(c2)
    ltgr.Draw('alp')
    c2.SetLogx()
    c2.SetLogy()
    return accgr, mgesp, ltgr




def makeBRPlot(model):
    pp = physicsParameters()
    mz = pp.masses[pp.name2particle['Z']]
    pp.setNCoupling(model)
    brpie, brpimu, brpinu, brrhoe, brrhomu, brrhonu, br3nu, breenu, brmumunu, brtaus, lt, uCalc, inclusiveHadrBR, DetectableFractionBR, decLength, bBAU, mass2 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    m = np.linspace(1,80,400).tolist()
    m.remove(m[0])
    for mass in m:
        pp.setNMass(mass)
        #brpie.append(pp.Width_H_l('pi',1)/pp.NDecayWidth())
        #brpimu.append(pp.Width_H_l('pi',2)/pp.NDecayWidth())
        #brrhoe.append(pp.Width_H_l('rho',1)/pp.NDecayWidth())
        #brrhomu.append(pp.Width_H_l('rho',2)/pp.NDecayWidth())
        #brpinu.append(sum([pp.Width_H0_nu('pi',l) for l in [1,2,3]])/pp.NDecayWidth())
        #brrhonu.append(sum([pp.Width_H0_nu('rho',l) for l in [1,2,3]])/pp.NDecayWidth())
        #breenu.append(sum([pp.Width_l1_l2_nu(1,1,l) for l in [1,2,3]])/pp.NDecayWidth())
        #brmumunu.append(sum([pp.Width_l1_l2_nu(2,2,l) for l in [1,2,3]])/pp.NDecayWidth())
        #br3nu.append(pp.Width_3nu()/pp.NDecayWidth())
        brpie.append(pp.findBranchingRatio('N -> pi e'))
        brpimu.append(pp.findBranchingRatio('N -> pi mu'))
        brrhoe.append(pp.findBranchingRatio('N -> rho e'))
        brrhomu.append(pp.findBranchingRatio('N -> rho mu'))
        brpinu.append(pp.findBranchingRatio('N -> pi0 nu'))
        brrhonu.append(pp.findBranchingRatio('N -> rho nu'))
        breenu.append(pp.findBranchingRatio('N -> e e nu'))
        brmumunu.append(pp.findBranchingRatio('N -> mu mu nu'))
        br3nu.append(pp.findBranchingRatio('N -> nu nu nu'))
        brtaus.append(pp.findBranchingRatio('N -> mu tau nu') + pp.findBranchingRatio('N -> e tau nu') + pp.findBranchingRatio('N -> tau tau nu'))
        inclusiveHadrBR.append(pp.findBranchingRatio('N -> hadrons'))
        #print mass, (pp.findBranchingRatio('N -> e e nu') + pp.findBranchingRatio('N -> mu mu nu') + pp.findBranchingRatio('N -> e mu nu')
        #     + pp.findBranchingRatio('N -> e tau nu') + pp.findBranchingRatio('N -> mu tau nu')
        #     + pp.findBranchingRatio('N -> nu nu nu') + pp.findBranchingRatio('N -> hadrons'))
        bBAU.append(BAU(mass,pp))
        gamma = (mz/(2.*pp.MN)) + (pp.MN/(2.*mz))
        decLength.append(pp.computeNLifetime()*pp.c*gamma)
        lt.append(pp.computeNLifetime())
        lt1 = pp.computeNLifetime()*pp.U2[1]
        ps = 1. - (pp.MN/90.)**2.
        acc = (45.*pp.c*lt1)
        if acc>1.:
            acc = 1.
        fact = 2.3/(0.16e12)
        #sbr = sum([pp.Width_l1_l2_nu(2,2,l)+pp.Width_l1_l2_nu(1,1,l)+pp.Width_l1_l2_nu(1,2,l)+pp.Width_l1_l2_nu(2,1,l) for l in [1,2,3]])/pp.NDecayWidth() #+ sum([pp.Width_l1_l2_nu(1,1,l) for l in [1,2,3]])/pp.NDecayWidth()
        decList = pp.HNLAllowedDecays()
        DetectableFraction = pp.findBranchingRatio('N -> charged hadrons')
        for dec in decList:
            if decList[dec] == 'yes' and (dec == 'N -> e e nu' or dec == 'N -> mu mu nu' or dec == 'N -> e mu nu'):
                DetectableFraction += pp.findBranchingRatio(dec)
        sbr = DetectableFraction
        DetectableFractionBR.append(DetectableFraction)
        #prova.append(lt1*45.*pp.c/pp.U2[1])
        #print sbr, acc, ps
        uc = math.sqrt(fact*(1./sbr)*acc*(1./ps))
        if (45.*pp.c*lt1/uc) > 5.e-3:
            uCalc.append(uc)
            mass2.append(mass)
        #if mass < 10.:
        #    inclusiveHadrBR.append(sum([pp.Width_H_l('pi',l) + pp.Width_H_l('rho',l) + pp.Width_H0_nu('pi0',l) + pp.Width_H0_nu('rho',l) + pp.Width_H0_nu('eta',l) + pp.Width_H0_nu('eta1',l) for l in [1,2,3]])/pp.NDecayWidth())
        #else:
        #    inclusiveHadrBR.append(sum([pp.Width_l1_l2_nu(a,b,g) for a in range(4,10) for b in range(4,10) for g in [1,2,3]])/pp.NDecayWidth())
    #brpieR = array('f',brpie)
    #brpinuR = array('f',brpinu)
    #brrhoeR = array('f',brrhoe)
    #br3nuR = array('f',br3nu)
    mR = array('f',m)
    mgr = r.TMultiGraph()
    grpie = r.TGraph(len(m),mR,array('f',brpie))
    grpinu = r.TGraph(len(m),mR,array('f',brpinu))
    grrhoe = r.TGraph(len(m),mR,array('f',brrhoe))
    gr3nu = r.TGraph(len(m),mR,array('f',br3nu))
    grrhonu = r.TGraph(len(m),mR,array('f',brrhonu))
    grpimu = r.TGraph(len(m),mR,array('f',brpimu))
    grrhomu = r.TGraph(len(m),mR,array('f',brrhomu))
    greenu = r.TGraph(len(m),mR,array('f',breenu))
    grmumunu = r.TGraph(len(m),mR,array('f',brmumunu))
    grtaus = r.TGraph(len(m),mR,array('f',brtaus))
    grHincl = r.TGraph(len(m),mR,array('f',inclusiveHadrBR))
    grDet = r.TGraph(len(m),mR,array('f',DetectableFractionBR))
    grpie.SetLineColor(r.kBlue)
    grpie.SetLineWidth(3)
    grpinu.SetLineColor(r.kRed)
    grpinu.SetLineWidth(3)
    grrhonu.SetLineColor(r.kYellow)
    grrhonu.SetLineWidth(3)
    grrhomu.SetLineColor(r.kGray)
    grrhomu.SetLineWidth(3)
    grpimu.SetLineColor(r.kOrange)
    grpimu.SetLineWidth(3)
    grrhoe.SetLineColor(r.kGreen)
    grrhoe.SetLineWidth(3)
    gr3nu.SetLineColor(r.kBlack)
    gr3nu.SetLineWidth(3)
    greenu.SetLineColor(r.kViolet)
    greenu.SetLineWidth(3)
    grmumunu.SetLineColor(r.kMagenta)
    grmumunu.SetLineWidth(3)
    grtaus.SetLineColor(r.kMagenta+3)
    grtaus.SetLineWidth(3)
    grHincl.SetLineColor(r.kBlack)
    grHincl.SetLineWidth(3)
    grHincl.SetLineStyle(7)
    grDet.SetLineColor(r.kRed)
    grDet.SetLineWidth(3)
    grDet.SetLineStyle(7)
    mgr.Add(gr3nu)
    mgr.Add(grpie)
    mgr.Add(grpimu)
    mgr.Add(grrhoe)
    mgr.Add(grrhomu)
    mgr.Add(grpinu)
    mgr.Add(grrhonu)
    mgr.Add(greenu)
    mgr.Add(grmumunu)
    mgr.Add(grtaus)
    mgr.Add(grHincl)
    mgr.Add(grDet)
    mgr.SetTitle("Branching ratios for HNL (model: U^{2} = %s)"%str(model))
    c1 = r.TCanvas(str(model),str(model))
    cSaver.append(c1)
    mgr.Draw("alp")
    mgr.GetYaxis().SetRangeUser(0.00001,1.)
    c1.SetLogx()
    c1.SetLogy()
    mgr.GetXaxis().SetRangeUser(m[0],m[-1])
    mgr.GetXaxis().SetTitle("HNL mass (GeV)")
    leg = r.TLegend(0.6,0.1,0.9,0.7)
    leg.SetFillColor(r.kWhite)
    leg.AddEntry(gr3nu,'N#rightarrow#nu#nu#nu','l')
    leg.AddEntry(grpinu,'N#rightarrow#pi#nu','l')
    leg.AddEntry(grpie,'N#rightarrow#pie','l')
    leg.AddEntry(grpimu,'N#rightarrow#pi#mu','l')
    leg.AddEntry(grrhonu,'N#rightarrow#rho#nu','l')
    leg.AddEntry(grrhoe,'N#rightarrow#rhoe','l')
    leg.AddEntry(grrhomu,'N#rightarrow#rho#mu','l')
    leg.AddEntry(greenu,'N#rightarrowe^{+}e^{-}#nu','l')
    leg.AddEntry(grmumunu,'N#rightarrow#mu^{+}#mu^{-}#nu','l')
    leg.AddEntry(grtaus,'N#rightarrow#tau X','l')
    leg.AddEntry(grHincl,'N#rightarrow hadrons','l')
    leg.AddEntry(grDet,'Detectable fraction','l')
    leg.Draw()

    c2 = r.TCanvas("lifetime","lifetime")
    cSaver.append(c2)
    grlt = r.TGraph(len(m),mR,array('f',lt))
    grlt.SetTitle("HNL lifetime for model %s"%str(model))
    grlt.GetXaxis().SetTitle("HNL mass (GeV)")
    grlt.GetYaxis().SetTitle("#tau_{N} (s)")
    grlt.SetLineWidth(3)
    grlt.Draw("ac")
    #c2.SetLogy()
    grlt.GetXaxis().SetRangeUser(m[0],m[-1])

    c3 = r.TCanvas("TLEP sensitivity")
    cSaver.append(c3)
    grtlep = r.TGraph(len(mass2),array('f',mass2),array('f',uCalc))
    #grtlep = r.TGraph(len(m),mR,array('f',uCalc))
    grtlep.SetLineColor(r.kBlue)
    grtlep.SetLineWidth(3)
    grbau = r.TGraph(len(m),mR,array('f',bBAU))
    grbau.SetLineColor(r.kBlack)
    grbau.SetLineWidth(1004)
    grbau.SetFillStyle(3002)
    mmgr = r.TMultiGraph()
    mmgr.Add(grbau)
    mmgr.Add(grtlep)
    mmgr.SetTitle('TLEP sensitivity')
    c3.SetLogy()
    mmgr.Draw("ac")
    mmgr.GetXaxis().SetRangeUser(1.,75.)
    mmgr.GetYaxis().SetRangeUser(1.e-12,5.e-7)
    #mmgr.GetYaxis().SetRangeUser(min([ui for ui in bBAU if ui>0]),max(uCalc))
    #mmgr.GetYaxis().SetRangeUser(min([ui for ui in bBAU if ui>0]),max(uCalc))
    mmgr.GetXaxis().SetTitle("HNL mass (GeV)")
    mmgr.GetYaxis().SetTitle("U_{#mu}^{2}")
    #fBau.Draw('same')

    c4 = r.TCanvas('N decay width')
    cSaver.append(c4)
    grtw = r.TGraph(len(m),mR,array('f',decLength))
    grtw.SetTitle('HNL decay lenght at TLEP')
    grtw.GetXaxis().SetTitle('HNL mass (GeV)')
    grtw.GetYaxis().SetTitle('c#gamma#tau_{N} (m)')
    grtw.Draw('alp')
    c4.SetLogy()
    c4.SetLogx()
    #c4 = r.TCanvas("prova")
    #cSaver.append(c4)
    #grsbr = r.TGraph(len(m),mR,array('f',prova))
    #c4.SetLogy()
    #grsbr.Draw("ac")
    return mgr,leg,grlt,mmgr,grtw#grtlep#,grsbr
