# Usage: ipython -i CMS.py
# bot, top, shipfile, labels, leg = sensitivityScanCMS(2, 0.1, 1., 1.e11, 220, 400)
from __future__ import division
import sys
import timeit
from options import *
from scan import *
from BauLimits import *

cSaver = [] # Keep canvases open

# All these functions are needed to look separately for the upper and lower sensitivity contours...

# Sensitivity scan mass range
start = [math.log10(0.1)]*3 #GeV
stop = [math.log10(43.)]*3 #GeV
# Sensitivity scan U^2 range
yStop = [math.log10(1.e-2)]*3
yStart = [math.log10(1.e-13)]*3

#a = math.log10(0.1) # works for 10^12 Z
a = math.log10(0.1)
aa = math.log10(20.)
e = math.log10(80.)
upper = "upper"
lower = "lower"
#sl1 = (-10.+7.)/(math.log10(70.)-math.log10(21.)) #works for 10^12 Z
sl1 = (-10.-math.log10(0.9e-5))/(math.log10(50.)-math.log10(7.))
#sl2 = (-11.+6.)/(math.log10(60.)-math.log10(0.7)) #works for 10^12 Z
sl2 = (-11.-math.log10(1.2e-7))/(math.log10(20.)-math.log10(1.))
# ystart (upper, lower)
#newy1 = [-6., math.log10(3.e-8)] # works for 10^12 Z
newy1 = [math.log10(3.e-5), math.log10(3.e-9)]
newy2 = [math.log10(5.e-4), math.log10(8.e-6)]

def inBeltSegments((x,y)):
    if (a < x < e) and (beltSegment(lower, newy1, aa, sl1, x) < y < beltSegment(upper, newy1, aa, sl1, x)): # sopra
        return True
    elif (a < x < e) and (beltSegment(lower, newy2, a, sl2, x) < y < beltSegment(upper, newy2, a, sl2, x)): # sotto
        return True
    return False

def inTopHalf(p1,p2,datum):
    if datum[1] > lineForTwoPoints(p1,p2,datum[0]):
        return True
    return False

def make2Countours(data,m):
    p1 = (math.log10(3.), -3.) #works for 1cm
    p2 = (math.log10(30.5), math.log10(1.5e-9)) # works for 1cm
    #p1 = (math.log10(2.), math.log10(1.e-3)) # works for 10cm
    #p2 = (math.log10(10.), math.log10(1.e-7)) # works for 10cm
    xAxis = []
    for datum in data:
        xAxis.append(datum[0])
    xAxis = list(set(xAxis))
    topContour = [datum for datum in data if inTopHalf(p1,p2,datum)]
    botContour = [datum for datum in data if not inTopHalf(p1,p2,datum)]
    print 'len(contours): ', len(topContour), len(botContour)
    # Now for every mass value select only the eps value with N closest to m
    topContour = sorted(topContour, key=lambda x: x[0]) #sorted by mass
    botContour = sorted(botContour, key=lambda x: x[0]) #sorted by mass
    top = []
    bot = []
    right = []
    for x in xAxis:
        tempTop = []
        for datum in topContour:
            if datum[0] == x:
                tempTop.append(datum)
        tempBot = []
        for datum in botContour:
            if datum[0] == x:
                tempBot.append(datum)
        if tempTop:
            bestPointTop=closestCut(tempTop, m)
            if bestPointTop:
                top.append((x, bestPointTop[1], bestPointTop[2]))
        if tempBot:
            bestPointBot=closest(tempBot, m)
            if bestPointBot:
                bot.append((x, bestPointBot[1], bestPointBot[2]))
    top = list(set(top))
    bot = list(set(bot))
    return top, bot

def closestCut(data,m): # Find the point of the parameter space where N(evts) is closest to 2.3 and is > 0
    useful = [dat for dat in data if dat[2]>0]
    if useful:
        result = min(useful, key=lambda x: math.fabs(x[2]-m))
        return result
    else:
        return False
######################################################################


def BAU(m,pp):
    return 2.e-6*(1/m**2.)*((1.-(m**2./pp.masses[pp.name2particle['W']]**2.))**2.)

def makeSensitivityTLEP(existingData, model, ndivx, ndivy, Rmin, Rmax, nW, verbose=0):
    points = [(x,y) for x in np.linspace(start[model-1], stop[model-1], ndivx) for y in np.linspace(yStart[model-1], yStop[model-1], ndivy)]# if inBeltSegments((x,y))]
    data = []
    pp = physicsParameters()
    ep = experimentParams(pp, 'TLEP')
    ep.Rmin = Rmin
    ep.Rmax = Rmax
    ep.nW = nW
    # Next two lines are useless. It is faster if you recompute the data each time.
    print 'Loaded %s previous data points.'%len(existingData)
    outFilePath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model,ep.Rmin,ep.Rmax,ep.nW)
    tic = timeit.default_timer()
    for i,point in enumerate(points):
        mass = roundToN( pow(10.,point[0]),3 )
        eps  = roundToN( pow(10.,point[1]),3 )
        n = roundToN( computeNEventsTLEP(pp, ep, model, mass, eps),3 )
        logmass = roundToN(point[0],3)
        logeps = roundToN(point[1],3)
        datum = ( logmass, logeps, n)
        if verbose:
            if not i%1000:
                print "Point %s of %s: \t log10(mass) %s \t log10(U2) %s \t\t Events: %s"%(i, len(points), logmass, logeps, n)
                gc.collect()
        data.append(datum)
    toc = timeit.default_timer()
    print 'Execution time: %s s'%(toc-tic)
    return data

def computeNEventsTLEP(pp, ep, model, mass, coupling):
    """ Choose model 1, 2 or 3 """
    pp.setNMass(mass)
    if pp.MN > pp.masses[pp.name2particle['Z']]:
        return 0.
    model = model - 1
    if model == 0:
        couplings = [coupling, pp.models[model][1]*coupling, pp.models[model][2]*coupling]
    elif model == 1:
        couplings = [pp.models[model][0]*coupling, coupling, pp.models[model][2]*coupling]
    elif model == 2:
        couplings = [pp.models[model][0]*coupling, pp.models[model][1]*coupling, coupling]
    pp.setNCoupling(couplings)
    decList = pp.HNLAllowedDecays()
    DetectableFraction = pp.findBranchingRatio('N -> TLEP visible') #(2l + nu) and (2jet + l)
    acc = ep.TLEPacceptance(pp,'W')
    NEv = ep.BRWlnu * pp.factors[model] * 2 * ep.nW * pp.U2[model] * acc * DetectableFraction # pp.factors scales U2_dom to U2_tot
    return NEv

def loadDataFileTLEP(model, Rmin, Rmax, nW):
    filepath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model,Rmin,Rmax,nW)
    if not os.path.isfile(filepath):
        return []
    data = []
    with open(filepath,"r") as ifile:
        for line in ifile:
            line = line.split()
            data.append( ( roundToN(float(line[0]),3),
                roundToN(float(line[1]),3),
                roundToN(float(line[-1]),3) ) )
    return data

#if __name__ == '__main__':
def sensitivityScanTLEP(model=2, Rmin=1.e-3, Rmax=1., nW=1.e12, ndivx=200, ndivy=400):
    verbose = True
    existingData = []#loadDataFileTLEP(model, Rmin, Rmax, nW)
    data = makeSensitivityTLEP(existingData, model, ndivx, ndivy, Rmin, Rmax, nW, verbose)
    #existingData = convertToLog(existingData)
    #data.extend(existingData)
    gc.collect()
    data = list(set(data))
    data.sort(key=lambda x: x[1])
    data.sort(key=lambda x: x[0])
    top, bot = make2Countours(data,2.3)
    top.sort(key=lambda x: x[1])
    top.sort(key=lambda x: x[0])
    bot.sort(key=lambda x: x[1])
    bot.sort(key=lambda x: -x[0])
    pp = physicsParameters()
    grRaw = r.TGraph(len(bot)+len(top))
    for i in xrange(len(top)):
        grRaw.SetPoint(i,10.**top[i][0],pp.factors[model-1]*10.**top[i][1])
    for i in xrange(len(bot)):
        grRaw.SetPoint(i+len(top),10.**bot[i][0],pp.factors[model-1]*10.**bot[i][1])
    grRaw.SetLineWidth(4)
    grRaw.SetLineColor(r.kOrange-3)
    grRaw.SetMarkerColor(r.kOrange-3)
    grRaw.SetTitle('TLEP')
    grRaw.SetName('TLEP')

    pp = physicsParameters()

    # Import BAU limits
    # Depends: BauLimits.py
    # Depends: content of subfolder Limits/existing/data
    modelUtot=4 #U^2 tot
    mmin, mmax = 0.2, 85.
    h = 'normal'
    if model == 1:
        h = 'inverted'
    grBauLower, grBauUpper, useless1, useless2 = importBAU(modelUtot, h, mmin, mmax)
    grBauUpper.SetTitle('BAUupper')
    grBauUpper.SetName('BAUupper')
    grBauUpper.SetLineWidth(6004)
    grBauUpper.SetFillStyle(3002)
    grBauUpper.SetLineColor(r.kBlack)
    grBauLower.SetTitle('Seesaw')
    grBauLower.SetName('Seesaw')
    grBauLower.SetLineWidth(-6004)
    grBauLower.SetFillStyle(3002)
    grBauLower.SetLineColor(r.kBlack)

    # Import Seesaw limits
    # Depends: BauLimits.py
    # Depends: content of subfolder Limits/existing/data
    grbbn = importBBN(modelUtot,h, mmin, mmax)
    grbbn.SetTitle('BBN')
    grbbn.SetName('BBN')
    grbbn.SetLineWidth(-6004)
    grbbn.SetFillColor(r.kBlack)
    grbbn.SetLineColor(r.kBlack)
    grbbn.SetLineStyle(7)
    grbbn.SetFillStyle(3002)

    # Import SHiP sensitivity
    # Depends: content of subfolder out/plots
    grShipFile = r.TFile('out/plots/grModel%s.root'%(model), 'read')
    #grShip = r.TGraph(grShipFile.Get('Graph')) # Raw graph
    grShip = r.TGraph(grShipFile.Get('Super')) # Smoothed graph
    grShip.SetTitle('SHiP')
    grShip.SetName('SHiP')
    grShip.SetLineWidth(4)
    grShip.SetLineColor(r.kBlue)
    grShip.SetMarkerColor(r.kBlue)

    # Draw all together
    r.gStyle.SetPadTopMargin(0.10);
    r.gStyle.SetPadRightMargin(0.05);
    r.gStyle.SetPadBottomMargin(0.13);
    r.gStyle.SetPadLeftMargin(0.10);
    c1 = r.TCanvas('c1','c1',852,380)
    cSaver.append(c1)
    mgr = r.TMultiGraph('mgr_%s_%s'%(h,nW),'mgr_%s_%s'%(h,nW))
    mgr.Add(grBauLower)
    mgr.Add(grbbn)
    mgr.Add(grBauUpper)
    mgr.Add(grShip)
    mgr.Add(grRaw)
    c1.SetLogy()
    c1.SetLogx()
    if h == 'inverted':        
        mgr.SetTitle('#bf{Expected sensitivity to HNL (Inverted Hierarchy)}')
    elif h =='normal':
        mgr.SetTitle('#bf{Expected sensitivity to HNL (Normal Hierarchy)}')
    r.gStyle.SetTitleFontSize(0.065)
    mgr.Draw('alp')
    mgr.GetXaxis().SetMoreLogLabels()
    mgr.GetYaxis().SetRangeUser(1.e-12, 2.e-6)
    mgr.GetXaxis().SetLimits(mmin, mmax)
    mgr.GetXaxis().SetTitle('HNL mass (GeV)')
    mgr.GetYaxis().SetTitle('U_{e}^{2} + U_{#mu}^{2} + U_{#tau}^{2}')
    mgr.GetXaxis().CenterTitle()
    mgr.GetYaxis().CenterTitle()
    mgr.GetXaxis().SetTitleSize(0.06)
    mgr.GetXaxis().SetTitleOffset(0.97)
    mgr.GetYaxis().SetTitleSize(0.06)
    mgr.GetYaxis().SetTitleOffset(0.78)
    mgr.GetXaxis().SetLabelSize(0.05)
    mgr.GetXaxis().SetLabelOffset(-0.004)
    mgr.GetYaxis().SetLabelSize(0.05)
    leg = r.TLegend(0.114,0.148,0.330,0.332)
    if nW > 1.e12:
        leg = r.TLegend(0.114,0.148,0.330,0.332)
    leg.SetFillColor(r.kWhite)
    leg.AddEntry(grShip,'SHiP','lp')
    if nW <= 1.e12:
        leg.AddEntry(grRaw,'CMS 10^{11} W','lp')
    else:
        leg.AddEntry(grRaw,'TLEP 10^{13} Z^{0}','lp')
    leg.SetTextSize(0.06)
    leg.Draw()
    labels = r.TLatex()
    labels.SetTextSize(0.06)
    if h == 'inverted':        
        labels.DrawLatexNDC(0.63, 0.70, '#it{BAU}')
        labels.DrawLatexNDC(0.43, 0.21, '#it{Seesaw}')
        labels.DrawLatexNDC(0.12, 0.52, '#it{BBN}')
    elif h == 'normal':
        labels.DrawLatexNDC(0.64, 0.66, '#it{BAU}')
        labels.DrawLatexNDC(0.45, 0.20, '#it{Seesaw}')
        labels.DrawLatexNDC(0.13, 0.55, '#it{BBN}')
    info = r.TLatex()
    info.SetTextSize(0.06)
    info.SetTextColor(r.kOrange-3)
    info.DrawLatexNDC(0.48, 0.45, 'r_{min}=1cm')
    info.DrawLatexNDC(0.48, 0.37, 'r_{max}=1m')
    #if nW <= 1.e12:
    #    info.DrawLatexNDC(0.48, 0.45, 'r_{min}=1mm')
    #    info.DrawLatexNDC(0.48, 0.37, 'r_{max}=1m')
    #if nW > 1.e12:
    #    info.DrawLatexNDC(0.68, 0.45, 'r_{min}=100#mum')
    #    info.DrawLatexNDC(0.68, 0.37, 'r_{max}=5m')
    c1.SetGrid()
    c1.SetTickx()
    c1.SetTicky()
    c1.Update()

    # Save everything in a rootfile
    # Depends: folder out/plots/perNicola
    filecontuttodentro = r.TFile('out/plots/perNicola/CMS-model%s_%s_%sZ_%s_%s.root'%(model,h,nW,Rmin,Rmax),'recreate')
    mgr.Write()
    c1.Write()
    grShip.Write()
    grbbn.Write()
    grBauLower.Write()
    grBauUpper.Write()
    grRaw.Write()
    leg.Write()
    labels.Write()
    info.Write()
    filecontuttodentro.Close()
    return mgr, bot, top, grShipFile, labels, leg


sensitivityScanCMS = sensitivityScanTLEP

# Merge graph for 10^12 and 10^13 Z0 @ TLEP
# (For now only normal hierarchy)
# Depends: plots in out/plots/perNicola
# generated by the above function
def merge1213nh():
    f13 = r.TFile('out/plots/perNicola/CMS-model2_normal_20000000000.0Z_0.01_1.0.root','read')
    f12 = r.TFile('out/plots/perNicola/CMS-model2_normal_20000000000.0Z_0.1_1.0.root','read')
    mgr = f12.Get('mgr_normal_20000000000.0')
    gr13 = r.TGraph(f13.Get('TLEP'))
    gr13.SetLineStyle(9)
    mgr.Add(gr13)
    r.gStyle.SetOptTitle(0)
    r.gStyle.SetPadTopMargin(0.05)
    r.gStyle.SetPadRightMargin(0.05)
    r.gStyle.SetPadBottomMargin(0.11)
    r.gStyle.SetPadLeftMargin(0.11)
    r.gStyle.SetTitleFontSize(0.065)
    c1 = r.TCanvas('merged','merged',712,615)
    cSaver.append(c1)
    c1.SetLogy()
    c1.SetLogx()
    mgr.Draw('alp')
    mgr.GetXaxis().SetLabelSize(0.04)
    mgr.GetYaxis().SetLabelSize(0.04)
    mgr.GetXaxis().SetTitleSize(0.05)
    mgr.GetYaxis().SetTitleSize(0.05)
    mgr.GetYaxis().SetTitleOffset(1.05)
    mgr.GetYaxis().SetRangeUser(1.e-12,1.3e-6)
    mgr.GetXaxis().SetTitleOffset(1.00)
    leg = r.TLegend(0.3347,0.2862,0.9265,0.4617)
    leg.AddEntry(f12.Get('SHiP'),'SHiP','lp')
    leg.AddEntry(f12.Get('TLEP'),'CMS 10^{11} W^{#pm}, 10cm < r < 1m','lp')
    leg.AddEntry(gr13,'CMS 10^{11} W^{#pm}, 1cm < r < 1m','lp')
    leg.SetFillColor(r.kWhite)
    leg.SetTextSize(0.043)
    leg.Draw()
    labels = r.TLatex()
    labels.SetTextSize(0.05)
    labels.DrawLatexNDC(0.67, 0.67, '#it{BAU}')
    labels.DrawLatexNDC(0.42, 0.21, '#it{Seesaw}')
    labels.DrawLatexNDC(0.14, 0.55, '#it{BBN}')
    c1.SetGrid()
    c1.SetTickx()
    c1.SetTicky()
    c1.Update()
    return mgr, leg, labels#, gr12out, g12smoother, gr12
