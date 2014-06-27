from __future__ import division
import sys
import timeit
from options import *
from scan import *
from BauLimits import *

cSaver = []
# Sensitivity scan mass range
start = [math.log10(0.1)]*3 #GeV
stop = [math.log10(75.)]*3 #GeV
# Sensitivity scan U^2 range
#yStart = [math.log10(5.e-06)]*3
yStart = [math.log10(40.)]*3
yStop = [math.log10(1.e-12)]*3

################### there is also an upper bound... ##################
a = math.log10(0.1)
aa = math.log10(20.)
e = math.log10(80.)
upper = "upper"
lower = "lower"
sl1 = (-10.+7.)/(math.log10(70.)-math.log10(21.))
sl2 = (-11.+6.)/(math.log10(60.)-math.log10(0.7))
# ystart (upper, lower)
newy1 = [-6., math.log10(3.e-8)]
newy2 = [math.log10(4.e-4), math.log10(4.e-7)]
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
    p1 = (math.log10(5.), -6.)
    p2 = (math.log10(80.), math.log10(3.e-11))
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
def closestCut(data,m):
    useful = [dat for dat in data if dat[2]>0]
    if useful:
        result = min(useful, key=lambda x: math.fabs(x[2]-m))
        return result
    else:
        return False
######################################################################

def BAU(m,pp):
    return 2.e-6*(1/m**2.)*((1.-(m**2./pp.masses[pp.name2particle['W']]**2.))**2.)

def makeSensitivityTLEP(existingData, model, ndivx, ndivy, Rmin, Rmax, nZ, verbose=0):
    points = [(x,y) for x in np.linspace(start[model-1], stop[model-1], ndivx) for y in np.linspace(yStart[model-1], yStop[model-1], ndivy) if inBeltSegments((x,y))]
    data = []
    pp = physicsParameters()
    ep = experimentParams(pp, 'TLEP')
    ep.Rmin = Rmin
    ep.Rmax = Rmax
    ep.nZ = nZ
    print 'Loaded %s previous data points.'%len(existingData)
    outFilePath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model,ep.Rmin,ep.Rmax,ep.nZ)
    tic = timeit.default_timer()
    for i,point in enumerate(points):
    #    found = False
        mass = roundToN( pow(10.,point[0]),3 )
        eps  = roundToN( pow(10.,point[1]),3 )
    #    for oldDatum in existingData:
    #        if eq(mass, oldDatum[0]) and eq(eps, oldDatum[1]):
    #            found = True
    #            n = oldDatum[2]
    #            break
    #    if not found:
    #    # indent next line
        n = roundToN( computeNEventsTLEP(outFilePath, pp, ep, model, mass, eps),3 )
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

def computeNEventsTLEP(outFilePath, pp, ep, model, mass, coupling):
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
    DetectableFraction = pp.findBranchingRatio('N -> charged hadrons')
    for dec in decList:
        if decList[dec] == 'yes' and (dec == 'N -> e e nu' or dec == 'N -> mu mu nu' or dec == 'N -> e mu nu'):
            DetectableFraction += pp.findBranchingRatio(dec)
    acc = ep.TLEPacceptance(pp)
    #print 'Acceptance: ', acc, pp.computeNLifetime(), ep.Rmin, DetectableFraction
    NEv = ep.BRZnunu * pp.factors[model] * 2 * ep.nZ * pp.U2[model] * acc * DetectableFraction # moltiplicare anche U2 per factor (U2 -> U2tot)???
    #with open(outFilePath,"a") as ofile:
    #    try:
    #        ofile.write('%s \t %s \t %s \t %s \t %s\n'%(mass, coupling,
    #            DetectableFraction, acc, NEv))
    #    except KeyboardInterrupt:
    #        pass
    #        #ofile.close()
    #        #sys.exit(-1)
    return NEv

def loadDataFileTLEP(model, Rmin, Rmax, nZ):
    filepath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model,Rmin,Rmax,nZ)
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
def sensitivityScanTLEP(model=2, Rmin=1.e-3, Rmax=1., nZ=1.e12, ndivx=200, ndivy=400):
    verbose = True
    existingData = []#loadDataFileTLEP(model, Rmin, Rmax, nZ)
    data = makeSensitivityTLEP(existingData, model, ndivx, ndivy, Rmin, Rmax, nZ, verbose)
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
    gr = r.TGraph(len(bot)+len(top))
    for i in xrange(len(top)):
        gr.SetPoint(i,10.**top[i][0],pp.factors[model-1]*10.**top[i][1])
    for i in xrange(len(bot)):
        gr.SetPoint(i+len(top),10.**bot[i][0],pp.factors[model-1]*10.**bot[i][1])
    gr.SetLineWidth(4)
    #gr.SetLineWidth(-1004)
    #gr.SetFillStyle(3002)
    gr.SetLineColor(r.kOrange-3)
    gr.SetMarkerColor(r.kOrange-3)
    gr.SetTitle('HNL model II')
    pp = physicsParameters()

    modelUtot=4 #U^2 tot
    mmin, mmax = 0.2, 83.
    h = 'normal'
    #grBauLower, grBauUpper, useless1, useless2 = importBAU(model, 'normal', min([10.**el[0] for el in bot]), max([10.**el[0] for el in bot]))
    grBauLower, grBauUpper, useless1, useless2 = importBAU(modelUtot, h, mmin, mmax)
    grBauUpper.SetLineWidth(6004)
    grBauUpper.SetFillStyle(3002)
    grBauUpper.SetLineColor(r.kBlack)
    grBauLower.SetLineWidth(-6004)
    grBauLower.SetFillStyle(3002)
    grBauLower.SetLineColor(r.kBlack)

    grbbn = importBBN(modelUtot,h, mmin, mmax)
    grbbn.SetLineWidth(-6004)
    grbbn.SetFillColor(r.kBlack)
    grbbn.SetLineColor(r.kBlack)
    grbbn.SetLineStyle(7)
    grbbn.SetFillStyle(3002)

    grShipFile = r.TFile('out/plots/grModel%s.root'%(model), 'read')
    grShip = r.TGraph(grShipFile.Get('Graph'))
    grShip.SetLineWidth(4)
    grShip.SetLineColor(r.kBlue)
    grShip.SetMarkerColor(r.kBlue)

    c1 = r.TCanvas()
    cSaver.append(c1)
    mgr = r.TMultiGraph()
    mgr.Add(grBauLower)
    mgr.Add(grbbn)
    mgr.Add(grBauUpper)
    mgr.Add(grShip)
    mgr.Add(gr)
    mgr.Draw('alp')
    mgr.GetXaxis().SetMoreLogLabels()
    mgr.GetYaxis().SetRangeUser(1.e-12, 2.e-6)
    mgr.GetXaxis().SetRangeUser(mmin, mmax)
    c1.SetLogy()
    c1.SetLogx()
    mgr.GetXaxis().SetTitle('HNL mass (GeV)')
    mgr.GetYaxis().SetTitle('U_{e}^{2} + U_{#mu}^{2} + U_{#tau}^{2}')
    c1.SetGrid()
    c1.SetTickx()
    c1.SetTicky()
    c1.Update()
    return mgr, bot, top, grShipFile

