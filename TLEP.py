from __future__ import division
#from HNLAcceptance import *
from options import *
# mettere HNLAllowedDecays in options.py e includere solo quello!!!!
from scan import *
from BauLimits import *

cSaver = []
# Sensitivity scan mass range
start = [math.log10(0.1)]*3 #GeV
#start = [math.log10(3.)]*3 #GeV
stop = [math.log10(75.)]*3 #GeV
# Sensitivity scan U^2 range
yStart = [math.log10(5.e-06)]*3
yStop = [math.log10(1.e-12)]*3

def BAU(m,pp):
    return 2.e-6*(1/m**2.)*((1.-(m**2./pp.masses[pp.name2particle['W']]**2.))**2.)

def makeSensitivityTLEP(existingData, model, ndivx, ndivy, Rmin, Rmax, nZ, verbose=0):
    points = [(x,y) for x in np.linspace(start[model-1], stop[model-1], ndivx) for y in np.linspace(yStart[model-1], yStop[model-1], ndivy)]
    data = []
    for i,point in enumerate(points):
        found = False
        mass = roundToN( pow(10.,point[0]) )
        eps  = roundToN( pow(10.,point[1]) )
        for oldDatum in existingData:
            if eq(mass, oldDatum[0]) and eq(eps, oldDatum[1]):
                found = True
                n = oldDatum[2]
                break
        if not found:
            n = roundToN( computeNEventsTLEP(model, mass, eps, Rmin, Rmax, nZ) )
        logmass = roundToN(point[0])
        logeps = roundToN(point[1])
        datum = ( logmass, logeps, n)
        if verbose:
            if not i%1000:
                print "Point %s of %s: \t log10(mass) %s \t log10(U2) %s \t\t Events: %s"%(i, len(points), logmass, logeps, n)
                gc.collect()
        data.append(datum)
        #gc.collect()
    return data

def computeNEventsTLEP(model, mass, coupling, Rmin, Rmax, nZ):
    """ Choose model 1, 2 or 3 """
    pp = physicsParameters()
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
    #print couplings
    pp.setNCoupling(couplings)
    ep = experimentParams(pp, 'TLEP')
    ep.Rmin = Rmin
    ep.Rmax = Rmax
    ep.nZ = nZ
    #decList = HNLAllowedDecays(pp)
    decList = pp.HNLAllowedDecays()
    DetectableFraction = pp.findBranchingRatio('N -> charged hadrons')
    for dec in decList:
        if decList[dec] == 'yes' and (dec == 'N -> e e nu' or dec == 'N -> mu mu nu' or dec == 'N -> e mu nu'):
            DetectableFraction += pp.findBranchingRatio(dec)
    acc = ep.TLEPacceptance(pp)
    #print ep.BRZnunu, ep.nZ, pp.U2[model], acc, DetectableFraction
    NEv = ep.BRZnunu * pp.factors[model] * 2 * ep.nZ * pp.U2[model] * acc * DetectableFraction # moltiplicare anche U2 per factor (U2 -> U2tot)???
    outFilePath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model+1,ep.Rmin,ep.Rmax,ep.nZ)
    with open(outFilePath,"a") as ofile:
        try:
            ofile.write('%s \t %s \t %s \t %s \t %s\n'%(mass, coupling,
                DetectableFraction, acc, NEv))
        except KeyboardInterrupt:
            pass
    return NEv

def loadDataFileTLEP(model, Rmin, Rmax, nZ):
    filepath = "out/TextData/sensitivityScan-HNLatTLEP-model%s-%s-%s-%s.txt"%(model+1,Rmin,Rmax,nZ)
    if not os.path.isfile(filepath):
        return []
    data = []
    with open(filepath,"r") as ifile:
        for line in ifile:
            line = line.split()
            data.append( ( roundToN(float(line[0])),
                roundToN(float(line[1])),
                roundToN(float(line[-1])) ) )
    return data

#if __name__ == '__main__':
def sensitivityScanTLEP(model=2, Rmin=1.e-3, Rmax=1., nZ=1.e12):
    verbose = True
    #model = 2
    existingData = loadDataFileTLEP(model, Rmin, Rmax, nZ)
    #data = makeSensitivityTLEP(existingData, model, 411, 200, verbose)
    data = makeSensitivityTLEP(existingData, model, 200, 400, Rmin, Rmax, nZ, verbose)
    existingData = convertToLog(existingData)
    data.extend(existingData)
    gc.collect()
    data = list(set(data))
    data.sort(key=lambda x: x[1])
    data.sort(key=lambda x: x[0])
    bot = makeCountours(data,2.3)
    bot.sort(key=lambda x: x[1])
    bot.sort(key=lambda x: x[0])
    pp = physicsParameters()
    gr = r.TGraph(len(bot))
    for i in xrange(len(bot)):
        gr.SetPoint(i,10.**bot[i][0],pp.factors[model-1]*10.**bot[i][1])
    gr.SetLineWidth(1004)
    gr.SetFillStyle(3002)
    gr.SetLineColor(r.kBlue)
    gr.SetMarkerColor(r.kBlue)
    gr.SetTitle('HNL model II')
    pp = physicsParameters()
    #bauUpper = [BAU(10.**data[0],pp) for data in bot]
    #grBauUpper = r.TGraph(len(bot))
    #for i in xrange(len(bot)):
    #    grBauUpper.SetPoint(i,10.**bot[i][0],bauUpper[i])
    #grBauUpper.SetLineWidth(1004)
    #grBauUpper.SetFillStyle(3002)
    #grBauUpper.SetLineColor(r.kBlack)

    model=4 #U^2 tot

    grBauLower, grBauUpper, useless1, useless2 = importBAU(model, 'normal', min([10.**el[0] for el in bot]), max([10.**el[0] for el in bot]))
    grBauUpper.SetLineWidth(1004)
    grBauUpper.SetFillStyle(3002)
    grBauUpper.SetLineColor(r.kBlack)
    grBauLower.SetLineWidth(-1004)
    grBauLower.SetFillStyle(3002)
    grBauLower.SetLineColor(r.kBlack)
    c1 = r.TCanvas()
    cSaver.append(c1)
    mgr = r.TMultiGraph()
    mgr.Add(grBauLower)
    mgr.Add(grBauUpper)
    mgr.Add(gr)
    mgr.Draw('alp')
    c1.SetLogy()
    c1.SetLogx()
    mgr.GetXaxis().SetMoreLogLabels()
    mgr.GetXaxis().SetTitle('HNL mass (GeV)')
    mgr.GetYaxis().SetTitle('U_{e}^{2} + U_{#mu}^{2} + U_{#tau}^{2}')
    return mgr

