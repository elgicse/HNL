from __future__ import division
import sys
import os
import struct
from HNLAcceptance import *

ppp = physicsParameters()
mTau = ppp.masses[ppp.name2particle['tau']]
mMu = ppp.masses[ppp.name2particle['mu']]
mDs = ppp.masses[ppp.name2particle['Ds']]
me = ppp.masses[ppp.name2particle['e']]
mB = ppp.masses['B']

start = [math.log10(0.1)]*3 #GeV
#start = [math.log10(6.e-03)]*3
#stop = [math.log10(mDs-mMu), math.log10(mDs-mMu), math.log10(mTau-mMu)]
stop = [math.log10(mB-me), math.log10(mB-mMu), math.log10(mB-mTau)]

y_start = [math.log10(1.e-12)]*3
y_stop = [math.log10(1.e-4)]*3

#yStart = [[math.log10(5.e-06), math.log10(1.e-07)]]*3 #model 3
#yStop = [[math.log10(1.e-09), math.log10(5.e-11)-1]]*3 #model 3
#yStart = [[math.log10(5.e-06), math.log10(1.e-07)]]*3 #model 1 & 2
#yStop = [[math.log10(1.e-09), math.log10(7.e-11)]]*3 #model 1 & 2
#print 25*(math.log10(mTau) - math.log10(5.e-03))/(math.log10(mTau) + 1.)

slope = (-6.+9.)/(-1)

def lineForTwoPoints(p1, p2, x):
    x1 = p1[0]
    x2 = p2[0]
    y1 = p1[1]
    y2 = p2[1]
    slop = (y2-y1)/(x2-x1)
    return y1 + slop*(x - x1)

def Vline(model,p1,p2,p3,x):
    x1 = p1[model][0]
    x2 = p2[model][0]
    x3 = p3[model][0]
    y1 = p1[model][1]
    y2 = p2[model][1]
    y3 = p3[model][1]
    slop1 = (y2-y1)/(x2-x1)
    slop2 = (y3-y2)/(x3-x2)
    if x < x2:
        y = y1 + slop1*(x - x1)
    else:
        y = y2 + slop2*(x - x2)
    return y


p1high = [(math.log10(0.1), math.log10(7.e-5) )]*3
p1low  = [(math.log10(0.1), math.log10(1.e-7) )]*3
p2high = [(math.log10(1.5), math.log10(3.e-8) )]*2 + [(math.log10(1.5), math.log10(1.e-6))]
p2low  = [(math.log10(1.5), math.log10(2.e-12))]*2 + [(math.log10(1.5), math.log10(8.e-11))]
p3high = [(math.log10(6.),  math.log10(5.e-5) )]*2 + [(math.log10(3.), math.log10(2.e-5))]
p3low  = [(math.log10(6.),  math.log10(1.e-7) )]*2 + [(math.log10(3.), math.log10(7.e-8))]


def roundToN(x, n=2):
    if x:
        result = round(x, -int(math.floor(math.log10(math.fabs(x)))) + (n - 1))
    else:
        result = 0.
    return result

def beltSegment(pos, yi, xi, slop, x):
    if pos == "upper":
        starty = yi[0]
    elif pos == "lower":
        starty = yi[1]
    else:
        print "ERROR: select upper or lower segment!"
    return starty - slop*xi + slop*x

def inBelt((x,y),model):
    #if (start[model] < x < stop[model]) and (lineForTwoPoints((start[model],yStart[model][1]),(stop[model],yStop[model][1]),x) < y < lineForTwoPoints((start[model],yStart[model][0]),(stop[model],yStop[model][0]),x)):
    if (start[model] < x < stop[model]) and (Vline(model-1,p1low,p2low,p3low,x) < y < Vline(model-1,p1high,p2high,p3high,x)):
        return True
    return False

def makeSensitivityBelt(root_dir_path, existingData, model, ndivx, ndivy, verbose=0):
    points = [(x,y) for x in np.linspace(start[model-1], stop[model-1], ndivx) for y in np.linspace(y_start[model-1], y_stop[model-1], ndivy) if inBelt((x,y),model-1)]
    print "Scanning number of events on %s phase-space points"%len(points)
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
            n = roundToN( computeNEvents(model, mass, eps, root_dir_path) )
        logmass = roundToN(point[0])
        logeps = roundToN(point[1])
        datum = ( logmass, logeps, n)
        if verbose:
            if not i%250:
                print "Point %s of %s: \t log10(mass) %s \t log10(U2) %s \t\t Events: %s"%(i, len(points), logmass, logeps, n)
                gc.collect()
        data.append(datum)
        gc.collect()
    return data
    #return points


def closest(data,m):
    useful = [dat for dat in data if (m/2. <= dat[2] <= 2.*m)]
    if useful:
        result = min(useful, key=lambda x: math.fabs(x[2]-m))
        return result
    else:
        return False

def makeCountours(data,m):
    xAxis = []
    for datum in data:
        xAxis.append(datum[0])
    xAxis = list(set(xAxis))
    # Now for every mass value select only the eps value with N closest to m
    botContour = sorted(data, key=lambda x: x[0]) #sorted by mass
    bot = []
    #right = []
    for x in xAxis:
        tempBot = []
        for datum in botContour:
            if datum[0] == x:
                tempBot.append(datum)
        if tempBot:
            bestPointBot=closest(tempBot, m)
            if bestPointBot:
                bot.append((x, bestPointBot[1], bestPointBot[2]))
    bot = list(set(bot))
    return bot

def loadDataFile(model):
    filepath = "out/TextData/sensitivityScan-HNL-model%s.txt"%model
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

def convertToLog(data):
    converted = []
    for datum in data:
        converted.append( ( roundToN(math.log10(datum[0])),
            roundToN(math.log10(datum[1])),
            roundToN(datum[2]) ) )
    return converted

def eq(n1, n2, toll = 1.e-01):
    if math.fabs(n2-n1) < toll*math.fabs(n2):
        return True
    return False

def readFile(filename):
    data = []
    with open('Limits/existing/'+filename,'r') as f:
        for line in f:
            line = line.split(',')
            data.append((float(line[0]), float(line[1])))
    return data



if __name__ == '__main__':
    model = int(sys.argv[1])
    root_dir_path = sys.argv[2]
    print 'Work directory is: %s'%root_dir_path
    verbose = True
    print "Scanning model %s..."%model
    existingData = loadDataFile(model)
    print 'Loaded %s previous data points.'%len(existingData)
    data = makeSensitivityBelt(root_dir_path, existingData, model, 200, 150, verbose)
    #data = makeSensitivityBelt(existingData, model, 4, 4, verbose)
    #data3 = []
    existingData = convertToLog(existingData)
    data.extend(existingData)
    gc.collect()
    data = list(set(data))

    #pp = physicsParameters()
    #rescaledData = []
    #outFilePath = "out/TextData/sensitivityScan-HNL-model%s-rescaled.txt"%(model)
    #for j,datum in enumerate(data3):
    #    if not j%100:
    #        print "Rescaling point %s of %s..."%(j,len(data3))
    #        gc.collect()
    #    mass = 10.**datum[0]
    #    coupling = 10.**datum[1]
    #    pp.setNMass(mass)
    #    if model == 1:
    #        couplings = [coupling, pp.models[model-1][1]*coupling, pp.models[model-1][2]*coupling]
    #    elif model == 2:
    #        couplings = [pp.models[model-1][0]*coupling, coupling, pp.models[model-1][2]*coupling]
    #    elif model == 3:
    #        couplings = [pp.models[model-1][0]*coupling, pp.models[model-1][1]*coupling, coupling]
    #    pp.setNCoupling(couplings)
    #    pp.computeProductionWeights('e')
    #    pp.computeProductionWeights('mu')
    #    adjustedBRs = ( ((pp.nDs + (pp.nD+pp.nD0)*pp.w3body['e']) / pp.nTotCharm) * pp.computeNProdBR(0),
    #        ((pp.nDs + (pp.nD+pp.nD0)*pp.w3body['mu']) / pp.nTotCharm) * pp.computeNProdBR(1),
    #        ((pp.nDs/pp.nTotCharm)*pp.BRDsToTau) * pp.computeNProdBR(2) )
    #    U2tot = sum(pp.U2)
    #    rescaledNevt = (( adjustedBRs[0]*pp.U2[0]/U2tot +
    #        adjustedBRs[1]*pp.U2[1]/U2tot +
    #        2*adjustedBRs[2]*pp.U2[2]/U2tot ) / adjustedBRs[model-1]) * datum[2]  # majorana stuff was not taken into account in model 3!! 30/07/2014
    #    if not j%100:
    #        print "%s -> %s"%(datum[2],rescaledNevt)
    #    #if model == 3:
    #    #    rescaledNevt *= 2 # majorana stuff was not taken into account!! 30/07/2014
    #    rescaledData.append( (roundToN(math.log10(mass)) ,
    #    roundToN(math.log10(coupling)),
    #    rescaledNevt ))
    #    with open(outFilePath,'a') as ofile:
    #        try:
    #            ofile.write('%s \t %s \t %s\n'%(mass, coupling, roundToN(rescaledNevt)))
    #        except KeyboardInterrupt:
    #            pass
    #data3 = rescaledData

    data.sort(key=lambda x: x[1])
    data.sort(key=lambda x: x[0])
    bot = makeCountours(data,2.3)
    bot.sort(key=lambda x: x[1])
    bot.sort(key=lambda x: x[0])
    num = len(bot)
    gr = r.TGraph(num+1)
    pp = physicsParameters()
    for i in xrange(len(bot)):
        gr.SetPoint(i,10.**bot[i][0],pp.factors[model-1]*10.**bot[i][1]) # plot as a function of U^2 tot
    gr.SetPoint(num,max([10.**x[0] for x in bot]),max([10.**x[1] for x in bot])) #close the plot

    # qui provo a smoothare
    logGraphTemp = r.TGraph(len(bot))
    for i in xrange(len(bot)):
        logGraphTemp.SetPoint(i, bot[i][0], bot[i][1])
    logGraphOut1 = r.TGraph(len(bot))
    logGraphOut2 = r.TGraph(len(bot))
    logGraphOut3 = r.TGraph(len(bot))
    logGraphSmoother1 = r.TGraphSmooth()
    logGraphSmoother2 = r.TGraphSmooth()
    logGraphSmoother3 = r.TGraphSmooth()
    logGraphOut1 = logGraphSmoother1.SmoothLowess(logGraphTemp)
    logGraphOut2 = logGraphSmoother2.SmoothKern(logGraphTemp)
    logGraphOut3 = logGraphSmoother3.SmoothSuper(logGraphTemp)
    cSmooth = r.TCanvas('cSmooth', 'cSmooth')
    cSmooth.Divide(2,3)
    cSmooth.cd(1)
    logGraphTemp.SetTitle('Log Original')
    logGraphTemp.Draw('alp')
    cSmooth.cd(2)
    logGraphNoErrs1 = r.TGraph(logGraphOut1.GetN(),logGraphOut1.GetX(),logGraphOut1.GetY())
    logGraphNoErrs1.SetTitle('Log Lowess')
    logGraphNoErrs1.Draw('alp')
    cSmooth.cd(3)
    logGraphNoErrs2 = r.TGraph(logGraphOut2.GetN(),logGraphOut2.GetX(),logGraphOut2.GetY())
    logGraphNoErrs2.SetTitle('Log Kern')
    logGraphNoErrs2.Draw('alp')
    cSmooth.cd(4)
    logGraphNoErrs3 = r.TGraph(logGraphOut3.GetN(),logGraphOut3.GetX(),logGraphOut3.GetY())
    logGraphNoErrs3.SetTitle('Log Super')
    logGraphNoErrs3.Draw('alp')
    cp5 = cSmooth.cd(5)
    GraphTemp = r.TGraph(len(bot)+1)
    for i in xrange(len(bot)):
        GraphTemp.SetPoint(i, 10.**bot[i][0], pp.factors[model-1]*10.**bot[i][1])
    GraphTemp.SetPoint(num,max([10.**x[0] for x in bot]),max([10.**x[1] for x in bot]))
    cp5.SetLogx()
    cp5.SetLogy()
    GraphTemp.SetTitle('Original')
    GraphTemp.Draw('alp')
    cp6 = cSmooth.cd(6)
    npGraph3 = logGraphOut3.GetN()
    GraphNoErrs3 = r.TGraph(npGraph3+1)
    tempx_buff = logGraphOut3.GetX()
    tempx_buff.SetSize(npGraph3)
    tempx = np.array(tempx_buff, copy=True)
    tempy_buff = logGraphOut3.GetY()
    tempy_buff.SetSize(npGraph3)
    tempy = np.array(tempy_buff, copy=True)
    gc.collect()
    for i in xrange(npGraph3):
        GraphNoErrs3.SetPoint(i, 10.**tempx[i], pp.factors[model-1]*10.**tempy[i])
    GraphNoErrs3.SetPoint(npGraph3, 10.**max([xi for xi in tempx]), 10.**max([yi for yi in tempy]))
    cp6.SetLogx()
    cp6.SetLogy()
    GraphNoErrs3.SetTitle('Super')
    GraphNoErrs3.SetName('Super')
    GraphNoErrs3.Draw('alp')




    #gr3 = r.TGraph(num3+2)
    #for i in xrange(len(bot3)):
    #    gr3.SetPoint(i,10.**bot3[i][0],10.**bot3[i][1])
    #gr3.SetPoint(num3,mTau,2.e-10)
    #gr3.SetPoint(num3+1,mTau,max([10.**x[1] for x in bot3]))
    gr.SetLineWidth(1004)
    gr.SetFillStyle(3002)
    gr.SetLineColor(r.kBlue)
    gr.SetMarkerColor(r.kBlue)
    gr.SetTitle('HNL model %s - ignore BAUs'%model)

    ascisse = [x[0] for x in bot]
    ordinate = [Vline(model-1,p1low,p2low,p3low,x) for x in ascisse]#[lineForTwoPoints((start[2],yStart[2][1]),(stop[2],yStop[2][1]),x) for x in ascisse]
    ordinate2 = [Vline(model-1,p1high,p2high,p3high,x) for x in ascisse]#[lineForTwoPoints((start[2],yStart[2][0]),(stop[2],yStop[2][0]),x) for x in ascisse]
    grline = r.TGraph(len(ascisse),array('d',[10.**x for x in ascisse]),array('d',[10.**y for y in ordinate]))
    grline2 = r.TGraph(len(ascisse),array('d',[10.**x for x in ascisse]),array('d',[10.**y for y in ordinate2]))

    if model == 1:
        baustring = 'ee-inverted'
    elif model == 2:
        baustring = 'mumu-normal'
    elif model == 3:
        baustring = 'tautau-normal'

    seesaw = readFile('seesaw-%s.csv'%baustring)
    grss = r.TGraph(len(seesaw))
    for i in xrange(len(seesaw)):
        grss.SetPoint(i,seesaw[i][0],seesaw[i][1])
    grss.SetTitle('Seesaw')
    grss.SetLineWidth(-4004)
    grss.SetLineStyle(2)
    grss.SetFillStyle(3002)
    grss.SetLineColor(r.kBlack)
    grss.SetMarkerColor(r.kBlack)

    baulow = readFile('bau-low-%s.csv'%baustring)
    grbl = r.TGraph(len(baulow))
    for i in xrange(len(baulow)):
        grbl.SetPoint(i,baulow[i][0],baulow[i][1])
    grbl.SetTitle('BAU low')
    grbl.SetLineWidth(-4004)
    grbl.SetLineStyle(1)
    grbl.SetFillStyle(3002)
    grbl.SetLineColor(r.kBlack)
    grbl.SetMarkerColor(r.kBlack)

    bauhigh = readFile('bau-high-%s.csv'%baustring)
    grbh = r.TGraph(len(bauhigh))
    for i in xrange(len(bauhigh)):
        grbh.SetPoint(i,bauhigh[i][0],bauhigh[i][1])
    grbh.SetTitle('BAU high')
    grbh.SetLineWidth(5004)
    grbh.SetLineStyle(1)
    grbh.SetFillStyle(3002)
    grbh.SetLineColor(r.kBlack)
    grbh.SetMarkerColor(r.kBlack)

    GraphNoErrs3.SetLineWidth(1004)
    GraphNoErrs3.SetFillStyle(3002)
    GraphNoErrs3.SetLineColor(r.kRed)
    GraphNoErrs3.SetMarkerColor(r.kRed)

    c3 = r.TCanvas('SHiP sensitivity (model %s)'%model, 'SHiP sensitivity (model %s)'%model)
    c3.SetLogx()
    c3.SetLogy()
    mgr = r.TMultiGraph('mgr-model%s'%model, 'mgr-model%s'%model)
    mgr.Add(grbh)
    mgr.Add(grbl)
    mgr.Add(grss)
    mgr.Add(gr)
    mgr.Add(GraphNoErrs3)
    mgr.Add(grline)
    mgr.Add(grline2)
    mgr.Draw('alp')
    mgr.GetXaxis().SetTitle('HNL mass (GeV)')
    mgr.GetYaxis().SetTitle('U_{tot}^{2}')
    mgr.GetXaxis().SetTitleSize(0.05)
    mgr.GetYaxis().SetTitleSize(0.05)
    mgr.GetXaxis().SetTitleOffset(0.90)
    mgr.GetYaxis().SetTitleOffset(0.90)
    
    graphFile = r.TFile('out/plots/SHiP-graphs-model%s.root'%model, 'recreate')
    gr.Write()
    cSmooth.Write()
    GraphTemp.Write()
    logGraphTemp.Write()
    logGraphNoErrs1.Write()
    logGraphNoErrs2.Write()
    logGraphNoErrs3.Write()
    GraphNoErrs3.Write()
    grbh.Write()
    grbl.Write()
    grss.Write()
    grline.Write()
    grline2.Write()
    mgr.Write()
    c3.Write()
    graphFile.Close()
