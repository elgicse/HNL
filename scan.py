from __future__ import division
import os
import struct
from HNLAcceptance import *

ppp = physicsParameters()
mTau = ppp.masses[ppp.name2particle['tau']]
mMu = ppp.masses[ppp.name2particle['mu']]
mDs = ppp.masses[ppp.name2particle['Ds']]
me = ppp.masses[ppp.name2particle['e']]

start = [math.log10(0.1)]*3 #GeV
#start = [math.log10(6.e-03)]*3
stop = [math.log10(mDs-mMu), math.log10(mDs-mMu), math.log10(mTau-mMu)]

yStart = [[math.log10(5.e-06), math.log10(1.e-07)]]*3
yStop = [[math.log10(1.e-09), math.log10(7.e-11)]]*3
print 25*(math.log10(mTau) - math.log10(5.e-03))/(math.log10(mTau) + 1.)

slope = (-6.+9.)/(-1)

def lineForTwoPoints(p1, p2, x):
    x1 = p1[0]
    x2 = p2[0]
    y1 = p1[1]
    y2 = p2[1]
    slop = (y2-y1)/(x2-x1)
    return y1 + slop*(x - x1)

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
    if (start[model] < x < stop[model]) and (lineForTwoPoints((start[model],yStart[model][1]),(stop[model],yStop[model][1]),x) < y < lineForTwoPoints((start[model],yStart[model][0]),(stop[model],yStop[model][0]),x)):
    #if (start[model] < x < stop[model]) and (beltSegment('lower',yStart[model],start[model],slope,x) < y < beltSegment('upper',yStart[model],start[model],slope,x)):
        return True
    return False

def makeSensitivityBelt(existingData, model, ndivx, ndivy, verbose=0):
    points = [(x,y) for x in np.linspace(start[model-1], stop[model-1], ndivx) for y in np.linspace(yStop[model-1][1], yStart[model-1][0], ndivy) if inBelt((x,y),model-1)]
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
            n = roundToN( computeNEvents(model, mass, eps) )
        logmass = roundToN(point[0])
        logeps = roundToN(point[1])
        datum = ( logmass, logeps, n)
        if verbose:
            if not i%50:
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
    right = []
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
    verbose = True
    model = 1
    print "Scanning model %s..."%model
    existingData3 = loadDataFile(model)
    data3 = makeSensitivityBelt(existingData3, model, 102, 100, verbose)
    existingData3 = convertToLog(existingData3)
    data3.extend(existingData3)
    gc.collect()
    data3 = list(set(data3))
    data3.sort(key=lambda x: x[1])
    data3.sort(key=lambda x: x[0])
    bot3 = makeCountours(data3,2.3)
    bot3.sort(key=lambda x: x[1])
    bot3.sort(key=lambda x: x[0])
    num3 = len(bot3)
    gr3 = r.TGraph(num3+1)
    pp = physicsParameters()
    for i in xrange(len(bot3)):
        gr3.SetPoint(i,10.**bot3[i][0],pp.factors[model-1]*10.**bot3[i][1])
    gr3.SetPoint(num3,max([10.**x[0] for x in bot3]),max([10.**x[1] for x in bot3]))

    # qui provo a smoothare
    logGraphTemp = r.TGraph(len(bot3))
    for i in xrange(len(bot3)):
        logGraphTemp.SetPoint(i, bot3[i][0], bot3[i][1])
    logGraphOut1 = r.TGraph(len(bot3))
    logGraphOut2 = r.TGraph(len(bot3))
    logGraphOut3 = r.TGraph(len(bot3))
    logGraphSmoother1 = r.TGraphSmooth()
    logGraphSmoother2 = r.TGraphSmooth()
    logGraphSmoother3 = r.TGraphSmooth()
    logGraphOut1 = logGraphSmoother1.SmoothLowess(logGraphTemp)
    logGraphOut2 = logGraphSmoother2.SmoothKern(logGraphTemp)
    logGraphOut3 = logGraphSmoother3.SmoothSuper(logGraphTemp)
    cSmooth = r.TCanvas()
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
    GraphTemp = r.TGraph(len(bot3)+1)
    for i in xrange(len(bot3)):
        GraphTemp.SetPoint(i, 10.**bot3[i][0], pp.factors[model-1]*10.**bot3[i][1])
    GraphTemp.SetPoint(num3,max([10.**x[0] for x in bot3]),max([10.**x[1] for x in bot3]))
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
    gr3.SetLineWidth(1004)
    gr3.SetFillStyle(3002)
    gr3.SetLineColor(r.kBlue)
    gr3.SetMarkerColor(r.kBlue)
    gr3.SetTitle('HNL model III - ignore BAUs')
    graphFile = r.TFile('out/plots/grModel%s.root'%model, 'recreate')
    gr3.Write()
    GraphNoErrs3.Write()
    graphFile.Close()

    ascisse = [x[0] for x in bot3]
    ordinate = [lineForTwoPoints((start[2],yStart[2][1]),(stop[2],yStop[2][1]),x) for x in ascisse]
    ordinate2 = [lineForTwoPoints((start[2],yStart[2][0]),(stop[2],yStop[2][0]),x) for x in ascisse]
    grline = r.TGraph(len(ascisse),array('d',[10.**x for x in ascisse]),array('d',[10.**y for y in ordinate]))
    grline2 = r.TGraph(len(ascisse),array('d',[10.**x for x in ascisse]),array('d',[10.**y for y in ordinate2]))

    seesaw = readFile('seesaw-tautau-normal.csv')
    grss = r.TGraph(len(seesaw))
    for i in xrange(len(seesaw)):
        grss.SetPoint(i,seesaw[i][0],seesaw[i][1])
    grss.SetTitle('Seesaw')
    grss.SetLineWidth(-4004)
    grss.SetLineStyle(2)
    grss.SetFillStyle(3002)
    grss.SetLineColor(r.kBlack)
    grss.SetMarkerColor(r.kBlack)

    baulow = readFile('bau-low-tautau-normal.csv')
    grbl = r.TGraph(len(baulow))
    for i in xrange(len(baulow)):
        grbl.SetPoint(i,baulow[i][0],baulow[i][1])
    grbl.SetTitle('Seesaw')
    grbl.SetLineWidth(-4004)
    grbl.SetLineStyle(1)
    grbl.SetFillStyle(3002)
    grbl.SetLineColor(r.kBlack)
    grbl.SetMarkerColor(r.kBlack)

    bauhigh = readFile('bau-high-tautau-normal.csv')
    grbh = r.TGraph(len(bauhigh))
    for i in xrange(len(bauhigh)):
        grbh.SetPoint(i,bauhigh[i][0],bauhigh[i][1])
    grbh.SetTitle('Seesaw')
    grbh.SetLineWidth(5004)
    grbh.SetLineStyle(1)
    grbh.SetFillStyle(3002)
    grbh.SetLineColor(r.kBlack)
    grbh.SetMarkerColor(r.kBlack)

    GraphNoErrs3.SetLineWidth(1004)
    GraphNoErrs3.SetFillStyle(3002)
    GraphNoErrs3.SetLineColor(r.kRed)
    GraphNoErrs3.SetMarkerColor(r.kRed)

    c3 = r.TCanvas()
    c3.SetLogx()
    c3.SetLogy()
    gr = r.TMultiGraph()
    gr.Add(grbh)
    gr.Add(grbl)
    gr.Add(grss)
    gr.Add(gr3)
    gr.Add(GraphNoErrs3)
    gr.Add(grline)
    gr.Add(grline2)
    gr.Draw('alp')
    gr.GetXaxis().SetTitle('HNL mass (GeV)')
    gr.GetYaxis().SetTitle('U_{#tau}^{2}')
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetTitleOffset(0.90)
    gr.GetYaxis().SetTitleOffset(0.90)
