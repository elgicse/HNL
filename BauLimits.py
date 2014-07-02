from __future__ import division
import math
import numpy as np
from array import array
import ROOT as r

cSaver = []
modelNames = [None, 'e', 'mu', 'tau', ''] # last is U2tot
# normal + inverted
logPosRightOfKinks = [None, -0.4, 0., -0.3, -0.1] + [0.2, 0., -0.6, 0.]
logPosLeftOfKinks = [None, -0.6, -0.6, -0.5, -0.25] + [-0.2, 0., -0.78, -0.6]
logPosOfKinks = [None, -0.55, -0.4, -0.4, -0.15] + [0., 0., -0.7, -0.4]

#def frange(x, y, jump):
#  while x < y:
#    yield x
#    x += jump

def importBBN(model, hierarchy, mmin, mmax):
    filename = 'Limits/existing/data/u2%s_%s_bbn.csv'%(modelNames[model], hierarchy)
    mass, u2 = [], []
    with open(filename,'r') as f:
        for line in f:
            line = line.split(', ')
            if mmin < float(line[0]) < mmax:
                mass.append(float(line[0]))
                u2.append(float(line[1]))
    grbbn = r.TGraph(len(mass), array('f',mass), array('f',u2))
    return grbbn

def importSeesaw(model, hierarchy, mmin, mmax, npoints=100):
    filename = 'Limits/existing/data/u2%s_%s_seesaw.csv'%(modelNames[model], hierarchy)
    rawMSSSup, rawUSSSup = [], []
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(', ')
            rawMSSSup.append(float(line[0]))
            rawUSSSup.append(float(line[1]))
    logGr = r.TGraph(len(rawMSSSup))
    for i in xrange(len(rawMSSSup)):
        logGr.SetPoint(i, math.log10(rawMSSSup[i]), math.log10(rawUSSSup[i]))
    fitRes = logGr.Fit('pol1')
    func = r.TF1(logGr.GetFunction('pol1'))
    extrapLogM = np.arange(math.log10(mmin), math.log10(mmax), (math.log10(mmax)-math.log10(mmin))/npoints).tolist()
    extrapLogU = []
    for i in xrange(len(extrapLogM)):
        extrapLogU.append(func.Eval(extrapLogM[i]))
    extrapGr = r.TGraph(len(extrapLogM))
    for i in xrange(len(extrapLogM)):
        extrapGr.SetPoint(i, 10.**extrapLogM[i], 10.**extrapLogU[i])
    return extrapGr, logGr # ritorno anche questo senno' fa segfault

def importBAU(model, hierarchy, mmin, mmax, npoints=100):
    """ Returns lower and upper BAU limits """
    filenameSup = 'Limits/existing/data/u2%s_%s_bau_sup.csv'%(modelNames[model], hierarchy)
    filenameInf = 'Limits/existing/data/u2%s_%s_bau_inf.csv'%(modelNames[model], hierarchy)
    rawMBauSup, rawUBauSup = [], []
    with open(filenameSup, 'r') as f:
        for line in f:
            line = line.split(', ')
            rawMBauSup.append(float(line[0]))
            rawUBauSup.append(float(line[1]))
    #rawGrSup = r.TGraph(len(rawMBauSup))
    logGrSup = r.TGraph(len(rawMBauSup))
    for i in xrange(len(rawMBauSup)):
    #    rawGrSup.SetPoint(i, rawMBauSup[i], rawUBauSup[i])
        logGrSup.SetPoint(i, math.log10(rawMBauSup[i]), math.log10(rawUBauSup[i]))
    rawMBauInf, rawUBauInf = [], []
    with open(filenameInf, 'r') as f:
        for line in f:
            line = line.split(', ')
            rawMBauInf.append(float(line[0]))
            rawUBauInf.append(float(line[1]))
    #rawGrInf = r.TGraph(len(rawMBauInf))
    logGrInf = r.TGraph(len(rawMBauInf))
    for i in xrange(len(rawMBauInf)):
    #    rawGrInf.SetPoint(i, rawMBauInf[i], rawUBauInf[i])
        logGrInf.SetPoint(i, math.log10(rawMBauInf[i]), math.log10(rawUBauInf[i]))
    #mgr = r.TMultiGraph()
    #mgr.Add(rawGrSup)
    #mgr.Add(rawGrInf)
    #c = r.TCanvas()
    #cSaver.append(c)
    #mgr.Draw('alp')
    #c.SetLogx()
    #c.SetLogy()
    #cLogSup = r.TCanvas()
    #cSaver.append(cLogSup)
    #logGrSup.Draw('alp')
    #cLogInf = r.TCanvas()
    #cSaver.append(cLogInf)
    #logGrInf.Draw('alp')
    # Extrapolation
    if hierarchy == 'inverted':
        model = model + 4
    fitResRight = logGrInf.Fit('pol1','','',logPosRightOfKinks[model],1.)
    funcRight = r.TF1(logGrInf.GetFunction('pol1'))
    funcRight.SetLineColor(r.kRed)
    fitResLeft = logGrInf.Fit('pol1','','',-1.,logPosLeftOfKinks[model])
    funcLeft = r.TF1(logGrInf.GetFunction('pol1'))
    funcLeft.SetLineColor(r.kBlue)
    extrapLogMInf = np.arange(math.log10(mmin), math.log10(mmax), (math.log10(mmax)-math.log10(mmin))/npoints).tolist()
    extrapLogUInf = []
    for i in xrange(len(extrapLogMInf)):
        if extrapLogMInf[i] >= logPosOfKinks[model]:
            extrapLogUInf.append(funcRight.Eval(extrapLogMInf[i]))
        else:
            extrapLogUInf.append(funcLeft.Eval(extrapLogMInf[i]))
    # Uncomment and return 'mgr2' to check results
    #logExtrapGrInf = r.TGraph(len(extrapLogMInf), array('f',extrapLogMInf), array('f',extrapLogUInf))
    #logExtrapGrInf.SetLineColor(r.kGreen)
    #mgr2 = r.TMultiGraph()
    #mgr2.Add(logExtrapGrInf)
    #mgr2.Add(logGrInf)
    #newcanv = r.TCanvas()
    #cSaver.append(newcanv)
    #mgr2.Draw('alp')
    extrapGrInf = r.TGraph(len(extrapLogMInf))
    for i in xrange(len(extrapLogMInf)):
        mTemp = 10.**extrapLogMInf[i]
        extrapGrInf.SetPoint(i, mTemp, 10.**extrapLogUInf[i])
    extrapGrInf.SetLineWidth(3)
    extrapGrInf.SetLineColor(r.kOrange)

    fitResSup = logGrSup.Fit('pol1')
    funcSup = r.TF1(logGrSup.GetFunction('pol1'))
    extrapLogUSup = []
    for i in xrange(len(extrapLogMInf)):
        extrapLogUSup.append(funcSup.Eval(extrapLogMInf[i]))
    extrapGrSup = r.TGraph(len(extrapLogMInf))
    for i in xrange(len(extrapLogMInf)):
        mTemp = 10.**extrapLogMInf[i]
        extrapGrSup.SetPoint(i, mTemp, 10.**extrapLogUSup[i]*((1 - (mTemp**2./80.**2.))**2.))
    extrapGrSup.SetLineColor(r.kBlue)

    #cAltra = r.TCanvas()
    #cSaver.append(cAltra)
    #mgr3 = r.TMultiGraph()
    #mgr3.Add(extrapGrInf)
    #mgr3.Add(extrapGrSup)
    #logGrInf.SetLineWidth(3)
    #mgr3.Add(rawGrInf)
    #mgr3.Draw('alp')
    #cAltra.SetLogx()
    #cAltra.SetLogy()
    return extrapGrInf, extrapGrSup, logGrSup, logGrInf # ritorno anche questi senno' fa segfault

