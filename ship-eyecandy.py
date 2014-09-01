import sys
import ROOT as r
import numpy as np
from BauLimits import *
from options import *

if __name__ == '__main__':
    model = int(sys.argv[1])
    f = r.TFile('out/plots/SHiP-graphs-model%s.root'%model,'read')
    f.ls()
    gr = f.Get('Super')
    h = 'normal'
    if model == 1:
        h = 'inverted'
    mmin = gr.GetXaxis().GetXmin()
    mmax = gr.GetXaxis().GetXmax()
    print mmin, mmax
    #grBauLower, grBauUpper, useless1, useless2 = importBAU(4, h, mmin, mmax) #Utot
    grBauLower, grBauUpper, useless1, useless2 = importBAU(4, h, 0.09, 5.5030595699) #Utot

    grBauUpper.SetTitle('BAUupper')
    grBauUpper.SetName('BAUupper')
    grBauUpper.SetLineWidth(5004)
    #grBauUpper.SetFillStyle(3002)
    grBauUpper.SetFillStyle(1001)
    grBauUpper.SetFillColor(r.kGray)
    grBauUpper.SetLineColor(r.kBlack)
    grBauLower.SetTitle('Seesaw')
    grBauLower.SetName('Seesaw')
    grBauLower.SetLineWidth(-6004)
    #grBauLower.SetFillStyle(3002)
    grBauLower.SetFillStyle(1001)
    grBauLower.SetFillColor(r.kGray)
    grBauLower.SetLineColor(r.kBlack)

    #grbbn = importBBN(4,h, mmin, mmax) #Utot
    grbbn = importBBN(4, h, 0.09, 5.5030595699) #Utot
    grbbn.SetTitle('BBN')
    grbbn.SetName('BBN')
    grbbn.SetLineWidth(-6004)
    grbbn.SetLineColor(r.kBlack)
    grbbn.SetLineStyle(7)
    #grbbn.SetFillStyle(3002)
    grbbn.SetFillStyle(1001)
    grbbn.SetFillColor(r.kGray)

    # Convert gr to Utot
    pp = physicsParameters()
    temp_n = gr.GetN()
    temp_x = gr.GetX()
    temp_x.SetSize(temp_n)
    temp_y = gr.GetY()
    temp_y.SetSize(temp_n)
    x = np.array(temp_x, copy=True)
    y = np.array(temp_y, copy=True)
    grUtot = r.TGraph(temp_n)
    for i in xrange(temp_n):
        grUtot.SetPoint(i, temp_x[i], pp.factors[model-1]*temp_y[i])


    grUtot.SetTitle('SHiP')
    grUtot.SetName('SHiP')
    grUtot.SetLineWidth(4)
    grUtot.SetLineColor(r.kBlue)
    grUtot.SetMarkerColor(r.kBlue)

    c1 = r.TCanvas('SHiP 90%% C.L. - model %s'%model,'SHiP 90%% C.L. - model %s'%model,359,295)

    mgr = r.TMultiGraph()
    mgr.Add(grbbn)
    mgr.Add(grBauLower)
    mgr.Add(grBauUpper)
    mgr.Add(grUtot)
    c1.SetLogy()
    c1.SetLogx()
    #c1.SetGrid()
    mgr.SetTitle('SHiP 90%% C.L. - model %s'%model)
    mgr.Draw('alp')

    mgr.GetXaxis().SetTitle("HNL mass (GeV)")
    mgr.GetXaxis().CenterTitle()
    mgr.GetYaxis().SetTitle("HNL coupling to SM U^{2}")
    mgr.GetYaxis().CenterTitle()
    mgr.GetXaxis().SetTitleSize(0.05)
    mgr.GetXaxis().SetTitleOffset(0.8)
    mgr.GetYaxis().SetTitleSize(0.05)
    mgr.GetYaxis().SetTitleOffset(0.95)
    mgr.GetYaxis().SetRangeUser(1.e-11, 1.e-6)

    info = r.TLatex()
    info.SetTextSize(0.05)
    if model == 1:
        info.DrawLatexNDC(0.20, 0.55, '#it{BBN}')
        info.DrawLatexNDC(0.42, 0.21, '#it{Seesaw}')
        info.DrawLatexNDC(0.69, 0.83, '#it{BAU}')
    else:
        info.DrawLatexNDC(0.20, 0.58, '#it{BBN}')
        info.DrawLatexNDC(0.38, 0.17, '#it{Seesaw}')
        info.DrawLatexNDC(0.61, 0.80, '#it{BAU}')

    r.gPad.RedrawAxis('g')