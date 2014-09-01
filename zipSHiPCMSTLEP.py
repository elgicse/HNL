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
    mmin, mmax = 0.2, 84.
    grBauLower, grBauUpper, useless1, useless2 = importBAU(4, h, mmin, mmax) #Utot

    grBauUpper.SetTitle('BAUupper')
    grBauUpper.SetName('BAUupper')
    grBauUpper.SetLineWidth(5004)
    grBauUpper.SetFillStyle(3002)
    grBauUpper.SetFillColor(r.kBlack)
    grBauUpper.SetLineColor(r.kBlack)
    grBauLower.SetTitle('Seesaw')
    grBauLower.SetName('Seesaw')
    grBauLower.SetLineWidth(-6004)
    grBauLower.SetFillStyle(3002)
    grBauLower.SetFillColor(r.kBlack)
    grBauLower.SetLineColor(r.kBlack)

    grbbn = importBBN(4, h, mmin, mmax) #Utot
    grbbn.SetTitle('BBN')
    grbbn.SetName('BBN')
    grbbn.SetLineWidth(-6004)
    grbbn.SetLineColor(r.kBlack)
    grbbn.SetLineStyle(7)
    grbbn.SetFillStyle(3002)
    grbbn.SetFillColor(r.kBlack)

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

    # Now import CMS
    CMScolor = r.kGreen+1
    fC13 = r.TFile('out/plots/perNicola/CMS-SHiPupdated-model2_normal_1e+11W_0.01_1.0.root','read')
    fC12 = r.TFile('out/plots/perNicola/CMS-SHiPupdated-model2_normal_1e+11W_0.1_1.0.root','read')
    grC12 = r.TGraph(fC12.Get('TLEP'))
    grC13 = r.TGraph(fC13.Get('TLEP'))
    grC13.SetLineStyle(7)
    grC12.SetLineColor(CMScolor)
    grC13.SetLineColor(CMScolor)
    grC12.SetMarkerColor(CMScolor)
    grC13.SetMarkerColor(CMScolor)

    # Now import TLEP
    TLEPcolor = r.kOrange-3
    fT13 = r.TFile('out/plots/perNicola/SHiPupdated-model2_normal_1e+13Z_0.0001_5.0.root','read')
    fT12 = r.TFile('out/plots/perNicola/SHiPupdated-model2_normal_1e+12Z_0.001_1.0.root','read')
    grT12 = r.TGraph(fT12.Get('TLEP'))
    grT13 = r.TGraph(fT13.Get('TLEP'))
    grT13.SetLineStyle(7)
    grT12.SetLineColor(TLEPcolor)
    grT13.SetLineColor(TLEPcolor)
    grT12.SetMarkerColor(TLEPcolor)
    grT13.SetMarkerColor(TLEPcolor)

    r.gStyle.SetOptTitle(0)
    r.gStyle.SetPadTopMargin(0.05)
    r.gStyle.SetPadRightMargin(0.30)
    r.gStyle.SetPadBottomMargin(0.13)
    r.gStyle.SetPadLeftMargin(0.11)
    r.gStyle.SetTitleFontSize(0.065)
    c1 = r.TCanvas('merged','merged',1134,554)
    c1.SetLogy()
    c1.SetLogx()
    c1.SetTickx()
    c1.SetTicky()
    c1.SetGridx()
    c1.SetGridy()

    mgr = r.TMultiGraph()
    mgr.Add(grbbn)
    mgr.Add(grBauLower)
    mgr.Add(grBauUpper)
    mgr.Add(grC12)
    mgr.Add(grC13)
    mgr.Add(grT12)
    mgr.Add(grT13)
    mgr.Add(grUtot)
    mgr.Draw('al')
    mgr.SetMinimum(1.e-12)
    mgr.SetMaximum(6.e-6)
    mgr.GetXaxis().SetLimits(mmin, mmax)
    mgr.GetXaxis().SetTitle('HNL mass (GeV)')
    mgr.GetYaxis().SetTitle('U_{e}^{2} + U_{#mu}^{2} + U_{#tau}^{2}')
    mgr.GetXaxis().CenterTitle()
    mgr.GetYaxis().CenterTitle()
    mgr.GetXaxis().SetTitleSize(0.06)
    mgr.GetXaxis().SetTitleOffset(1.0)
    mgr.GetYaxis().SetTitleSize(0.06)
    mgr.GetYaxis().SetTitleOffset(0.78)
    mgr.GetXaxis().SetLabelSize(0.05)
    mgr.GetXaxis().SetLabelOffset(-0.004)
    mgr.GetYaxis().SetLabelSize(0.05)
    mgr.GetXaxis().SetMoreLogLabels()

    leg = r.TLegend(0.7257,0.1300,0.9929,0.9506)
    leg.SetFillColor(r.kWhite)
    leg.SetTextSize(0.055)
    leg.SetTextAlign(13)
    legship = leg.AddEntry(grUtot,'#splitline{SHiP}{}','lp')
    legship.SetTextAlign(12)
    leg.AddEntry(grC12,'#splitline{CMS 10^{11} W^{#pm}}{10cm < r < 1m}','lp')
    leg.AddEntry(grC13,'#splitline{CMS 10^{11} W^{#pm}}{1cm < r < 1m}','lp')
    leg.AddEntry(grT12,'#splitline{TLEP 10^{12} Z^{0}}{1mm < r < 1m}','lp')
    leg.AddEntry(grT13,'#splitline{TLEP 10^{13} Z^{0}}{100#mum < r < 5m}','lp')
    leg.Draw()

    labels = r.TLatex()
    labels.SetTextSize(0.06)
    if model == 1:
        pass
    else:
        labels.DrawLatexNDC(0.46, 0.65, '#it{BAU}')
        labels.DrawLatexNDC(0.32, 0.19, '#it{Seesaw}')
        labels.DrawLatexNDC(0.12, 0.51, '#it{BBN}')
    
    r.gPad.RedrawAxis('g')