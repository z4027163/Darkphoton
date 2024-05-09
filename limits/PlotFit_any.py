import os
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")
                                                      
file = "2.012_463"
#file = "4.201_611"
mass = file[:5]
year = '2018'
lumi=0
if year =='2017':
    lumi = 35.3
else:
    lumi = 61.3

os.system("combine -M MultiDimFit -d output/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
os.system("combine -M MultiDimFit -d output/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")



f = TFile("higgsCombineSB.MultiDimFit.mH"+mass+".root","READ")
fb = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+".root","READ")

w = f.Get("w")
w.loadSnapshot("MultiDimFit")

wb = fb.Get("w")
wb.loadSnapshot("MultiDimFit")

m2mu = w.var("m2mu")
m2muvars = RooArgSet(m2mu)
xmin = m2mu.getMin()
xmax = m2mu.getMax()
nbins = m2mu.getBinning().numBins()

data = w.data("data_obs")
data = data.reduce(RooFit.SelectVars(m2muvars))  
rdh = data.binnedClone()
testhisto = rdh.createHistogram("testhisto", m2mu, RooFit.Binning(nbins, xmin, xmax));
#testhisto2 = rdh.createHistogram("testhisto2", m2mu, RooFit.Binning(nbins, xmin, xmax));
ymin=testhisto.GetMinimum()
ymax=testhisto.GetMaximum()


gStyle.SetOptTitle(0); 
gStyle.SetOptStat(0);

gStyle.SetPadTopMargin(0.08);
gStyle.SetPadBottomMargin(0.39);
gStyle.SetPadLeftMargin(0.15);
gStyle.SetPadRightMargin(0.07);

gStyle.SetNdivisions(508, "X");
gStyle.SetNdivisions(508, "Y");

c1 = TCanvas("c1", "c1", 650, 720);
c1.SetLogx(0);
c1.SetLogy(0);
c1.SetBottomMargin(0.29);
c1.SetFillStyle(4000)

hframe = TH2F("hframe","hframe",500,xmin,xmax,500,0.8*ymin,1.2*ymax);
hframe.SetXTitle("m_{#mu#mu} [GeV]");

hframe.GetXaxis().SetLabelSize(18);
hframe.GetXaxis().SetLabelFont(43);
hframe.GetXaxis().SetTitleFont(63);
hframe.GetXaxis().SetTitleSize(22);
hframe.GetXaxis().SetTitleOffset(5.0);                                                                                                                                                                                                     
hframe.GetXaxis().SetLabelOffset(0.215);         
hframe.GetXaxis().SetNdivisions(505);


hframe.GetYaxis().SetLabelSize(20);
hframe.GetYaxis().SetLabelFont(43);
hframe.GetYaxis().SetTitleFont(63);
hframe.GetYaxis().SetTitleSize(22);
hframe.GetYaxis().SetNdivisions(505);
hframe.GetYaxis().SetTitleOffset(1.4);
hframe.GetYaxis().SetMaxDigits(2)
hframe.GetYaxis().SetTitle("Events / 0.003 GeV")
#hframe.GetYaxis().SetLabelOffset(0.007);

      
frame = m2mu.frame();
data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(1));
totalpdf = w.pdf("pdf_binCatAB_nuis")
totalpdf.plotOn(frame, RooFit.LineColor(kOrange-3), RooFit.Name("total_pdf"))


bkgpdf = wb.pdf("pdf_binCatAB_nuis")
#bkgpdf.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"));

bkgpdf1 = w.pdf("shapeBkg_bkg_mass_CatAB")
bkgpdf1.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"));

frame2=m2mu.frame()
sigpdf = w.pdf("shapeSig_signalModel_generic_CatAB")
nsig = w.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
#sigpdf.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig,RooAbsReal.NumEvent))
print("number of events under sig: " + str(nsig))


testpdf   = bkgpdf1 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf.Scale(testhisto.Integral(1,-1)/testpdf.Integral(1,-1));

sbpdf   = totalpdf .createHistogram("totalpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));

sigLegend = TH2F();
sigLegend.SetLineColor(kOrange-3)
sigLegend.SetLineWidth(3)

bkgLegend = TH2F();
bkgLegend.SetLineColor(4)
bkgLegend.SetLineWidth(3)

datLegend = TH2F()
datLegend.SetMarkerSize(1)
datLegend.SetMarkerStyle(20)

leg= TLegend(0.50, 0.75,0.85, 0.85)
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite)
leg.AddEntry( sigLegend, "Signal + Background Fit",  "L" )
leg.AddEntry( bkgLegend, "Background Only Fit",  "L" )
leg.AddEntry( datLegend, "Data ("+year+")",  "P" )



hframe.Draw("");



leg.Draw("same")
frame .Draw("same");
frame.SetYTitle("Events");

cmsTag=TLatex(0.22,0.93,"#scale[1.2]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()
#cmsTag2=TLatex(0.33,0.87,"#scale[1.0]{#bf{#it{Preliminary}}}")
#cmsTag2.SetNDC()
#cmsTag2.SetTextAlign(11)
#cmsTag2.Draw()
cmsTag3=TLatex(0.92,0.93,"#scale[0.65]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
cmsTag3.Draw()


#frame.SetXTitle("m_{#mu#mu} [GeV]");    


ratioPad = TPad("BottomPad","",0.,0.08,1.,0.285);
ratioPad.SetFillStyle(4000);
ratioPad.SetBottomMargin(0);
ratioPad.Draw();

ratioPad.SetFillStyle(4000)
ratioPad.cd();

 
for i in range(1,testhisto.GetNbinsX()+1): testhisto.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))

Ratio = testhisto.Clone("Ratio");
RatioSB = testhisto.Clone("RatioSB");
for i in range(1,testhisto.GetNbinsX()+1):
    Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i))/TMath.Sqrt(testhisto.GetBinContent(i)));
    Ratio.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
    RatioSB.SetBinContent(i, (testhisto.GetBinContent(i) - sbpdf.GetBinContent(i))/TMath.Sqrt(testhisto.GetBinContent(i)));
    RatioSB.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))

Ratio.SetFillColor(4);
Ratio.SetBarWidth(0.8);
Ratio.SetMarkerColor(1);
Ratio.SetMarkerSize(1);
Ratio.SetMarkerStyle(20);
Ratio.SetLineColor(1);
Ratio.SetStats(0)

RatioSB.SetBarWidth(0.8);
RatioSB.SetMarkerColor(1);
RatioSB.SetMarkerSize(1);
RatioSB.SetMarkerStyle(20);
RatioSB.SetLineColor(kOrange-3);
RatioSB.SetLineWidth(2);
RatioSB.IsTransparent();
RatioSB.SetStats(0)

#Ratio.SetXTitle("m_{#mu#mu} [GeV]");

Ratio.GetXaxis().SetLabelSize(0);
Ratio.GetXaxis().SetLabelFont(43);
Ratio.GetXaxis().SetLabelOffset(0.012);
Ratio.GetXaxis().SetTitleFont(0);
Ratio.GetXaxis().SetTitleSize(22);


Ratio.GetYaxis().SetLabelOffset(0.012);
Ratio.GetYaxis().SetLabelSize(20);
Ratio.GetYaxis().SetLabelFont(43);
Ratio.GetYaxis().SetTitleFont(63);
Ratio.GetYaxis().SetTitleSize(22);
Ratio.GetYaxis().SetNdivisions(505);
Ratio.GetYaxis().SetRangeUser(-3.9, 3.9);
Ratio.GetYaxis().SetTitleOffset(1.3);
Ratio.GetYaxis().SetTitle("Pull")

Ratio.Draw("hist b");
RatioSB.Draw("hist same");



#TLine *line1 = TLine(xmin, 1, xmax, 1);
#line1.SetLineColor(kRed);
line1 = TLine(xmin, 0, xmax, 0);
line1.SetLineColor(kBlack);
line1.SetLineWidth(1);
line1.Draw("same");


c1.SaveAs("fit_m"+mass+"_"+year+".pdf");
