import os,sys
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")

year = sys.argv[1]
idx_CR = sys.argv[2]

if year == "2018":
        lumi = 61.3

if year == "2017":
        lumi = 35.3

if idx_CR == "1":
        file = "1.584_207"
        backname="dkk"
if idx_CR == "2":
        file = "1.716_215"        
        backname="dkpi"
#file = "1.716_215"
mass = file[:5]

os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.1,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")


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

data = w.data("data_obs").reduce("CMS_channel==1")
data = data.reduce(RooFit.SelectVars(m2muvars))  
rdh = data.binnedClone()
testhisto = rdh.createHistogram("testhisto", m2mu, RooFit.Binning(nbins, xmin, xmax));
#testhisto2 = rdh.createHistogram("testhisto2", m2mu, RooFit.Binning(nbins, xmin, xmax));
ymin=testhisto.GetMinimum()
ymax=testhisto.GetMaximum()


gStyle.SetOptTitle(0); 
gStyle.SetOptStat(0);

gStyle.SetPadTopMargin(0.08);
gStyle.SetPadBottomMargin(0.45);
gStyle.SetPadLeftMargin(0.15);
gStyle.SetPadRightMargin(0.07);

gStyle.SetNdivisions(508, "X");
gStyle.SetNdivisions(508, "Y");

c1 = TCanvas("c1", "c1", 650, 720);
c1.SetLogx(0);
c1.SetLogy(0);
c1.SetBottomMargin(0.45);
c1.SetFillStyle(4000)


hframe = TH2F("hframe","hframe",500,xmin,xmax,500,0.8*ymin,1.2*ymax);

hframe.SetYTitle("Events / 0.002 GeV");
hframe.SetXTitle("m_{#mu#mu} [GeV]");

hframe.GetXaxis().SetLabelSize(20);
hframe.GetXaxis().SetLabelFont(43);
hframe.GetXaxis().SetTitleFont(43);
hframe.GetXaxis().SetTitleSize(22);
hframe.GetXaxis().SetTitleOffset(8.0);

hframe.GetXaxis().SetLabelOffset(0.360);
hframe.GetXaxis().SetNdivisions(505);

hframe.GetYaxis().SetLabelSize(20);
hframe.GetYaxis().SetLabelFont(43);
hframe.GetYaxis().SetTitleFont(43);
hframe.GetYaxis().SetTitleSize(22);
hframe.GetYaxis().SetTitleOffset(2.0);
hframe.GetYaxis().SetNdivisions(508);
hframe.GetYaxis().SetMaxDigits(2)

frame = m2mu.frame();
data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(1));

totalpdf = w.pdf("pdf_binhighip_nuis")

bkgpdf = wb.pdf("pdf_binhighip_nuis")
bkgpdf.plotOn(frame, RooFit.LineColor(4),  RooFit.Name("new_pdf"));

bkgpdf1 = wb.pdf("shapeBkg_bkg_mass_highip")
#bkgpdf1.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"));
nbkg=wb.function("n_exp_binhighip_proc_bkg_mass").getVal()

#bkgpdf1.plotOn(frame, RooFit.LineColor(20), RooFit.Name("new_pdf"),RooFit.Normalization(nbkg,RooAbsReal.NumEvent));

if idx_CR == "1":
    bkgpdf2 = wb.pdf("shapeBkg_dkk_mass_highip")
    ndkk = wb.function("n_exp_binhighip_proc_dkk_mass").getVal()

if idx_CR == "2":
    bkgpdf2 = wb.pdf("shapeBkg_dkpi_mass_highip")
    ndkk = wb.function("n_exp_binhighip_proc_dkpi_mass").getVal()

frame2=m2mu.frame()

sigpdf = w.pdf("shapeSig_signalModel_generic_highip")
nsig = w.function("n_exp_binhighip_proc_signalModel_generic").getVal()

if idx_CR == "1":
    bkgpdf2.plotOn(frame2, RooFit.LineColor(8), RooFit.LineStyle(6), RooFit.LineWidth(4), RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent));
if idx_CR == "2":
    bkgpdf2.plotOn(frame2, RooFit.LineColor(222),RooFit.LineStyle(2), RooFit.LineWidth(4), RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent));
#####/my chi2 calculation
sbpdf   = bkgpdf .createHistogram("bkgpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
sbpdfnorm=sbpdf.Integral(1,-1)
sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));

testpdf   = bkgpdf1 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf.Scale(nbkg/testpdf.Integral(1,-1))

testpdf2   = bkgpdf .createHistogram("testpdf2"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf2.Scale(testhisto.Integral(1,-1)/testpdf2.Integral(1,-1));

hframe.Draw("");
frame .Draw("same");

sigLegend = TH2F();
sigLegend.SetLineColor(kOrange-3)
sigLegend.SetLineWidth(3)

bkgLegend = TH2F();
bkgLegend.SetLineColor(4)
bkgLegend.SetLineWidth(3)

datLegend = TH2F()
datLegend.SetMarkerSize(1)
datLegend.SetMarkerStyle(20)

dkkLegend = TH2F()
dkkLegend.SetLineColor(8)
dkkLegend.SetLineWidth(3)
dkkLegend.SetLineStyle(6)

dkpLegend = TH2F()
dkpLegend.SetLineColor(222)
dkpLegend.SetLineWidth(3)
dkpLegend.SetLineStyle(2)

leg= TLegend(0.45, 0.79,0.85, 0.89)
leg.SetBorderSize( 0 )
leg.SetFillStyle( 1001 )
leg.SetFillColor(kWhite)
leg.AddEntry( bkgLegend, "Background Only Fit",  "L" )
leg.AddEntry( datLegend, "Data ("+year+")",  "P" )
leg.SetTextSize(0.035)

leg2= TLegend(0.45, 0.74,0.85, 0.89)
leg2.SetBorderSize( 0 )
leg2.SetFillStyle( 1001 )
leg2.SetTextSize(0.09)
leg2.SetFillColor(kWhite)

if idx_CR == "1":
    leg2.AddEntry( dkkLegend, "D#rightarrow KK", "L")
if idx_CR == "2":
    leg2.AddEntry( dkpLegend, "D#rightarrow K#pi","L")


leg.Draw("same")

cmsTag=TLatex(0.22,0.93,"#scale[1.2]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()

cmsTag3=TLatex(0.92,0.93,"#scale[0.65]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
cmsTag3.Draw()

latex4 = TLatex();
latex4.SetNDC();
latex4.SetTextSize(0.04);
latex4.SetTextAlign(33);
h_string4 = "Control region";
latex4.DrawLatex(0.5,0.57, h_string4);

latex5 = TLatex();
latex5.SetNDC();
latex5.SetTextSize(0.04);
latex5.SetTextAlign(33);
h_string5 = "3.5#sigma_{L} < L < 11#sigma_{L}";
latex5.DrawLatex(0.5,0.52, h_string5);


ratioPad = TPad("BottomPad","",0.,0.08,1.,0.45);
ratioPad.SetFillStyle(4000);
ratioPad.SetBottomMargin(0.05);
#ratioPad.SetLogx(1);
ratioPad.Draw();
ratioPad.cd();

for i in range(1,testhisto.GetNbinsX()+1): testhisto.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
Ratio = testhisto.Clone("Ratio");
for i in range(1,testhisto.GetNbinsX()+1):
    Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i)));
    Ratio.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))


Ratio.GetXaxis().SetLabelSize(0);
Ratio.GetXaxis().SetLabelFont(43)
Ratio.GetXaxis().SetLabelOffset(0.012);
Ratio.GetXaxis().SetTitleFont(0);
Ratio.GetXaxis().SetTitleSize(22);


#Ratio.GetYaxis().SetLabelOffset(0.005);
Ratio.GetYaxis().SetLabelSize(20);
Ratio.GetYaxis().SetLabelFont(43);
Ratio.GetYaxis().SetTitleFont(43);
Ratio.GetYaxis().SetTitleSize(22);
Ratio.GetYaxis().SetNdivisions(505);
Ratio.GetYaxis().SetTitleOffset(2.0);
Ratio.GetYaxis().SetTitle("Data - Comb. Bkg.")
Ratio.SetMaximum(Ratio.GetMaximum()*1.8)
#Ratio.SetMinimum(-1200)


Ratio.SetMarkerColor(1);
Ratio.SetLineColor(1);
Ratio.SetMarkerSize(1)
Ratio.SetMarkerStyle(20)


Ratio.Draw("ep");

frame2.Draw("same")

line1 = TLine(xmin, 0, xmax, 0);
line1.SetLineColor(kBlack);
line1.SetLineWidth(1);
line1.Draw("same");

leg2.Draw("same")

c1.SaveAs("bump_CR_"+backname+"_"+year+".pdf");
