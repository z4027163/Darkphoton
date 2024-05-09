import os
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")
                                                      

file = "1.584_207"
#file = "1.733_216"
mass = file[:5]


lumi = 96.6

year="2017"

os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.1,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")

os.system("mv higgsCombineSB.MultiDimFit.mH"+mass+".root higgsCombineSB.MultiDimFit.mH"+mass+"_2017.root")
os.system("mv higgsCombineBonly.MultiDimFit.mH"+mass+".root higgsCombineBonly.MultiDimFit.mH"+mass+"_2017.root")

year="2018"

os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.1,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")

os.system("mv higgsCombineSB.MultiDimFit.mH"+mass+".root higgsCombineSB.MultiDimFit.mH"+mass+"_2018.root")
os.system("mv higgsCombineBonly.MultiDimFit.mH"+mass+".root higgsCombineBonly.MultiDimFit.mH"+mass+"_2018.root")


f1 = TFile("higgsCombineSB.MultiDimFit.mH"+mass+"_2017.root","READ")
fb1 = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+"_2017.root","READ")
f2 = TFile("higgsCombineSB.MultiDimFit.mH"+mass+"_2018.root","READ")
fb2 = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+"_2018.root","READ")


w1 = f1.Get("w")
w1.loadSnapshot("MultiDimFit")

wb1 = fb1.Get("w")
wb1.loadSnapshot("MultiDimFit")

w2 = f2.Get("w")
w2.loadSnapshot("MultiDimFit")

wb2 = fb2.Get("w")
wb2.loadSnapshot("MultiDimFit")



m2mu = w2.var("m2mu")
m2muvars = RooArgSet(m2mu)
xmin = m2mu.getMin()
xmax = m2mu.getMax()
nbins = m2mu.getBinning().numBins()

data = w1.data("data_obs").reduce("CMS_channel==1")
data = data.reduce(RooFit.SelectVars(m2muvars))  
rdh1 = data.binnedClone()
testhisto = rdh1.createHistogram("testhisto", m2mu, RooFit.Binning(nbins, xmin, xmax));

nentry1=data.sumEntries()

data2 = w2.data("data_obs").reduce("CMS_channel==1")
data2 = data2.reduce(RooFit.SelectVars(m2muvars))
rdh2 = data2.binnedClone()
testhisto2 = rdh2.createHistogram("testhisto2", m2mu, RooFit.Binning(nbins, xmin, xmax));

nentry2=data2.sumEntries()

testhisto.Add(testhisto2)

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
hframe.SetYTitle("Events / 0.002 GeV");
hframe.SetXTitle("m_{#mu#mu} [GeV]");    

hframe.GetXaxis().SetLabelSize(18);
hframe.GetXaxis().SetLabelFont(43);
hframe.GetXaxis().SetTitleFont(63);
hframe.GetXaxis().SetTitleSize(22);
hframe.GetXaxis().SetTitleOffset(5.2);

hframe.GetXaxis().SetLabelOffset(0.215);
hframe.GetXaxis().SetNdivisions(505);

hframe.GetYaxis().SetLabelSize(15);
hframe.GetYaxis().SetLabelFont(43);
hframe.GetYaxis().SetTitleFont(63);
hframe.GetYaxis().SetTitleSize(22);
hframe.GetYaxis().SetTitleOffset(1.8);
hframe.GetYaxis().SetNdivisions(508);
hframe.GetYaxis().SetMaxDigits(2)

      
frame = m2mu.frame();
data.append(data2);

data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(1));

frac11=w1.function("n_exp_binhighip_proc_signalModel_generic")
frac12=w1.function("n_exp_binhighip_proc_dkk_mass")
#frac12=w1.function("n_exp_binhighip_proc_dkpi_mass")
frac13=w1.function("n_exp_binhighip_proc_bkg_mass")

frac21=w2.function("n_exp_binhighip_proc_signalModel_generic")
frac22=w2.function("n_exp_binhighip_proc_dkk_mass")
#frac22=w2.function("n_exp_binhighip_proc_dkpi_mass")
frac23=w2.function("n_exp_binhighip_proc_bkg_mass")

nsig1=w1.function("n_exp_binhighip_proc_signalModel_generic").getVal()
ndkk1=w1.function("n_exp_binhighip_proc_dkk_mass").getVal()
#ndkk1=w1.function("n_exp_binhighip_proc_dkpi_mass").getVal()
nbkg1=w1.function("n_exp_binhighip_proc_bkg_mass").getVal()

nsig2=w2.function("n_exp_binhighip_proc_signalModel_generic").getVal()
ndkk2=w2.function("n_exp_binhighip_proc_dkk_mass").getVal()
#ndkk2=w2.function("n_exp_binhighip_proc_dkpi_mass").getVal()
nbkg2=w2.function("n_exp_binhighip_proc_bkg_mass").getVal()


nsigv1 = RooRealVar("nsigv1","nsigv1",1000,0,1e9)
ndkkv1 = RooRealVar("ndkkv1","ndkkv1",1000,0,1e9)
nbkgv1 = RooRealVar("nbkgv1","nbkgv1",1000,0,1e9)
nsigv1.setVal(nsig1)
ndkkv1.setVal(ndkk1)
nbkgv1.setVal(nbkg1)
nsigv1.setConstant(True)
ndkkv1.setConstant(True)
nbkgv1.setConstant(True)

nsigv2 = RooRealVar("nsigv2","nsigv2",1000,0,1e9)
ndkkv2 = RooRealVar("ndkkv2","ndkkv2",1000,0,1e9)
nbkgv2 = RooRealVar("nbkgv2","nbkgv2",1000,0,1e9)
nsigv2.setVal(nsig2)
ndkkv2.setVal(ndkk2)
nbkgv2.setVal(nbkg2)
nsigv2.setConstant(True)
ndkkv2.setConstant(True)
nbkgv2.setConstant(True)

sigpdf1=w1.pdf("shapeSig_signalModel_generic_highip")
bkgpdf11=w1.pdf("shapeBkg_dkk_mass_highip")
#bkgpdf11=w1.pdf("shapeBkg_dkpi_mass_highip")
bkgpdf12=w1.pdf("shapeBkg_bkg_mass_highip")

sigpdf2=w2.pdf("shapeSig_signalModel_generic_highip")
bkgpdf21=w2.pdf("shapeBkg_dkk_mass_highip")
#bkgpdf21=w2.pdf("shapeBkg_dkpi_mass_highip")
bkgpdf22=w2.pdf("shapeBkg_bkg_mass_highip")


print "total=",nsig1+ndkk1+nbkg1+nsig2+ndkk2+nbkg2
#totalpdf = RooAddPdf("total_pdf", "total_pdf", RooArgList(sigpdf1,bkgpdf11,bkgpdf12,sigpdf2,bkgpdf21,bkgpdf22),  RooArgList(nsigv1,ndkkv1,nbkgv1,nsigv2,ndkkv2,nbkgv2))
#totalpdf.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"),RooFit.Normalization(nsig1+ndkk1+nbkg1+nsig2+ndkk2+nbkg2,RooAbsReal.NumEvent));
#totalpdf = RooAddPdf("total_pdf", "total_pdf", RooArgList(bkgpdf12,bkgpdf22),  RooArgList(nbkgv1,nbkgv2))
#totalpdf.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"),RooFit.Normalization(nbkg1+nbkg2,RooAbsReal.NumEvent));


ntotal1 = RooRealVar("ntotal1","ntotal1",1000,0,1e9)
ntotal2 = RooRealVar("ntotal2","ntotal2",1000,0,1e9)



ntotal1.setVal(nentry1)
ntotal2.setVal(nentry2*2)
ntotal1.setConstant(True)
ntotal2.setConstant(True)

print "ndkk1=",ndkk1
print "nbkg2=",nbkg2
print "total1=",nsig1+ndkk1+nbkg1
print "entry1=",nentry1
print "ndkk2=",ndkk2
print "nbkg2=",nbkg2
print "total2=",nsig2+ndkk2+nbkg2
print "entry2=",nentry2


frac00 = RooRealVar("frac00","frac00",0.57316,0.5,0.6)

pdf111 = w1.pdf("pdf_binhighip_nuis")

pdf222 = w2.pdf("pdf_binhighip_nuis")
print "combining"
totalpdf = RooAddPdf("totalpdf","totalpdf",pdf222,pdf111,frac00)
print "finish combine"

totalpdf.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"),RooFit.Normalization(nentry1+nentry2,RooAbsReal.NumEvent));

frac11_b=wb1.function("n_exp_binhighip_proc_signalModel_generic")
frac12_b=wb1.function("n_exp_binhighip_proc_dkk_mass")
#frac12_b=wb1.function("n_exp_binhighip_proc_dkpi_mass")
frac13_b=wb1.function("n_exp_binhighip_proc_bkg_mass")

frac21_b=wb2.function("n_exp_binhighip_proc_signalModel_generic")
frac22_b=wb2.function("n_exp_binhighip_proc_dkk_mass")
#frac22_b=wb2.function("n_exp_binhighip_proc_dkpi_mass")
frac23_b=wb2.function("n_exp_binhighip_proc_bkg_mass")

nsig1_b=wb1.function("n_exp_binhighip_proc_signalModel_generic").getVal()
ndkk1_b=wb1.function("n_exp_binhighip_proc_dkk_mass").getVal()
#ndkk1_b=wb1.function("n_exp_binhighip_proc_dkpi_mass").getVal()

nbkg1_b=1
nbkg1_b=wb1.function("n_exp_binhighip_proc_bkg_mass").getVal()



nsig2_b=wb2.function("n_exp_binhighip_proc_signalModel_generic").getVal()
ndkk2_b=wb2.function("n_exp_binhighip_proc_dkk_mass").getVal()
#ndkk2_b=wb2.function("n_exp_binhighip_proc_dkpi_mass").getVal()
nbkg2_b=wb2.function("n_exp_binhighip_proc_bkg_mass").getVal()


ntotal1_b = RooRealVar("ntotal1_b","ntotal1_b",1000,0,1e9)
ntotal2_b = RooRealVar("ntotal2_b","ntotal2_b",1000,0,1e9)

ntotal1_b.setVal(ndkk1_b+nbkg1_b)
ntotal2_b.setVal(ndkk2_b+nbkg2_b)
ntotal1_b.setConstant(True)
ntotal2_b.setConstant(True)


bkgpdf1 = wb1.pdf("pdf_binhighip_nuis")
bkgpdf2 = wb2.pdf("pdf_binhighip_nuis")
bkgpdf = RooAddPdf("bkgpdf","bkgpdf",RooArgList(bkgpdf1,bkgpdf2),RooArgList(ntotal1_b,ntotal2_b))

#bkgpdf.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"),RooFit.Normalization(ndkk1_b+nbkg1_b+ndkk2_b+nbkg2_b,RooAbsReal.NumEvent));

bkgpdf11 = wb1.pdf("shapeBkg_bkg_mass_highip")
bkgpdf21 = wb2.pdf("shapeBkg_bkg_mass_highip")
bkgpdf01 = RooAddPdf("bkgpdf01","bkgpdf01",RooArgList(bkgpdf11,bkgpdf21),RooArgList(frac13_b,frac23_b))

#bkgpdf01.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"));

bkgpdf12 = wb1.pdf("shapeBkg_dkk_mass_highip")
bkgpdf22 = wb2.pdf("shapeBkg_dkk_mass_highip")
bkgpdf02 = RooAddPdf("bkgpdf02","bkgpdf02",RooArgList(bkgpdf12,bkgpdf22),RooArgList(frac12_b,frac22_b))
#bkgpdf2.plotOn(frame, RooFit.LineColor(8), RooFit.Name("new_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent));


frame2=m2mu.frame()

#sigpdf = w.pdf("shapeSig_signalModel_generic_highip")
#nsig = w.function("n_exp_binhighip_proc_signalModel_generic").getVal()
#sigpdf.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig,RooAbsReal.NumEvent))

#addpeak = RooAddPdf("peak", "sig+dkk", RooArgList(sigpdf,bkgpdf),  RooArgList(nsig,ndkk))
#addpeak.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig+ndkk,RooAbsReal.NumEvent))

bkgpdf02.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk1_b+ndkk2_b,RooAbsReal.NumEvent));
#####/my chi2 calculation
sbpdf   = totalpdf.createHistogram("totalpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
sbpdfnorm=sbpdf.Integral(1,-1)
sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));

testpdf   = bkgpdf .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
#testpdf.Scale(nbkg/testpdf.Integral(1,-1))
testpdf.Scale(testhisto.Integral(1,-1)/testpdf.Integral(1,-1))

hframe.Draw("");
frame .Draw("same");
testhisto.SetLineColor(kRed)
testhisto.Draw("HIST same")

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
leg.AddEntry( datLegend, "Data",  "P" )

leg.Draw("same")

cmsTag=TLatex(0.22,0.93,"#scale[1.2]{CMS}")
cmsTag.SetNDC()
cmsTag.SetTextAlign(11)
cmsTag.Draw()

cmsTag3=TLatex(0.92,0.93,"#scale[0.65]{#bf{"+str(lumi)+" fb^{-1} (13 TeV)}}")
cmsTag3.SetNDC()
cmsTag3.SetTextAlign(31)
cmsTag3.Draw()


ratioPad = TPad("BottomPad","",0.,0.08,1.,0.285);
ratioPad.SetFillStyle(4000);
ratioPad.SetBottomMargin(0);
#ratioPad.SetLogx(1);
ratioPad.Draw();
ratioPad.cd();

#Ratio.Divide(testhisto, testpdf);
for i in range(1,testhisto.GetNbinsX()+1): testhisto.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
Ratio = testhisto.Clone("Ratio");
for i in range(1,testhisto.GetNbinsX()+1):
    Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i)));
    Ratio.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))

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
Ratio.GetYaxis().SetTitleOffset(1.8);
Ratio.GetYaxis().SetTitle("Data-Fit")



Ratio.SetMarkerColor(1);
Ratio.SetLineColor(1);
Ratio.SetMarkerSize(1)
Ratio.SetMarkerStyle(20)
Ratio.Draw("ep");
#Ratio2.Draw("ep same")

#hframe2.Draw("");
#frame2.Draw("same")

#line1.SetLineColor(kRed);
line1 = TLine(xmin, 0, xmax, 0);
line1.SetLineColor(kBlack);
line1.SetLineWidth(1);
line1.Draw("same");


c1.SaveAs("test2.png");
#c1.SaveAs("test_m1.584.pdf");
c1.SaveAs("test2_m"+mass+"_combinedBKG.pdf");
