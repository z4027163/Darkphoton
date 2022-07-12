import os
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")
                                                      
file = "1.858_223"
#file = "1.733_216"

mass = file[:5]

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

gStyle.SetPadTopMargin(0.07);
gStyle.SetPadBottomMargin(0.1);
gStyle.SetPadLeftMargin(0.13);
gStyle.SetPadRightMargin(0.07);

gStyle.SetNdivisions(508, "X");
gStyle.SetNdivisions(508, "Y");

c1 = TCanvas("c1", "c1", 650, 720);
c1.SetLogx(0);
c1.SetLogy(0);
c1.SetBottomMargin(0.32);


hframe = TH2F("hframe","hframe",500,xmin,xmax,500,0.8*ymin,1.2*ymax);
hframe.SetYTitle("Events");
hframe.SetXTitle("m_{#mu#mu} [GeV/c^{2}]");    

hframe.GetXaxis().SetLabelSize(18);
hframe.GetXaxis().SetLabelFont(43);
hframe.GetXaxis().SetTitleFont(63);
hframe.GetXaxis().SetTitleSize(22);
# hframe.GetXaxis().SetTitleOffset(1.5);
# hframe.GetXaxis().SetLabelOffset(3.5);
hframe.GetYaxis().SetLabelSize(20);
hframe.GetYaxis().SetLabelFont(43);
hframe.GetYaxis().SetTitleFont(63);
hframe.GetYaxis().SetTitleSize(22);
hframe.GetYaxis().SetTitleOffset(1.5);
hframe.GetYaxis().SetNdivisions(505);
hframe.GetYaxis().SetTitleOffset(2.1);
#hframe.GetYaxis().SetLabelOffset(0.007);

hframe2= TH2F("hframe2","hframe2",500, xmin, xmax, 500, -4, 4);
hframe2.SetYTitle("(Data-Fit)/Unc.");

hframe2.GetXaxis().SetLabelOffset(1);
hframe2.GetXaxis().SetLabelSize(0.1);  
hframe2.GetYaxis().SetLabelOffset(0.012);
hframe2.GetYaxis().SetLabelSize(0.1);  
hframe2.GetYaxis().SetTitleOffset(0.4);
hframe2.GetYaxis().SetTitleSize(0.1);  
#hframe2.GetYaxis().SetNdivisions(503); 
      
      
frame = m2mu.frame();
data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(0));
totalpdf = w.pdf("pdf_binhighip_nuis")
totalpdf.plotOn(frame, RooFit.LineColor(kOrange+1), RooFit.Name("total_pdf"))

bkgpdf = wb.pdf("pdf_binhighip_nuis")
bkgpdf.plotOn(frame, RooFit.LineColor(8), RooFit.Name("new_pdf"));

bkgpdf1 = wb.pdf("shapeBkg_bkg_mass_highip")
bkgpdf1.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf"));
nbkg=wb.function("n_exp_binhighip_proc_bkg_mass").getVal()

bkgpdf2 = wb.pdf("shapeBkg_dmm_mass_highip")
ndmm = wb.function("n_exp_binhighip_proc_dmm_mass").getVal()


frame2=m2mu.frame()

sigpdf = w.pdf("shapeSig_signalModel_generic_highip")
nsig = w.function("n_exp_binhighip_proc_signalModel_generic").getVal()

frac1=w.function("n_exp_binhighip_proc_signalModel_generic")
frac2=w.function("n_exp_binhighip_proc_dmm_mass")
#bkgpdf2.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("new_pdf"),RooFit.Normalization(ndmm,RooAbsReal.NumEvent))
addpeak = RooAddPdf("peak", "sig+dmm", RooArgList(sigpdf,bkgpdf2),  RooArgList(frac1,frac2))
addpeak.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig+ndmm,RooAbsReal.NumEvent))
bkgpdf2.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("dmm_pdf"),RooFit.Normalization(ndmm,RooAbsReal.NumEvent));

#####/my chi2 calculation
sbpdf   = totalpdf .createHistogram("totalpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
sbpdfnorm=sbpdf.Integral(1,-1)
sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));

testpdf   = bkgpdf1 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf.Scale(nbkg/testpdf.Integral(1,-1))

testpdf2   = bkgpdf .createHistogram("testpdf2"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf2.Scale(testhisto.Integral(1,-1)/testpdf2.Integral(1,-1));

#testpdf.Scale(testhisto.Integral(1,-1)/testpdf.Integral(1,-1))
print "sbpdfnorm=",sbpdf.Integral(1,-1)," pdfnorm=",testpdf.Integral(1,-1)," nexp=",nbkg 
print "dmmnorm=",ndmm
print "data=", testhisto.Integral(1,-1)
#testpdf.Scale(testhisto.Integral(1,-1)/sbpdfnorm);

# std::cout<<"testhisto bins = "<<testhisto.GetNbinsX()<<" and testpdf bins = "<<testpdf.GetNbinsX()<<"\n"; */
# std::cout<<"testhisto integral = "<<testhisto.Integral(1,-1)<<" and testpdf integral = "<<testpdf.Integral(1,-1)<<"\n"; */

mychi2 = 0.;
mychi2sb = 0.;
for i in range(1,testhisto.GetNbinsX()+1):
    tmp_mychi2 = ((testpdf2.GetBinContent(i) - testhisto.GetBinContent(i)) * (testpdf2.GetBinContent(i) - testhisto.GetBinContent(i))) / testpdf2.GetBinContent(i);
    mychi2 += tmp_mychi2;

    tmp_mychi2 = ((sbpdf.GetBinContent(i) - testhisto.GetBinContent(i)) * (sbpdf.GetBinContent(i) - testhisto.GetBinContent(i))) / sbpdf.GetBinContent(i);
    mychi2sb += tmp_mychi2;
#    print "bin",i," has chi2=",mychi2sb, " data=", testhisto.GetBinContent(i), " hist=", testpdf.GetBinContent(i)
    # std::cout<<"bin"<<i<<" has "<<testhisto.GetBinContent(i)<<" and "<<testpdf.GetBinContent(i)<<" entries\n"; *

noofparm=4
mychi2_final = mychi2/(nbins - (noofparm + 1));
mychi2_final_sb = mychi2sb/(nbins - (noofparm + 1));
RooPlot_chi2 = frame.chiSquare("new_pdf", "new_hist", noofparm + 1);
RooPlot_chi2_sb = frame.chiSquare("total_pdf", "new_hist", noofparm + 2);
print "doFit:  mychi^2 = "+str(mychi2)
print "doFit:  mychi^2 s+b = "+str(mychi2sb)
print "doFit:  mychi^2/(nbins-noofparm-1) = "+str(mychi2_final)
print "doFit:  RooPlot chi^2/(nbins-noofparm-1) = "+str(RooPlot_chi2)
print "doFit:  RooPlot s+b chi^2/(nbins-noofparm-1) = "+str(RooPlot_chi2_sb)
print "dmmnorm=",ndmm
#binmax = testhisto.GetMaximumBin(); 
#max_entry = testhisto.GetBinContent(binmax);
#max_entry += max_entry/10.;
#hframe.GetYaxis().SetLimits(0., max_entry);

hframe.Draw("");
frame .Draw("same");

latex2 = TLatex();
latex2.SetNDC();
latex2.SetTextSize(0.03);
latex2.SetTextAlign(33);
#TString h_string = "Chi2/ndf = " + std::to_string(RooPlot_chi2);
#h_string = "mass = " + std::to_string(mh) + "     #chi^{2} = " + std::to_string(RooPlot_chi2) + "     order = " + std::to_string(order);
h_string = "#chi^{2}_{bkg only} = " + str(mychi2)[:6] +  ", dof = " + str(nbins - (noofparm + 1));
latex2.DrawLatex(0.92,0.905, h_string);

latex3 = TLatex();
latex3.SetNDC();
latex3.SetTextSize(0.03);
latex3.SetTextAlign(33);
#TString h_string = "Chi2/ndf = " + std::to_string(RooPlot_chi2);
#h_string = "mass = " + std::to_string(mh) + "     #chi^{2} = " + std::to_string(RooPlot_chi2) + "     order = " + std::to_string(order);
h_string3 = "#chi^{2}_{sig+bkg} = " + str(mychi2sb)[:5] + ", dof = " + str(nbins - (noofparm + 2));
latex3.DrawLatex(0.92,0.855, h_string3);


#paramList = pdf.getParameters(m2mu);
#paramList.printLatex();

ratioPad = TPad("BottomPad","",0.,0.03,1.,0.23);
ratioPad.SetBottomMargin(2.1);
#ratioPad.SetLogx(1);
ratioPad.Draw();
ratioPad.cd();

#Ratio.Divide(testhisto, testpdf);
for i in range(1,testhisto.GetNbinsX()+1): testhisto.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
#Ratio.SetMarkerStyle(20);
#Ratio.SetMarkerColor(kBlack);
#Ratio.Draw("samePE");

#testhisto.Rebin(6);
#testpdf.Rebin(6);
Ratio = testhisto.Clone("Ratio");
for i in range(1,testhisto.GetNbinsX()+1):

    Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i)));
    #Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - sbpdf.GetBinContent(i)));

    Ratio.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
#Ratio.SetFillColor(2);
#Ratio.SetLineColor(2);
Ratio.SetMarkerColor(1);
#Ratio.SetMarkerStyle(22);
Ratio.SetLineColor(1);
#Ratio.SetLineSize(2);
#Ratio.Draw("histsame");

Ratio2 = testhisto.Clone("Ratio2");
for i in range(1,testhisto.GetNbinsX()+1):

    Ratio2.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i)));
    Ratio2.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))

Ratio2.SetMarkerColor(1);
Ratio2.SetLineColor(1);

Ratio.Draw("ep");
#Ratio2.Draw("ep same")

#hframe2.Draw("");
frame2.Draw("same")

#TLine *line1 = TLine(xmin, 1, xmax, 1);
#line1.SetLineColor(kRed);
line1 = TLine(xmin, 0, xmax, 0);
line1.SetLineColor(kBlack);
line1.SetLineWidth(1);
line1.Draw("same");


c1.SaveAs("test4_m"+mass+".png");
#c1.SaveAs("test_m1.584.pdf");
c1.SaveAs("test4_m"+mass+"_combinedBKG.pdf");
