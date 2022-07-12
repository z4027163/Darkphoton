import os
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")
                                                      
year="2017"
file = "1.584_207"
#file = "1.733_216"
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

data = w.data("data_obs").reduce("CMS_channel==0")
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
totalpdf = w.pdf("pdf_binCatAB_nuis")
totalpdf.plotOn(frame, RooFit.LineColor(kOrange+1), RooFit.Name("total_pdf"))

bkgpdf = wb.pdf("pdf_binCatAB_nuis")
bkgpdf.plotOn(frame, RooFit.LineColor(8), RooFit.Name("new_pdf"));

bkgpdf1 = w.pdf("shapeBkg_bkg_mass_CatAB")
bkgpdf1.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf_2"));

bkgpdf2 = w.pdf("shapeBkg_dkk_mass_CatAB")
ndkk = w.function("n_exp_binCatAB_proc_dkk_mass").getVal()
#bkgpdf2.plotOn(frame, RooFit.LineColor(8), RooFit.Name("new_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent));


print "nkk_catAB=", ndkk
print "nkk_highip=",w.function("n_exp_binhighip_proc_dkk_mass").getVal()

frame2=m2mu.frame()
sigpdf = w.pdf("shapeSig_signalModel_generic_CatAB")
nsig = w.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
#sigpdf.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig,RooAbsReal.NumEvent))

frac1=w.function("n_exp_binCatAB_proc_signalModel_generic")
frac2=w.function("n_exp_binCatAB_proc_dkk_mass")
#bkgpdf2.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("new_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent))
addpeak = RooAddPdf("peak", "sig+dkk", RooArgList(sigpdf,bkgpdf2),  RooArgList(frac1,frac2))
addpeak.plotOn(frame2, RooFit.LineColor(kOrange+1), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig+ndkk,RooAbsReal.NumEvent))
bkgpdf2.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk,RooAbsReal.NumEvent));


testpdf   = bkgpdf1 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
#testpdf.Scale(testhisto.Integral(1,-1)/testpdf.Integral(1,-1));
testpdf.Scale(w.function("n_exp_binCatAB_proc_bkg_mass").getVal()/testpdf.Integral(1,-1));

testpdf2   = bkgpdf .createHistogram("testpdf2"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
testpdf2.Scale(testhisto.Integral(1,-1)/testpdf2.Integral(1,-1));

sbpdf   = totalpdf .createHistogram("totalpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));

# std::cout<<"testhisto bins = "<<testhisto.GetNbinsX()<<" and testpdf bins = "<<testpdf.GetNbinsX()<<"\n"; */
# std::cout<<"testhisto integral = "<<testhisto.Integral(1,-1)<<" and testpdf integral = "<<testpdf.Integral(1,-1)<<"\n"; */

mychi2 = 0.;
mychi2sb = 0.;
for i in range(1,testhisto.GetNbinsX()+1):
    tmp_mychi2 = ((testpdf2.GetBinContent(i) - testhisto.GetBinContent(i)) * (testpdf2.GetBinContent(i) - testhisto.GetBinContent(i))) / testpdf2.GetBinContent(i);
    mychi2 += tmp_mychi2;

    tmp_mychi2 = ((sbpdf.GetBinContent(i) - testhisto.GetBinContent(i)) * (sbpdf.GetBinContent(i) - testhisto.GetBinContent(i))) / sbpdf.GetBinContent(i);
    mychi2sb += tmp_mychi2;
    print "bin",i," has chi2=",mychi2sb, " data=", testhisto.GetBinContent(i), " hist=", testpdf.GetBinContent(i)
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
ndkk_2 = w.function("n_exp_binhighip_proc_dkk_mass").getVal()
print "nkk_highip=", ndkk_2
print "nkk_catAB=", ndkk
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

latex4 = TLatex();
latex4.SetNDC();
latex4.SetTextSize(0.03);
latex4.SetTextAlign(33);
h_string4 = year+" IPsig<3.5,PVd<0.2";
latex4.DrawLatex(0.5,0.87, h_string4);

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
    #Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i))/testhisto.GetBinError(i));
    #print testhisto.GetBinContent(i),testpdf.GetBinContent(i),testhisto.GetBinError(i),Ratio.GetBinContent(i)

    Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - testpdf.GetBinContent(i)));
    #Ratio.SetBinContent(i, (testhisto.GetBinContent(i) - sbpdf.GetBinContent(i)));

    Ratio.SetBinError(i, TMath.Sqrt(testhisto.GetBinContent(i)))
#Ratio.SetFillColor(2);
#Ratio.SetLineColor(2);
Ratio.SetMarkerColor(1);
Ratio.SetLineColor(1);
#Ratio.Draw("histsame");
Ratio.Draw("ep");

#hframe2.Draw("");
frame2.Draw("same")

#TLine *line1 = TLine(xmin, 1, xmax, 1);
#line1.SetLineColor(kRed);
line1 = TLine(xmin, 0, xmax, 0);
line1.SetLineColor(kBlack);
line1.SetLineWidth(1);
line1.Draw("same");


c1.SaveAs("test3_m"+mass+".png");
#c1.SaveAs("test_m1.584.pdf");
c1.SaveAs("test3_m"+mass+"_combinedBKG.pdf");
