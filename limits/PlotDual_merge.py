import os
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")

year="2018"
#fname = ["1.600_208","1.665_212","1.682_213"]
#fname = ["1.508_202","1.569_206","1.616_209","1.682_213","1.523_203","1.584_207","1.632_210","1.538_204","1.649_211","1.553_205","1.600_208","1.665_212","1.682_213"]
#fname = ["1.508_202","1.569_206","1.523_203","1.584_207","1.538_204","1.553_205"]
#fname = ["1.616_209","1.682_213","1.632_210","1.649_211","1.600_208","1.665_212","1.699_214"]

lumi = 96.6

fname = ["1.649_211"]
for file in fname:

    mass = file[:5]
    
    os.system("combine -M MultiDimFit -d output_dual/dpCard_2017IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.05,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
    os.system("combine -M MultiDimFit -d output_dual/dpCard_2017IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")

    os.system("mv higgsCombineSB.MultiDimFit.mH"+mass+".root higgsCombineSB.MultiDimFit.mH"+mass+"_2017.root")
    os.system("mv higgsCombineBonly.MultiDimFit.mH"+mass+".root higgsCombineBonly.MultiDimFit.mH"+mass+"_2017.root")

    os.system("combine -M MultiDimFit -d output_dual/dpCard_2018IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.05,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
    os.system("combine -M MultiDimFit -d output_dual/dpCard_2018IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")

    os.system("mv higgsCombineSB.MultiDimFit.mH"+mass+".root higgsCombineSB.MultiDimFit.mH"+mass+"_2018.root")
    os.system("mv higgsCombineBonly.MultiDimFit.mH"+mass+".root higgsCombineBonly.MultiDimFit.mH"+mass+"_2018.root")


    if mass[-1] == '0':
        mass = mass[:-1]
    if mass[-1] == '0':
        mass = mass[:-1]
                                                          
#2017    
    f1 = TFile("higgsCombineSB.MultiDimFit.mH"+mass+"_2017.root","READ")
    fb1 = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+"_2017.root","READ")
    
    w1 = f1.Get("w")
    w1.loadSnapshot("MultiDimFit")
    
    wb1 = fb1.Get("w")
    wb1.loadSnapshot("MultiDimFit")
    
    m2mu = w1.var("m2mu")
    m2muvars = RooArgSet(m2mu)
    xmin = m2mu.getMin()
    xmax = m2mu.getMax()
    print "min=",xmin," max=",xmax
    nbins = m2mu.getBinning().numBins()
    
    data = w1.data("data_obs").reduce("CMS_channel==0")
    data = data.reduce(RooFit.SelectVars(m2muvars))  
    rdh1 = data.binnedClone()
    testhisto = rdh1.createHistogram("testhisto", m2mu, RooFit.Binning(nbins, xmin, xmax));


#2018
    f2 = TFile("higgsCombineSB.MultiDimFit.mH"+mass+"_2018.root","READ")
    fb2 = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+"_2018.root","READ")

    w2 = f2.Get("w")
    w2.loadSnapshot("MultiDimFit")

    wb2 = fb2.Get("w")
    wb2.loadSnapshot("MultiDimFit")


    data2 = w2.data("data_obs").reduce("CMS_channel==0")
    data2 = data2.reduce(RooFit.SelectVars(m2muvars))
    rdh2 = data2.binnedClone()
    testhisto2 = rdh2.createHistogram("testhisto2", m2mu, RooFit.Binning(nbins, xmin, xmax));

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
    hframe.GetYaxis().SetNdivisions(505);
    hframe.GetYaxis().SetMaxDigits(2)
    #hframe.GetYaxis().SetLabelOffset(0.007);
    
          
          
    frame = m2mu.frame();

    #mark1
    data.append(data2);

    data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(1));


    frac11=w1.function("n_exp_binCatAB_proc_signalModel_generic")
    frac12=w1.function("n_exp_binCatAB_proc_dkk_mass")
    frac13=w1.function("n_exp_binCatAB_proc_dkpi_mass")
    frac14=w1.function("n_exp_binCatAB_proc_bkg_mass")


    frac21=w2.function("n_exp_binCatAB_proc_signalModel_generic")
    frac22=w2.function("n_exp_binCatAB_proc_dkk_mass")
    frac23=w2.function("n_exp_binCatAB_proc_dkpi_mass")
    frac24=w2.function("n_exp_binCatAB_proc_bkg_mass")

    nsig1 = w1.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk1 = w1.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi1 = w1.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg1 = w1.function("n_exp_binCatAB_proc_bkg_mass").getVal()

    nsig2 = w2.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk2 = w2.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi2 = w2.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg2 = w2.function("n_exp_binCatAB_proc_bkg_mass").getVal()


    ntotal1 = RooRealVar("ntotal1","ntotal1",1000,0,1e9)
    ntotal2 = RooRealVar("ntotal2","ntotal2",1000,0,1e9)

    ntotal1.setVal(nsig1+ndkk1+ndkpi1+nbkg1)
    ntotal2.setVal(nsig2+ndkk2+ndkpi2+nbkg2)
    totalpdf1 = w1.pdf("pdf_binCatAB_nuis")
    totalpdf2 = w2.pdf("pdf_binCatAB_nuis")
    totalpdf = RooAddPdf("totalpdf","totalpdf",RooArgList(totalpdf1,totalpdf2),RooArgList(ntotal1,ntotal2))

    totalpdf.plotOn(frame, RooFit.LineColor(kOrange-3), RooFit.Name("total_pdf"))


    #mark2
    nsig1_b = wb1.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk1_b = wb1.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi1_b = wb1.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg1_b = wb1.function("n_exp_binCatAB_proc_bkg_mass").getVal()

    nsig2_b = wb2.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk2_b = wb2.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi2_b = wb2.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg2_b = wb2.function("n_exp_binCatAB_proc_bkg_mass").getVal()

    ntotal1_b = RooRealVar("ntotal1_b","ntotal1_b",1000,0,1e9)
    ntotal2_b = RooRealVar("ntotal2_b","ntotal2_b",1000,0,1e9)

    ntotal1_b.setVal(nsig1_b+ndkk1_b+ndkpi1_b+nbkg1_b)
    ntotal2_b.setVal(nsig2_b+ndkk2_b+ndkpi2_b+nbkg2_b)

    bkgpdf1 = wb1.pdf("pdf_binCatAB_nuis")
    bkgpdf2 = wb2.pdf("pdf_binCatAB_nuis")
    bkgpdf = RooAddPdf("bkgpdf","bkgpdf",RooArgList(bkgpdf1,bkgpdf2),RooArgList(ntotal1_b,ntotal2_b))
    bkgpdf.plotOn(frame, RooFit.LineColor(4),RooFit.LineWidth(2), RooFit.Name("new_pdf"));
    
    #mark3
    bkgpdf11 = w1.pdf("shapeBkg_bkg_mass_CatAB")
    bkgpdf21 = w2.pdf("shapeBkg_bkg_mass_CatAB")
    bkgpdf01 = RooAddPdf("bkgpdf01","bkgpdf01",RooArgList(bkgpdf11,bkgpdf21),RooArgList(frac11,frac21))
#    bkgpdf01.plotOn(frame, RooFit.LineColor(4), RooFit.Name("new_pdf_2"));
    
    #mark4
    bkgpdf12 = w1.pdf("shapeBkg_dkk_mass_CatAB")
    bkgpdf22 = w2.pdf("shapeBkg_dkk_mass_CatAB")
    bkgpdf02 = RooAddPdf("bkgpdf02","bkgpdf02",RooArgList(bkgpdf12,bkgpdf22),RooArgList(frac12,frac22))

    #mark5
    bkgpdf13 = w1.pdf("shapeBkg_dkpi_mass_CatAB")
    bkgpdf23 = w2.pdf("shapeBkg_dkpi_mass_CatAB")
    bkgpdf03 = RooAddPdf("bkgpdf03","bkgpdf03",RooArgList(bkgpdf13,bkgpdf23),RooArgList(frac13,frac23))    

    frame2=m2mu.frame()
    sigpdf1 = w1.pdf("shapeSig_signalModel_generic_CatAB")
    sigpdf2 = w2.pdf("shapeSig_signalModel_generic_CatAB")
 
    
    nsigv1 = RooRealVar("nsigv1","nsigv1",1000,0,1e9)
    ndkkv1 = RooRealVar("ndkkv1","ndkkv1",1000,0,1e9)
    ndkpiv1 = RooRealVar("ndkpiv1","ndkpiv1",1000,0,1e9) 

    nsigv2 = RooRealVar("nsigv2","nsigv2",1000,0,1e9)
    ndkkv2 = RooRealVar("ndkkv2","ndkkv2",1000,0,1e9)
    ndkpiv2 = RooRealVar("ndkpiv2","ndkpiv2",1000,0,1e9)


    nsigv1.setVal(nsig1)
    ndkkv1.setVal(ndkk1)
    ndkpiv1.setVal(ndkpi1)
    nsigv1.setConstant(True)
    ndkkv1.setConstant(True)
    ndkpiv1.setConstant(True)

    nsigv2.setVal(nsig2)
    ndkkv2.setVal(ndkk2)
    ndkpiv2.setVal(ndkpi2)
    nsigv2.setConstant(True)
    ndkkv2.setConstant(True)
    ndkpiv2.setConstant(True)

    addpeak = RooAddPdf("peak", "sig+dkk+dkpi", RooArgList(sigpdf1,bkgpdf12,bkgpdf13,sigpdf2,bkgpdf22,bkgpdf23),  RooArgList(nsigv1,ndkkv1,ndkpiv1,nsigv2,ndkkv2,ndkpiv2))


    testpeak   = addpeak.createHistogram("testpeak"  , m2mu, RooFit.Binning(nbins, xmin, xmax));

    print "all=",nsig1+ndkk1+ndkpi1+nsig2+ndkk2+ndkpi2
    addpeak.plotOn(frame2, RooFit.LineColor(kOrange-3), RooFit.Name("sig_pdf"), RooFit.Normalization(nsig1+ndkk1+ndkpi1+nsig2+ndkk2+ndkpi2,RooAbsReal.NumEvent))

    bkgpdf02.plotOn(frame2, RooFit.LineColor(8), RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk1+ndkk2,RooAbsReal.NumEvent));
    bkgpdf03.plotOn(frame2, RooFit.LineColor(9), RooFit.Name("dkpi_pdf"),RooFit.Normalization(ndkpi1+ndkpi2,RooAbsReal.NumEvent));
    
    testpdf   = bkgpdf01 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
    testpdf.Scale((nbkg1+nbkg2)/testpdf.Integral(1,-1));
    
    testpdf2   = bkgpdf .createHistogram("testpdf2"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
    testpdf2.Scale(testhisto.Integral(1,-1)/testpdf2.Integral(1,-1));
    
    sbpdf   = totalpdf .createHistogram("totalpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
    sbpdf.Scale(testhisto.Integral(1,-1)/sbpdf.Integral(1,-1));
    
    
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
    #Ratio.SetXTitle("m_{#mu#mu} [GeV]");
    #Ratio.Draw("histsame");
    #frame2.Draw("same")
    Ratio.Draw("ep");
    
    frame2.Draw("same")
    
    #TLine *line1 = TLine(xmin, 1, xmax, 1);
    #line1.SetLineColor(kRed);
    line1 = TLine(xmin, 0, xmax, 0);
    line1.SetLineColor(kBlack);
    line1.SetLineWidth(1);
    line1.Draw("same");
    
    c1.SaveAs("test_"+mass+".png");
    #c1.SaveAs("test_m1.584.pdf");
    c1.SaveAs("test_m"+mass+"_combinedBKG.pdf");
