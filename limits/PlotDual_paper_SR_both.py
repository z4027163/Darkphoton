import os,sys
from ROOT import *
gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
gSystem.Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so")
gSystem.AddIncludePath("-I$ROOFITSYS/include")
gSystem.AddIncludePath("-Iinclude/")

year = sys.argv[1]

if year == "2018":
        lumi = 61.3

if year == "2017":
        lumi = 35.3

#fname = ["1.600_208","1.665_212","1.682_213"]
#fname = ["1.508_202","1.569_206","1.616_209","1.682_213","1.523_203","1.584_207","1.632_210","1.538_204","1.649_211","1.553_205","1.600_208","1.665_212","1.682_213"]
#fname = ["1.508_202","1.569_206","1.523_203","1.584_207","1.538_204","1.553_205"]
#fname = ["1.616_209","1.682_213","1.632_210","1.649_211","1.600_208","1.665_212","1.699_214"]


fname = ["1.649_211"]
for file in fname:

    mass = file[:5]
    
    os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.05,5 --cminDefaultMinimizerStrategy 0  -n SB -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")
    os.system("combine -M MultiDimFit -d output_dual/dpCard_"+year+"IterV3_m"+file+".txt --algo none --setParameters r=0.0 --setParameterRanges r=0.001,0.002 --cminDefaultMinimizerStrategy 0  -n Bonly -m "+mass+" --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --cminDefaultMinimizerTolerance=0.001 --saveWorkspace")


    if mass[-1] == '0':
        mass = mass[:-1]
    if mass[-1] == '0':
        mass = mass[:-1]
                                                          
    f1 = TFile("higgsCombineSB.MultiDimFit.mH"+mass+".root","READ")
    fb1 = TFile("higgsCombineBonly.MultiDimFit.mH"+mass+".root","READ")
    
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
    
    hframe.GetXaxis().SetLabelSize(20);
    hframe.GetXaxis().SetLabelFont(43);
    hframe.GetXaxis().SetTitleFont(43);
    hframe.GetXaxis().SetTitleSize(22);
    hframe.GetXaxis().SetTitleOffset(5.2);

    hframe.GetXaxis().SetLabelOffset(0.215);
    hframe.GetXaxis().SetNdivisions(505);

    hframe.GetYaxis().SetLabelSize(20);
    hframe.GetYaxis().SetLabelFont(43);
    hframe.GetYaxis().SetTitleFont(43);
    hframe.GetYaxis().SetTitleSize(22);
    hframe.GetYaxis().SetTitleOffset(2.0);
    hframe.GetYaxis().SetNdivisions(505);
    hframe.GetYaxis().SetMaxDigits(2)
    #hframe.GetYaxis().SetLabelOffset(0.007);
    
          
          
    frame = m2mu.frame();

    #mark1
    data.plotOn(frame, RooFit.Binning(nbins, xmin, xmax), RooFit.Name("new_hist"), RooFit.MarkerSize(1));


    frac11=w1.function("n_exp_binCatAB_proc_signalModel_generic")
    frac12=w1.function("n_exp_binCatAB_proc_dkk_mass")
    frac13=w1.function("n_exp_binCatAB_proc_dkpi_mass")
    frac14=w1.function("n_exp_binCatAB_proc_bkg_mass")


    nsig1 = w1.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk1 = w1.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi1 = w1.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg1 = w1.function("n_exp_binCatAB_proc_bkg_mass").getVal()

    ntotal1 = RooRealVar("ntotal1","ntotal1",1000,0,1e9)
    ntotal2 = RooRealVar("ntotal2","ntotal2",1000,0,1e9)

    ntotal1.setVal(nsig1+ndkk1+ndkpi1+nbkg1)
    totalpdf = w1.pdf("pdf_binCatAB_nuis")

    totalpdf.plotOn(frame, RooFit.LineColor(kOrange-3), RooFit.Name("total_pdf"))


    nsig1_b = wb1.function("n_exp_binCatAB_proc_signalModel_generic").getVal()
    ndkk1_b = wb1.function("n_exp_binCatAB_proc_dkk_mass").getVal()
    ndkpi1_b = wb1.function("n_exp_binCatAB_proc_dkpi_mass").getVal()
    nbkg1_b = wb1.function("n_exp_binCatAB_proc_bkg_mass").getVal()

    ntotal1_b = RooRealVar("ntotal1_b","ntotal1_b",1000,0,1e9)

    ntotal1_b.setVal(nsig1_b+ndkk1_b+ndkpi1_b+nbkg1_b)

    bkgpdf = wb1.pdf("pdf_binCatAB_nuis")
    bkgpdf.plotOn(frame, RooFit.LineColor(4),RooFit.LineWidth(2), RooFit.Name("new_pdf"));
    
    bkgpdf01 = w1.pdf("shapeBkg_bkg_mass_CatAB")
    bkgpdf02 = w1.pdf("shapeBkg_dkk_mass_CatAB")

    bkgpdf03 = w1.pdf("shapeBkg_dkpi_mass_CatAB")
    frame2=m2mu.frame()
    sigpdf1 = w1.pdf("shapeSig_signalModel_generic_CatAB")
    
    nsigv1 = RooRealVar("nsigv1","nsigv1",1000,0,1e9)
    ndkkv1 = RooRealVar("ndkkv1","ndkkv1",1000,0,1e9)
    ndkpiv1 = RooRealVar("ndkpiv1","ndkpiv1",1000,0,1e9) 

    nsigv1.setVal(nsig1)
    ndkkv1.setVal(ndkk1)
    ndkpiv1.setVal(ndkpi1)
    nsigv1.setConstant(True)
    ndkkv1.setConstant(True)
    ndkpiv1.setConstant(True)

    addpeak = RooAddPdf("peak", "sig+dkk+dkpi", RooArgList(sigpdf1,bkgpdf02,bkgpdf03),  RooArgList(nsigv1,ndkkv1,ndkpiv1))

    testpeak   = addpeak.createHistogram("testpeak"  , m2mu, RooFit.Binning(nbins, xmin, xmax));

    addpeak.plotOn(frame2, RooFit.LineColor(kOrange-3), RooFit.LineWidth(4),RooFit.Name("sig_pdf"), RooFit.Normalization(nsig1+ndkk1+ndkpi1,RooAbsReal.NumEvent))

    bkgpdf02.plotOn(frame2, RooFit.LineColor(8), RooFit.LineStyle(6), RooFit.LineWidth(4),RooFit.Name("dkk_pdf"),RooFit.Normalization(ndkk1,RooAbsReal.NumEvent));
    bkgpdf03.plotOn(frame2, RooFit.LineColor(14), RooFit.LineStyle(2),RooFit.LineWidth(4),RooFit.Name("dkpi_pdf"),RooFit.Normalization(ndkpi1,RooAbsReal.NumEvent));
    
    testpdf   = bkgpdf01 .createHistogram("testpdf"  , m2mu, RooFit.Binning(nbins, xmin, xmax));
    testpdf.Scale((nbkg1)/testpdf.Integral(1,-1));
    
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

    dkkLegend = TH2F()
    dkkLegend.SetLineColor(8)
    dkkLegend.SetLineWidth(3)
    dkkLegend.SetLineStyle(6)

    dkpLegend = TH2F()
    dkpLegend.SetLineColor(14)
    dkpLegend.SetLineWidth(3)
    dkpLegend.SetLineStyle(2)

    leg= TLegend(0.55, 0.65,0.9, 0.85)
    leg.SetBorderSize( 0 )
    leg.SetFillStyle( 1001 )
    leg.SetFillColor(kWhite)
    leg.AddEntry( sigLegend, "Signal + Background Fit",  "L" )
    leg.AddEntry( bkgLegend, "Background Only Fit",  "L" )
    leg.AddEntry( datLegend, "Data ("+year+")",  "P" )
    leg.AddEntry( dkkLegend, "D->KK", "L")
    leg.AddEntry( dkpLegend, "D->K#pi","L")
 
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
    latex4.SetTextSize(0.03);
    latex4.SetTextAlign(33);
    h_string4 = "Signal region";
    latex4.DrawLatex(0.5,0.87, h_string4);

    latex5 = TLatex();
    latex5.SetNDC();
    latex5.SetTextSize(0.03);
    latex5.SetTextAlign(33);
    h_string5 = "L < 3.5#sigma_{L}";
    latex5.DrawLatex(0.5,0.82, h_string5);
    
    ratioPad = TPad("BottomPad","",0.,0.08,1.,0.285);
    ratioPad.SetFillStyle(4000);
    ratioPad.SetBottomMargin(0.05);

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


    #Ratio.GetYaxis().SetLabelOffset(0.012);
    Ratio.GetYaxis().SetLabelSize(20);
    Ratio.GetYaxis().SetLabelFont(43);
    Ratio.GetYaxis().SetTitleFont(43);
    Ratio.GetYaxis().SetTitleSize(22);
    Ratio.GetYaxis().SetNdivisions(505);
    Ratio.GetYaxis().SetTitleOffset(2.0);
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
    
    c1.SaveAs("bump_SR_"+year+".pdf");
