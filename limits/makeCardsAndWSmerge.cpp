#include <iostream>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include <THStack.h>
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFractionFitter.h"
#include <string>
#include <vector>
#include <math.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "TF1.h"
#include "TEfficiency.h"

#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <iostream>
#include <valarray>

#include <RooAbsPdf.h> 
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooBernstein.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooMsgService.h> 
#include <RooAddPdf.h> 
#include <TROOT.h> 

//#include "pdfs.h"
//#include <RooDoubleCB.h>


using namespace std;

void makeCardsAndWSmerge(){

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");
  gSystem->Load("$CMSSW_BASE/lib/slc7_amd64_gcc700/libHiggsAnalysisCombinedLimit.so");
  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gSystem->AddIncludePath("-Iinclude/");

TString year[2] = {"2017","2018"};
for(int y = 0; y < 2; y++){ //year        
  //WHICH YEAR
  TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  TFile* file[2] = {NULL, NULL};  // Above and below 3 GeV 


  if (year[y] == "2017"){
    file[0]=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2017.root"); //38.7                                                                                                                         
    file[1]=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2017.root"); //38.7
  }
  else if (year[y] == "2018"){
    file[0]=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2018.root"); //61.3 fb -1
    file[1]=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
  }

  double m[11] 		= {2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   9.0,   10.0,  12.5,  20.0};
  double xSec[11] 	= {8541, 7514, 3323, 2055, 1422, 1043, 793.6, 621.1, 484.3, 292.8, 98.95};//[pb], model-dependent
  TGraph* xsecgraph 	= new TGraph(11,m,xSec);

  //ACCEPTANCE
  TFile* acc_file = TFile::Open("acceptances_dy.root");
  TEfficiency* acc_teff = (TEfficiency*)acc_file->Get("cmsacc");
  int nbins_acc=acc_teff->GetPassedHistogram()->GetNbinsX();
  double acceptances[nbins_acc];
  double m_acceptances[nbins_acc];
  for (int j=1; j<=nbins_acc; j++){
      acceptances[j-1] = acc_teff->GetEfficiency(j);
      m_acceptances[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
  }
  /*TFile* acc_file = TFile::Open("acc_dyturbo.root");
	TH1D* acc_teff = (TH1D*)acc_file->Get("s_m");
	int nbins_acc=acc_teff->GetNbinsX();
	double acceptances[nbins_acc];
	double m_acceptances[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
		acceptances[j-1] = acc_teff->GetBinContent(j);
		cout<<"The acceptence is \n"<<acceptances[j-1]<<endl;		
		m_acceptances[j-1] = acc_teff->GetBinCenter(j);
		}*/
  TGraph* accgraph 	= new TGraph(nbins_acc,m_acceptances,acceptances);
	
  //LUMINOSITY
  double lumi17 = 35300;
  double lumi18 = 61300;
  double luminosity = lumi17+lumi18; 


  //EFFICIENCY
  //get acceptance from hist
  TFile* eff_file_2017 = TFile::Open("modifiedMllEff"+year[0]+".root");
  TFile* eff_file_2018 = TFile::Open("modifiedMllEff"+year[1]+".root");
  //TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_2018_mll_dR_wieghted.root");
  TEfficiency *teff2017 = ((TEfficiency*)eff_file_2017->Get("honemllD_clone"));
  teff2017->Draw();
  teff2017->Paint("");
  TGraphAsymmErrors* effgraph2017 = teff2017->GetPaintedGraph();

  TEfficiency *teff2018 = ((TEfficiency*)eff_file_2018->Get("honemll_clone"));
  teff2018->Draw();
  teff2018->Paint("");
  TGraphAsymmErrors* effgraph2018 = teff2018->GetPaintedGraph();

  //EFFICIENCY SYSTEMATIC
  TFile* trigEffSystFile2017 = TFile::Open("SystEst"+year[0]+".root");
  TH1F *tsys2017 = ((TH1F*)trigEffSystFile2017->Get("triggerSys"));

  TFile* trigEffSystFile2018 = TFile::Open("SystEst"+year[1]+".root");
  TH1F *tsys2018 = ((TH1F*)trigEffSystFile2018->Get("triggerSys"));

  int nbins_tsys=tsys2017->GetNbinsX();

	// int nbins_eff=eff_hist->GetNbinsX();
	// double effA[nbins_eff];
	// double m_effA[nbins_eff];
	// for (int j=1; j<=nbins_eff; j++){
	// 	effA[j-1] = eff_hist->GetBinContent(j);
	// 	m_effA[j-1] = eff_hist->GetBinCenter(j);
	// }
	// TGraph* effgraph 	= new TGraph(nbins_eff,m_effA,effA);
	//double effcuts = 0.64; //from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/329 
  
  effgraph2017->SaveAs("output/effgraph2017.root");
  effgraph2018->SaveAs("output/effgraph2018.root");
  accgraph->SaveAs("output/accF.root");
  xsecgraph->SaveAs("output/xsecgraph.root");
  	
  //scale
  double eps2scale = 1;//0.01;//1;//0.1;//0.002; //this scales eps 0.02 -> 
  
  //*****----->>>>> nSignal = xsecgraph->Eval(mass)*eff*luminosity*acceptances[i]
  
  //define unfittable ranges
  //float unfittable_mins[9] = {0,    0.28, 0.5,   0.7,   0.95, 2.75,  3.55,  8.9};
  
  double unfittable_regions[8][2] = {{0,0.22}, {0.53,0.575}, {0.74,0.85}, {0.97,1.07}, {2.8,3.85}, {9.0,11}};
    


  //ID EFFICIENCY
  TFile* IDfileMVA = NULL;
  IDfileMVA=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_mva.root");	
  TFile* IDfileMVA2 = NULL;
  IDfileMVA2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_mvajpsi.root");	

  TFile* IDfileNO = NULL;
  IDfileNO=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_noid.root");	



   //LOOP OVER MASS INDICES AND MAKE THE CARDS/WORKSPACES
	double mass = -1.;
	TFile* f_ws = TFile::Open(("../mass_calibration/pdfs"+(string)year[0]+".root").c_str(), "READ");
	RooWorkspace *w = (RooWorkspace*)f_ws->Get("dpworkspace");
	w->loadSnapshot("calibrated");
	w->var("alpha1")->setConstant(true);//all this should actually be automatic. Check!!
	w->var("alpha2")->setConstant(true);
	w->var("n1")->setConstant(true);
	w->var("frac_gau")->setConstant(true);
	w->var("gau_reso_scale")->setConstant(true);

	w->Print();

	TGraph* effValues = new TGraph(190);
	TGraph* accValues = new TGraph(190);
	TGraph* plotValues = new TGraph(190);

	TGraph* lar1Values = new TGraph(190);
	TGraph* lar2Values = new TGraph(190);
	TGraph* lar3Values = new TGraph(190);
	TGraph* lar4Values = new TGraph(190);

	double rel_reso=0.013;//temporary

	char errorShape[200];
	sprintf(errorShape,"errorParams%s.txt",year[0].Data());
	ofstream errorparamShape;
	errorparamShape.open(errorShape);


	for(int i=0; i<400; i++){
	  //if (i!=168) continue;
	  	//get the histograms
        TH1D* catA;
        TH1D* catB;

        
        if (i > 290){ 
	  catA=(TH1D*)file[0]->Get(Form("massforLimit_CatA%d",i));
	  catB=(TH1D*)file[0]->Get(Form("massforLimit_CatB%d",i));
        }
	else{
	  catA=(TH1D*)file[1]->Get(Form("massforLimit_CatA%d",i));
	  catB=(TH1D*)file[1]->Get(Form("massforLimit_CatB%d",i));
        }

	//we're using only one category, so we sum all histos
	catA->Add(catB);
	delete catB;
	
		//repeat for both the histograms w/ and w/o the MVA ID
	  	//get the histograms
		TH1D* catAMVA;
                TH1D* catBMVA;
		if (i > 290){
		  catAMVA=(TH1D*)IDfileMVA->Get(Form("massforLimit_CatA%d",i));
		  catBMVA=(TH1D*)IDfileMVA->Get(Form("massforLimit_CatB%d",i));
		} 
		else{
		  catAMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatA%d",i));
                  catBMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatB%d",i));
		}
	  	catAMVA->Add(catBMVA);
	  	delete catBMVA;
		double countMVA = catAMVA->Integral();


	  	TH1D* catANO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatA%d",i));
	  	TH1D* catBNO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatB%d",i));
	  	catANO->Add(catBNO);
	  	delete catBNO;
		double countNO = catANO->Integral();
	
		
	  	double massLow  =  catA->GetXaxis()->GetXmin();
		double massHigh =  catA->GetXaxis()->GetXmax();
		double massBinWidth = massHigh-massLow;
		int nbins = catA->GetNbinsX();


		//compute mass point and define ROOFit variables
	  	bool dontfit=false;

	  	//reduce the mass bin width to avoid being inside the forbidden regions (resonances)
	  	// for (auto fbin : unfittable_regions){
	  	// 	if ((massLow>fbin[0] && massHigh<fbin[1]) || (massLow<fbin[0] && massHigh>fbin[1])) dontfit=true; //the current bin is completely inside a forbidden region, or viceversa
	  	// 	else if ((massHigh-fbin[0])*(massHigh-fbin[1])<=0){ //high edge of our bin is in a forbidden region
	  	// 		massHigh = fbin[0];
	  	// 	}
	  	// 	else if ((massLow-fbin[0])*(massLow-fbin[1])<=0){ //low edge of our bin is in a forbidden region
	  	// 		massLow = fbin[1];
	  	// 	}
	  	// 	if ((massHigh-massLow)<0.2*massBinWidth) dontfit=true; //skip if the mass bin after this reduction is too small (or negative, which would mean the original mass bin was all inside a forbidden region)
	  	// }
	  	// if (dontfit) continue;

	  	mass = 0.5*(massLow+massHigh);
		if (mass < 1.0) continue;
		if (mass >= 8.265) continue;
		if ((mass >= 2.658) && (mass <= 4.16)) continue;
	  	for (auto fbin : unfittable_regions){
	  		if ((mass>fbin[0] && mass<fbin[1])) dontfit=true; //the current point is inside a forbidden region
	  		else if ((massHigh-fbin[0])*(massHigh-fbin[1])<=-0.1){ //high edge of our bin is in a forbidden region
			  //massHigh = fbin[0];
			  //massLow = massHigh-massBinWidth;
			  dontfit=true;
	  		}
	  		else if ((massLow-fbin[0])*(massLow-fbin[1])<=-0.1){ //low edge of our bin is in a forbidden region
			  //massLow = fbin[1];
			  //massHigh = massLow+massBinWidth;
			  dontfit=true;
	  		}
			if ((mass-massLow)<4*rel_reso*mass || (massHigh-mass)<4*rel_reso*mass) dontfit=true; //too close to the bin edge
	  	}
	  	if (dontfit) continue;

		double effcuts = countMVA / countNO;
		//if (mass < 2.0) effcuts = 0.383615;
		if (mass < 2.0) effcuts = 0.05*mass+0.68;
		cout << "The ID efficiency is " << effcuts << " at mass " << mass ; 
		cout << ".  The numerator is " << countMVA << " and the denominator is " << countNO << "\n"; 

        effValues->SetPoint(i, mass, effcuts);
        accValues->SetPoint(i, mass, accgraph->Eval(mass,0,""));
        plotValues->SetPoint(i, mass, effcuts*accgraph->Eval(mass,0,"S")*((lumi17*effgraph2017->Eval(mass,0,"S")+lumi18*effgraph2018->Eval(mass,0,"S"))/(luminosity)));
        
		//Calculate log normal uncertainty for trigger efficiency
		double triggSysVal = tsys2018->GetBinContent(tsys2018->FindBin(mass));
		double triggSys = 1.00 + abs(triggSysVal); 
		cout << ".  The value of trigger syst is " <<  triggSys << "\n"; 

		//cout<<"Spline: "<<effAgraph->Eval(mass,0,"S")<<endl;
		//cout<<"Graph : "<<effAgraph->Eval(mass)<<endl;
		RooRealVar* m2mu = w->var("m2mu");
		m2mu->setMax(massHigh);
		m2mu->setMin(massLow);

		RooAddPdf* signalModel = (RooAddPdf*)w->pdf("signalModel_generic");

		//define the signal model
		w->var("M_generic")->setVal(mass);
		w->var("M_generic")->setConstant(true);
		w->var("res_rel_generic")->setVal(rel_reso);
		w->var("res_rel_generic")->setConstant(true);
		//in pdf.h aggiungi una pdf generica e salvala nel workspace con tutti i param già fissati. poi riprendila da qui, e usa dirett
		// la sua variabile massa osservabile come massa qui, semplicemente cambiandogli il range.

		RooDataHist data_obs("data_obs", "", RooArgList(*m2mu), catA);
		RooRealVar bkg_norm("bkg_norm", "",catA->Integral());


		/*
                RooRealVar par1("par1", "par1", 0.2, 0, 10);
                RooRealVar par2("par2", "par2", 1.5, 0, 10);
                RooRealVar par3("par3", "par3", 2.0, 0, 10);
                RooRealVar par4("par4", "par4", 2.0, 0, 10);
                RooArgList alist(par1, par2, par3, par4);
                RooBernstein bkg_model_bern("bkg_model_bern", "bkg_model_bern", *m2mu, alist);
                //Exponentials
                RooRealVar bar1("bar1", "bar1", -0.5, -7, 7);
                RooExponential bkg_model_exp4("bkg_model_exp4", "bkg_model_exp4", *m2mu, bar1);
                //Product of the two
                RooProdPdf bkg_model("bkg_model", "bkg_model", bkg_model_bern, bkg_model_exp4);
                bkg_model.chi2FitTo(data_obs);
		*/
		/*
                //4th order polynomial
                RooRealVar qar1("qar1", "qar1", 0.0, -10.0, 10.0);
                RooRealVar qar2("qar2", "qar2", 0.0, -10.0, 10.0);
                RooRealVar qar3("qar3", "qar3", 0.0, -10.0, 10.0);
                //RooRealVar qar4("qar4", "qar4", 0.0, -10.0, 10.0);
                RooArgList qlist(qar1, qar2, qar3);
                RooPolynomial bkg_model_line3("bkg_model_line3", "bkg_model_line3", *m2mu, qlist, 1);
                //Breit-Wignet
                RooRealVar war1("war1", "war1", -0.5, -20, 20);
                RooRealVar war2("war2", "war2", -0.5, -20, 20);
                RooBreitWigner bkg_model_BW("bkg_model_BW","bkg_model_BW", *m2mu, war1, war2);
                //Product of the two
                RooProdPdf bkg_model("bkg_model", "bkg_model", bkg_model_line3, bkg_model_BW);
                bkg_model.chi2FitTo(data_obs);
		*/

		

	RooRealVar par1_2017("par1_2017", "par1_2017", 0.2, 0, 10);
	RooRealVar par2_2017("par2_2017", "par2_2017", 1.5, 0, 10);
	RooRealVar par3_2017("par3_2017", "par3_2017", 2.0, 0, 10);
	RooRealVar par4_2017("par4_2017", "par4_2017", 2.0, 0, 10);
	RooArgList alist_2017(par1_2017, par2_2017, par3_2017, par4_2017);
        RooBernstein bern_2017("bern_2017", "bern_2017", *m2mu, alist_2017);

	RooRealVar par1_2018("par1_2018", "par1_2018", 0.2, 0, 10);
	RooRealVar par2_2018("par2_2018", "par2_2018", 1.5, 0, 10);
	RooRealVar par3_2018("par3_2018", "par3_2018", 2.0, 0, 10);
	RooRealVar par4_2018("par4_2018", "par4_2018", 2.0, 0, 10);
	RooArgList alist_2018(par1_2018, par2_2018, par3_2018, par4_2018);
        RooBernstein bern_2018("bern_2018", "bern_2018", *m2mu, alist_2018);

        RooRealVar* M_phi_2017 = w->var("M_phi");
        RooRealVar* res_rel_phi_2017 = w->var("res_rel_phi");
        RooFormulaVar* res_CB_phi_2017 = (RooFormulaVar*)w->function("res_CB_phi");
        RooRealVar* alpha1_2017 = (RooRealVar*)w->var("alpha1");
        RooRealVar* alpha2_2017 = (RooRealVar*)w->var("alpha2");
        RooRealVar* n1_2017 = (RooRealVar*)w->var("n1");
        alpha1_2017->setConstant(false);
        alpha2_2017->setConstant(false);
        n1_2017->setConstant(false);
        n1_2017->setMax(20.0);
        RooRealVar n2_2017("n2_2017","n2_2017",1.0,0.0,20.0);
        RooDoubleCB* signalModel_CB_phi_2017 = new RooDoubleCB("signalModel_CB_phi_2017", "signalModel_CB_phi_2017", *m2mu, *M_phi_2017,*res_CB_phi_2017,*alpha1_2017,*n1_2017,*alpha2_2017,n2_2017);
        RooFormulaVar* res_gau_phi_2017 = (RooFormulaVar*)w->function("res_gau_phi");
        RooGaussian* signalModel_gau_phi_2017 = new RooGaussian("signalModel_gau_phi_2017", "signalModel_gau_phi_2017", *m2mu, *M_phi_2017,*res_gau_phi_2017);
        RooRealVar* frac_gau_2017 = (RooRealVar*)w->var("frac_gau");
        RooAddPdf* phipdf_2017 = new RooAddPdf("phipdf_2017", "phipdf_2017", RooArgList(*signalModel_CB_phi_2017, *signalModel_gau_phi_2017), RooArgList(*frac_gau_2017));

        RooRealVar fbern_2017("fbern_2017", "fbern_2017", 0.5, 0, 1.0);
        RooAddPdf bkg_model_2017("bkg_model_2017","bkg_model_2017",bern_2017,*phipdf_2017,fbern_2017);

        RooRealVar* M_phi_2018 = w->var("M_phi");
        RooRealVar* res_rel_phi_2018 = w->var("res_rel_phi");
        RooFormulaVar* res_CB_phi_2018 = (RooFormulaVar*)w->function("res_CB_phi");
        RooRealVar* alpha1_2018 = (RooRealVar*)w->var("alpha1");
        RooRealVar* alpha2_2018 = (RooRealVar*)w->var("alpha2");
        RooRealVar* n1_2018 = (RooRealVar*)w->var("n1");
        alpha1_2018->setConstant(false);
        alpha2_2018->setConstant(false);
        n1_2018->setConstant(false);
        n1_2018->setMax(20.0);
        RooRealVar n2_2018("n2_2018","n2_2018",1.0,0.0,20.0);
        RooDoubleCB* signalModel_CB_phi_2018 = new RooDoubleCB("signalModel_CB_phi_2018", "signalModel_CB_phi_2018", *m2mu, *M_phi_2018,*res_CB_phi_2018,*alpha1_2018,*n1_2018,*alpha2_2018,n2_2018);
        RooFormulaVar* res_gau_phi_2018 = (RooFormulaVar*)w->function("res_gau_phi");
        RooGaussian* signalModel_gau_phi_2018 = new RooGaussian("signalModel_gau_phi_2018", "signalModel_gau_phi_2018", *m2mu, *M_phi_2018,*res_gau_phi_2018);
        RooRealVar* frac_gau_2018 = (RooRealVar*)w->var("frac_gau");
        RooAddPdf* phipdf_2018 = new RooAddPdf("phipdf_2018", "phipdf_2018", RooArgList(*signalModel_CB_phi_2018, *signalModel_gau_phi_2018), RooArgList(*frac_gau_2018));

        RooRealVar fbern_2018("fbern_2018", "fbern_2018", 0.5, 0, 1.0);
        RooAddPdf bkg_model_2018("bkg_model_2018","bkg_model_2018",bern_2018,*phipdf_2018,fbern_2018);

        if (i>176) {
	  //fbern_2017.setVal(0.999);
          //fbern_2017.setConstant(true);
	  //fbern_2018.setVal(0.999);
	  //fbern_2018.setConstant(true);
        }        
	
	if (year[y] == "2017"){
	  bkg_model_2017.fitTo(data_obs);		
	}
	else if (year[y] == "2018"){
	  bkg_model_2018.fitTo(data_obs);		
	}

        alpha1_2017->setConstant(true);
        alpha2_2017->setConstant(true);
        n1_2017->setConstant(true);
        M_phi_2017->setConstant(true);
        res_rel_phi_2017->setConstant(true);

        alpha1_2018->setConstant(true);
        alpha2_2018->setConstant(true);
        n1_2018->setConstant(true);
        M_phi_2018->setConstant(true);
        res_rel_phi_2018->setConstant(true);


	/*
		RooPlot *frame = m2mu->frame();
		data_obs.plotOn(frame);
        bkg_modelb.plotOn(frame,RooFit::LineColor(4));
		TCanvas c_all("c_all", "c_all", 800, 500);
		frame->Draw("goff");
		c_all.SaveAs(Form("output/catA_%d_Merged.png",i));
	*/
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	gStyle->SetPadTopMargin(0.07);
	gStyle->SetPadBottomMargin(0.3);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.07);

	gStyle->SetNdivisions(508, "X");
	gStyle->SetNdivisions(508, "Y");

	TCanvas c_all("c_all", "c_all", 800, 800);
	RooPlot *frame = m2mu->frame();
	frame->SetYTitle("Events");
	frame->SetXTitle("Dimuon Mass (GeV)");
	frame->GetXaxis()->SetLabelSize(0.025);
	frame->GetYaxis()->SetLabelSize(0.025);

	data_obs.plotOn(frame);
	if (year[y] == "2017"){
	  bkg_model_2017.plotOn(frame, RooFit::Name("pdf1"), LineColor(2));
	}else if (year[y] == "2018"){
	  bkg_model_2018.plotOn(frame, RooFit::Name("pdf1"), LineColor(2));
	}


	frame->Draw("goff");
	// Binned instances of the pdf and data for chisq calculation                                                                                                                                                                 
	TH1 *testhisto = data_obs.createHistogram("testhisto", *m2mu, RooFit::Binning(nbins, massLow, massHigh));
	TH1 *testpdf1;
	if (year[y] == "2017"){
	  testpdf1   = bkg_model_2017.createHistogram("testpdf1"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
	}else if (year[y] == "2018"){
	  testpdf1   = bkg_model_2018.createHistogram("testpdf1"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
	}
	testpdf1->Scale(testhisto->Integral(1,-1)/testpdf1->Integral(1,-1));

	int noofparm1 = 5;
	Double_t mychi21 = 0.;

	for(int i=1; i<=testhisto->GetNbinsX(); i++){
	  Double_t tmp_mychi21 = ((testpdf1->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf1->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf1->GetBinContent(i);
	  cout << "testpdf1 increment " << testpdf1->GetBinContent(i) << "\n";
	  mychi21 += tmp_mychi21;
	}
	Double_t mychi2_final1 = mychi21/(nbins - (noofparm1 + 1));
	Double_t RooPlot_chi21 = frame->chiSquare("new_pdf", "new_hist", noofparm1 + 1);

	TLatex latex1;
	latex1.SetNDC();
	latex1.SetTextSize(0.04);
	latex1.SetTextAlign(33);
	TString h_string1 = "Chi2/ndf = " + std::to_string(mychi2_final1);
	h_string1 = "#chi^{2}/(N_{bins}-N_{d.o.f}+1) =" + std::to_string(mychi2_final1);
	latex1.DrawLatex(0.92,0.980, h_string1);

	// Set up lower frame, to display the pull                                                                                                                                                                                    
	TH2F *hframe2= new TH2F("hframe2","hframe2",500, massLow, massHigh, 500, -4, 4);
	hframe2->SetYTitle("(Data-Fit)/Unc.");
	hframe2->GetXaxis()->SetLabelOffset(1);
	hframe2->GetXaxis()->SetLabelSize(0.1);
	hframe2->GetYaxis()->SetLabelOffset(0.012);
	hframe2->GetYaxis()->SetLabelSize(0.1);
	hframe2->GetYaxis()->SetTitleOffset(0.27);
	hframe2->GetYaxis()->SetTitleSize(0.15);

	// Set up and draw on lower pad                                                                                                                                                                                               
	TPad *ratioPad = new TPad("BottomPad","",0.,0.03,1.,0.23);
	ratioPad->SetBottomMargin(2.1);
	ratioPad->Draw();
	ratioPad->cd();

	hframe2->Draw("");

	TLine *line1 = new TLine(massLow, 0, massHigh, 0);
	line1->SetLineColor(kBlue);
	line1->SetLineWidth(1);
	line1->Draw("same");

	TH1F * Ratio1 = ((TH1F*)testhisto->Clone("Ratio1"));
	for(int i = 1; i<testhisto->GetNbinsX(); i++){
	  Ratio1->SetBinContent(i, (testhisto->GetBinContent(i) - testpdf1->GetBinContent(i))/testhisto->GetBinError(i));
	}
	Ratio1->SetFillColor(2);
	Ratio1->SetLineColor(2);
	Ratio1->Draw("histsame");
	c_all.SaveAs(Form("output/catA_%d_"+year[y]+".png",i));

		//save into ROO workspace
		RooWorkspace dpworkspace("dpworkspace", "");
		dpworkspace.import(data_obs);
		dpworkspace.import(*signalModel);
		if (year[y] == "2017"){
		  dpworkspace.import(bkg_model_2017);
		} else if (year[y] == "2018"){
		  dpworkspace.import(bkg_model_2018);
		}

        RooRealVar bkg_model_norm("bkg_model_norm","Number of background events",catA->Integral(),0.8*catA->Integral(),1.2*catA->Integral());
        dpworkspace.import(bkg_model_norm);
		dpworkspace.writeToFile(Form("output/dpWorkspace_Merged"+suff+"_%d.root",i));

		//write the datacard
		char inputShape[200];
		sprintf(inputShape,"output/dpCard_Merged"+suff+"_m%.3f_%d.txt",mass,i);
		ofstream newcardShape;
		newcardShape.open(inputShape);
		newcardShape << Form("imax * number of channels\n");
		newcardShape << Form("jmax * number of background\n");
		newcardShape << Form("kmax * number of nuisance parameters\n");
		newcardShape << Form("shapes data_obs   CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:data_obs\n",i);
		newcardShape << Form("shapes bkg_mass   CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:bkg_model_"+year[y]+"\n",i);
                newcardShape << Form("shapes signalModel_generic        CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:signalModel_generic\n",i);
		newcardShape << Form("bin		CatAB\n");
		newcardShape << Form("observation 	-1.0\n");
		newcardShape << Form("bin     		CatAB		CatAB		\n");
		newcardShape << Form("process 		signalModel_generic  	bkg_mass	\n");
		newcardShape << Form("process 		0    		1	   	\n");
		newcardShape << Form("rate    		%f  		%f		\n",
                             //effcuts*((lumi17*effgraph2017->Eval(mass,0,"S")+lumi18*effgraph2018->Eval(mass,0,"S"))/(luminosity))*luminosity, catA->Integral());
                             effcuts*((lumi17*effgraph2017->Eval(mass,0,"S")+lumi18*effgraph2018->Eval(mass,0,"S"))/(luminosity))*luminosity, 1.0);
		newcardShape << Form("lumi13TeV lnN 	1.026 	-\n");
		newcardShape << Form("id_eff_mva lnN	1.10 	-\n");
		newcardShape << Form("eff_trig lnN         %f        -\n", triggSys);
		newcardShape << Form("pdf_index discrete \n");
		//newcardShape << Form("eff_mu_13TeV_2017 lnN	1.015 	-\n");
		//newcardShape << Form("bkg_norm rateParam CatA bkg_mass %f\n",catA->Integral());
		//newcardShape << Form("resA param %f %f\n",resA.getValV(),resA.getValV()*0.1);
		newcardShape.close();
		/*
		double par1val = par1.getValV();
		double par2val = par2.getValV();
		double par3val = par3.getValV();
		double par4val = par4.getValV();
		double par1err = par1.getError();
		double par2err = par2.getError();
		double par3err = par3.getError();
		double par4err = par4.getError();

                //write the error params
                errorparamShape << Form("massPoint%.3f --setParameterRanges par1=%f,%f:", mass, par1val-5.0*par1err, par1val+5.0*par1err);
                errorparamShape << Form("par2=%f,%f:", par2val-5.0*par2err, par2val+5.0*par2err);
                errorparamShape << Form("par3=%f,%f:", par3val-5.0*par3err, par3val+5.0*par3err);
                errorparamShape << Form("par4=%f,%f\n", par4val-5.0*par4err, par4val+5.0*par4err);
		*/

	}
	errorparamShape.close();
	f_ws->Close();

	/*
	TCanvas c_fVal("c_fVal", "c_fVal", 950, 1020);
	//effValues->GetYaxis()->SetRangeUser(0.00001, 1);
	effValues->GetXaxis()->SetRangeUser(0.8, 9);
	effValues->GetXaxis()->SetTitle("m_{#mu#mu}");
	effValues->SetTitle("Efficiency (Muon MVA + PVd Cut)");
        effValues->Draw("a*");
	TMarker *point = new TMarker(1,0.73,8);
	point->SetMarkerColor(4);
	point->Draw("s");
	TLine *l = new TLine(2,0,2,1);
	//TLine *l = new TLine(2,0,2,c_fVal.GetY1());
	l->SetLineColor(2);
	l->Draw("s");
	auto legend = new TLegend(0.6,0.1,0.9,0.4);
	legend->AddEntry(effValues,"ID Efficiency","p");
	legend->AddEntry(point,"Scalar sample point","p");
	legend->AddEntry(l,"Extrapolation threshold","l");
	
	legend->Draw();
	//c_fVal.SetLogy();
        c_fVal.SaveAs("ID_EffBareDistribution.png");
	*/
	/*
	TCanvas c_fVal("c_fVal", "c_fVal", 1250, 1020);
        //effValues->GetYaxis()->SetRangeUser(0.00001, 1);
        lar1Values->GetXaxis()->SetRangeUser(0.8, 9);
        lar1Values->GetYaxis()->SetRangeUser(-5, 5);
        lar1Values->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
        c_fVal.cd(4);
	lar1Values->SetMarkerColor(2);
	lar1Values->SetMarkerStyle(18);
        lar1Values->Draw("ap");
	lar1Values->SetTitle("");

	lar2Values->SetMarkerColor(3);
	lar2Values->SetMarkerStyle(18);
        lar2Values->Draw("SP");
	lar3Values->SetMarkerColor(4);
	lar3Values->SetMarkerStyle(18);
        lar3Values->Draw("SP");
	lar4Values->SetMarkerColor(6);
	lar4Values->SetMarkerStyle(18);
        lar4Values->Draw("SP");
	//plotValues->Draw("aSAME");
	//auto legend = new TLegend(0.6,0.1,0.9,0.4);
        //legend->AddEntry(accValues,"Acceptance","p");
        //legend->AddEntry(plotValues,"Efficiency (ID + Trigger) #times Acceptance","p");
        //legend->Draw();
	auto cmsTag= new TLatex(0.13,0.917,"#color[2]{#scale[1.1]{lar1}} #color[3]{#scale[1.1]{lar2}} #color[4]{#scale[1.1]{lar3}} #color[6]{#scale[1.1]{lar4}}");
	cmsTag->SetNDC();
	cmsTag->SetTextAlign(11);
	cmsTag->Draw();
	//auto cmsTag2 = new TLatex(0.215,0.917,"#scale[0.825]{#bf{#it{Preliminary}}}");
	//cmsTag2->SetNDC();
	//cmsTag2->SetTextAlign(11);
	//cmsTag2->Draw();
        //c_fVal.SetLogy();
        c_fVal.SaveAs("carplot.png");


	for (int j=1; j<=nbins_tsys; j++){
	  double val = 1.00 + abs(tsys2018->GetBinContent(j));
	  cout <<  "Tris sys " << val << "\n"; 
	}
	*/
  }
	
}
