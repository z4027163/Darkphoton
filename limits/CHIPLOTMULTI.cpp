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
//#include "../../../HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

using namespace std;

void CHIPLOTMULTI(){

  TString year[2] = {"2017","2018"};
  for(int y = 0; y < 2; y++){ //year
  

  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  	TFile* file = NULL;  // Above 3 GeV
  	TFile* file2 = NULL; // Below 3 GeV
	//if (year == "2017") file=TFile::Open("/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2017/ScoutingRunD/mergedHistos_v1.root");
        if (year[y] == "2017"){
          file=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2017.root"); //38.7
          file2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2017.root"); //38.7
        }
        else if (year[y] == "2018"){
          file=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2018.root"); //61.3 fb -1
          file2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
        }
  //PREPARE EXPECTED NUMBER OF SIGNAL EVENTS PER CATEGORY
	//X-SECTION GRAPH
	//double m[17] 		= {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 14.0, 16.0, 18.0, 20.0};
	//double xSec[17] 	= {10,  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  10,  10,  10,  10,  10};//[pb]
	double m[11] 		= {2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   9.0,   10.0,  12.5,  20.0};
	double xSec[11] 	= {8541, 7514, 3323, 2055, 1422, 1043, 793.6, 621.1, 484.3, 292.8, 98.95};//[pb], model-dependent
	//TGraph* xsecgraph 	= new TGraph(17,m,xSec);
	TGraph* xsecgraph 	= new TGraph(11,m,xSec);



	//SECLECTION UNCERTAINTY
	TFile* uncFile = TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/unc/idunc_mass.root");
        TH1F* selUncH      =     NULL;
	if (year[y] == "2017"){
          selUncH = (TH1F*)uncFile->Get("unc2017");
        } else if (year[y] == "2018"){
          selUncH = (TH1F*)uncFile->Get("unc2017");
        }



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
	
	//TF1* accF = (TF1*)acc_file->Get("fit_func");
	//LUMINOSITY
	double luminosity = 0; //4000.;//pb-1
	if (year[y] == "2017") luminosity = 35300;
	//if (year[y] == "2017") luminosity = 4000;

	//else if (year[y] == "2018") luminosity = 6600;
	else if (year[y] == "2018") luminosity = 61300;
	//else if (year[y] == "2018") luminosity = 1183;

	//EFFICIENCY
	//get acceptance from hist
	//TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_"+year[y]+"_mll_dR_wieghted.root");
	TFile* eff_file = TFile::Open("modifiedMllEff"+year[y]+".root");
	//TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_2018_mll_dR_wieghted.root");
	TEfficiency *teff;
	if (year[y] == "2017"){
	  teff = ((TEfficiency*)eff_file->Get("honemllD_clone"));
	}
	else if (year[y] == "2018"){
	  teff = ((TEfficiency*)eff_file->Get("honemll_clone"));
	}

	//cout<<teff<<endl;
	teff->Draw();
	teff->Paint("");
	TGraphAsymmErrors* effgraph = teff->GetPaintedGraph();
	//EFFICIENCY SYSTEMATIC
	//get estimated efficiency from list
	TFile* trigEffSystFile = TFile::Open("SystEst"+year[y]+".root");
	TH1F *tsys = ((TH1F*)trigEffSystFile->Get("triggerSys"));
	int nbins_tsys=tsys->GetNbinsX();

	// int nbins_eff=eff_hist->GetNbinsX();
	// double effA[nbins_eff];
	// double m_effA[nbins_eff];
	// for (int j=1; j<=nbins_eff; j++){
	// 	effA[j-1] = eff_hist->GetBinContent(j);
	// 	m_effA[j-1] = eff_hist->GetBinCenter(j);
	// }
	// TGraph* effgraph 	= new TGraph(nbins_eff,m_effA,effA);
	//double effcuts = 0.64; //from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/329 

	effgraph->SaveAs("output/effgraph.root");
	accgraph->SaveAs("output/accF.root");
	xsecgraph->SaveAs("output/xsecgraph.root");

	
	//scale
	double eps2scale = 1;//0.01;//1;//0.1;//0.002; //this scales eps 0.02 -> 

	//*****----->>>>> nSignal = xsecgraph->Eval(mass)*eff*luminosity*acceptances[i]

	//define unfittable ranges
	//float unfittable_mins[9] = {0,    0.28, 0.5,   0.7,   0.95, 2.75,  3.55,  8.9};

	double unfittable_regions[8][2] = {{0,0.22}, {0.53,0.575}, {0.74,0.85}, {0.97,1.12}, {2.8,3.85}, {9.0,11}};



	//ID EFFICIENCY
	TFile* IDfileMVA = NULL;
	IDfileMVA=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_mva.root");	
	TFile* IDfileMVA2 = NULL;
	IDfileMVA2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_mvajpsi.root");	

	TFile* IDfileNO = NULL;
	IDfileNO=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/lowDY/lowDY_noid.root");	



   //LOOP OVER MASS INDICES AND MAKE THE CARDS/WORKSPACES
	double mass = -1.;
	double rel_reso=0.013;//temporary
	TFile* f_ws = TFile::Open(("../mass_calibration/pdfs"+(string)year[y]+".root").c_str(), "READ");
	RooWorkspace *w = (RooWorkspace*)f_ws->Get("dpworkspace");
	w->loadSnapshot("calibrated");
	w->var("alpha1")->setConstant(true);//all this should actually be automatic. Check!!
	w->var("alpha2")->setConstant(true);
	w->var("n1")->setConstant(true);
	w->var("frac_gau")->setConstant(true);
	w->var("gau_reso_scale")->setConstant(true);
	w->var("res_rel_unc")->setConstant(true);
    
	w->Print();


	RooWorkspace *upsilon = (RooWorkspace*)f_ws->Get("dpworkspace");
        upsilon->loadSnapshot("calibrated");
        upsilon->var("alpha1")->setConstant(true);//all this should actually be automatic. Check!!
        upsilon->var("alpha2")->setConstant(true);
        upsilon->var("n1")->setConstant(true);
        upsilon->var("frac_gau")->setConstant(true);
        upsilon->var("gau_reso_scale")->setConstant(true);
	upsilon->var("M_generic")->setVal(9.46);
	upsilon->var("M_generic")->setConstant(true);
	upsilon->var("res_rel_generic")->setVal(rel_reso);
	upsilon->var("res_rel_generic")->setConstant(true);

        TGraph* effValues = new TGraph(190);
        TGraph* accValues = new TGraph(190);
        TGraph* plotValues = new TGraph(190);


        TGraph* chi1Values = new TGraph(190);
        TGraph* chi2Values = new TGraph(190);
        TGraph* chi3Values = new TGraph(190);
        TGraph* chi4Values = new TGraph(190);
        TGraph* chi5Values = new TGraph(190);
        TGraph* chi6Values = new TGraph(190);


	for(int i=0; i<400; i++){
	  	//get the histograms
                if (i<177) continue;
		if (i>368) continue;
	        TH1D* catA;
	        TH1D* catB;
	        if (i > 290){ 
		  catA=(TH1D*)file->Get(Form("massforLimit_CatA%d",i));
		  catB=(TH1D*)file->Get(Form("massforLimit_CatB%d",i));
	        }
		else{
		  catA=(TH1D*)file2->Get(Form("massforLimit_CatA%d",i));
		  catB=(TH1D*)file2->Get(Form("massforLimit_CatB%d",i));
	        }
	  	//TH1D* catC=(TH1D*)file->Get(Form("massforLimit_CatC%d",i));

	  	//we're using only one category, so we sum all histos
	  	catA->Add(catB);
	  	//catA->Add(catC);
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
		catAMVA->Rebin(2);
	  	delete catBMVA;
		double countMVA = catAMVA->Integral();


	  	TH1D* catANO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatA%d",i));
	  	TH1D* catBNO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatB%d",i));
	  	catANO->Add(catBNO);
		catANO->Rebin(2);
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
		if ((mass >= 2.6) && (mass <= 4.16)) continue;
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
                plotValues->SetPoint(i, mass, effcuts*accgraph->Eval(mass,0,"S")*effgraph->Eval(mass,0,"S"));

		//Calculate log normal uncertainty for trigger efficiency
		double triggSysVal = tsys->GetBinContent(tsys->FindBin(mass));
		double triggSys = 1.00 + abs(triggSysVal); 
		cout << ".  The value of trigger syst is " <<  triggSys << "\n"; 


                double selSysVal = selUncH->GetBinContent(selUncH->FindBin(mass));
		double selSys = 1.00 + abs(selSysVal);
		cout << ".  The value of sel syst is " <<  selSys << "\n";



		//cout<<"Spline: "<<effAgraph->Eval(mass,0,"S")<<endl;
		//cout<<"Graph : "<<effAgraph->Eval(mass)<<endl;
		RooRealVar* m2mu = w->var("m2mu");
		m2mu->setMax(massHigh);
		m2mu->setMin(massLow);

		RooAddPdf* signalModel = (RooAddPdf*)w->pdf("signalModel_generic");


                RooAddPdf* upsilonBG = (RooAddPdf*)upsilon->pdf("signalModel_generic");

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
                RooRealVar lar1_2017("lar1_2017", "lar1_2017", 0.0, -20.0, 20.0);
                RooRealVar lar2_2017("lar2_2017", "lar2_2017", 0.0, -20.0, 20.0);
                RooRealVar lar3_2017("lar3_2017", "lar3_2017", 0.0, -20.0, 20.0);
                RooRealVar lar4_2017("lar4_2017", "lar4_2017", 0.0, -20.0, 20.0);
                RooArgList llist_2017(lar1_2017, lar2_2017, lar3_2017, lar4_2017);
                RooPolynomial bkg_model_line4_2017("bkg_model_line_2017", "bkg_model_line_2017", *m2mu, llist_2017, 1);
                //Exponentials
                RooRealVar car1_2017("car1_2017", "car1_2017", -0.5, -20, 20);
                RooExponential bkg_model_exp4_2017("bkg_model_exp4_2017", "bkg_model_exp4_2017", *m2mu, car1_2017);
                //Product of the two
                RooProdPdf bkg_model_pol4xexp_2017("bkg_model_pol4xexp_2017", "bkg_model_pol4xexp_2017", bkg_model_line4_2017, bkg_model_exp4_2017);
                bkg_model_pol4xexp_2017.chi2FitTo(data_obs);

                RooRealVar lar1_2018("lar1_2018", "lar1_2018", 0.0, -20.0, 20.0);
                RooRealVar lar2_2018("lar2_2018", "lar2_2018", 0.0, -20.0, 20.0);
                RooRealVar lar3_2018("lar3_2018", "lar3_2018", 0.0, -20.0, 20.0);
                RooRealVar lar4_2018("lar4_2018", "lar4_2018", 0.0, -20.0, 20.0);
                RooArgList llist_2018(lar1_2018, lar2_2018, lar3_2018, lar4_2018);
                RooPolynomial bkg_model_line4_2018("bkg_model_line_2018", "bkg_model_line_2018", *m2mu, llist_2018, 1);
                //Exponentials
                RooRealVar car1_2018("car1_2018", "car1_2018", -0.5, -20, 20);
                RooExponential bkg_model_exp4_2018("bkg_model_exp4_2018", "bkg_model_exp4_2018", *m2mu, car1_2018);
                //Product of the two
                RooProdPdf bkg_model_pol4xexp_2018("bkg_model_pol4xexp_2018", "bkg_model_pol4xexp_2018", bkg_model_line4_2018, bkg_model_exp4_2018);
                bkg_model_pol4xexp_2018.chi2FitTo(data_obs);
		*/

                RooRealVar lar1_2017("lar1_2017", "lar1_2017",0.5, 0., 1.);
                RooRealVar lar2_2017("lar2_2017", "lar2_2017",-0.2, -1, 1.);
                RooRealVar lar3_2017("lar3_2017", "lar3_2017",-0.2, -1, 1.);
                RooRealVar lar4_2017("lar4_2017", "lar4_2017",-0.2, -1, 1.);
                RooChebychev bkg_model_pol4xexp_2017("bkg_model_pol4xexp_2017", "bkg_model_pol4xexp_2017", *m2mu, RooArgSet(lar1_2017, lar2_2017, lar3_2017, lar4_2017));
                bkg_model_pol4xexp_2017.chi2FitTo(data_obs);

                RooRealVar lar1_2018("lar1_2018", "lar1_2018",0.5, 0., 1.);
                RooRealVar lar2_2018("lar2_2018", "lar2_2018",-0.2, -1, 1.);
                RooRealVar lar3_2018("lar3_2018", "lar3_2018",-0.2, -1, 1.);
                RooRealVar lar4_2018("lar4_2018", "lar4_2018",-0.2, -1, 1.);
                RooChebychev bkg_model_pol4xexp_2018("bkg_model_pol4xexp_2018", "bkg_model_pol4xexp_2018", *m2mu, RooArgSet(lar1_2018, lar2_2018, lar3_2018, lar4_2018));
		bkg_model_pol4xexp_2018.chi2FitTo(data_obs);


		
		
		/*
		RooRealVar lar1_2017("war1_2017", "war1_2017", 0.2, 0, 10);
		RooRealVar war2_2017("war2_2017", "war2_2017", 1.5, 0, 10);
		RooRealVar war3_2017("war3_2017", "war3_2017", 2.0, 0, 10);
		RooArgList wlist_2017(war1_2017, war2_2017, war3_2017);
		RooBernstein bkg_model_bern3_2017("bkg_model_bern3_2017", "bkg_model_bern3_2017", *m2mu, wlist_2017);
		bkg_model_bern3_2017.fitTo(data_obs);		


		RooRealVar war1_2018("war1_2018", "war1_2018", 0.2, 0, 10);
		RooRealVar war2_2018("war2_2018", "war2_2018", 1.5, 0, 10);
		RooRealVar war3_2018("war3_2018", "war3_2018", 2.0, 0, 10);
		RooArgList wlist_2018(war1_2018, war2_2018, war3_2018);
		RooBernstein bkg_model_bern3_2018("bkg_model_bern3_2018", "bkg_model_bern3_2018", *m2mu, wlist_2018);
		bkg_model_bern3_2018.fitTo(data_obs);
		*/
		
		
		RooRealVar par1_2017("par1_2017", "par1_2017", 0.2, 0, 20);
		RooRealVar par2_2017("par2_2017", "par2_2017", 1.5, 0, 20);
		RooRealVar par3_2017("par3_2017", "par3_2017", 2.0, 0, 20);
		RooRealVar par4_2017("par4_2017", "par4_2017", 2.0, 0, 20);
		RooArgList alist_2017(par1_2017, par2_2017, par3_2017, par4_2017);
		RooBernstein bkg_model_bern4_2017("bkg_model_bern4_2017", "bkg_model_bern4_2017", *m2mu, alist_2017);
		bkg_model_bern4_2017.fitTo(data_obs);		
		

		RooRealVar par1_2018("par1_2018", "par1_2018", 0.2, 0, 20);
		RooRealVar par2_2018("par2_2018", "par2_2018", 1.5, 0, 20);
		RooRealVar par3_2018("par3_2018", "par3_2018", 2.0, 0, 20);
		RooRealVar par4_2018("par4_2018", "par4_2018", 2.0, 0, 20);
		RooArgList alist_2018(par1_2018, par2_2018, par3_2018, par4_2018);
		RooBernstein bkg_model_bern4_2018("bkg_model_bern4_2018", "bkg_model_bern4_2018", *m2mu, alist_2018);
		bkg_model_bern4_2018.fitTo(data_obs);
		
		/*
		
		RooRealVar war1_2017("war1_2017", "war1_2017", 0.2, 0, 20);
		RooRealVar war2_2017("war2_2017", "war2_2017", 1.5, 0, 20);
		RooRealVar war3_2017("war3_2017", "war3_2017", 2.0, 0, 20);
		RooRealVar war4_2017("war4_2017", "war4_2017", 2.0, 0, 20);
		RooRealVar war5_2017("war5_2017", "war5_2017", 2.0, 0, 20);
		RooArgList wlist_2017(war1_2017, war2_2017, war3_2017, war4_2017, war5_2017);
		RooBernstein bkg_model_bern5_2017("bkg_model_bern5_2017", "bkg_model_bern5_2017", *m2mu, wlist_2017);
		bkg_model_bern5_2017.fitTo(data_obs);		
		RooRealVar war1_2018("war1_2018", "war1_2018", 0.2, 0, 20);
		RooRealVar war2_2018("war2_2018", "war2_2018", 1.5, 0, 20);
		RooRealVar war3_2018("war3_2018", "war3_2018", 2.0, 0, 20);
		RooRealVar war4_2018("war4_2018", "war4_2018", 2.0, 0, 20);
		RooRealVar war5_2018("war5_2018", "war5_2018", 2.0, 0, 20);
		RooArgList wlist_2018(war1_2018, war2_2018, war3_2018, war4_2018, war5_2018);
		RooBernstein bkg_model_bern5_2018("bkg_model_bern5_2018", "bkg_model_bern5_2018", *m2mu, wlist_2018);
		bkg_model_bern5_2018.fitTo(data_obs);		
		
		RooRealVar dar1_2017("dar1_2017", "dar1_2017", 0.2, 0, 20);
		RooRealVar dar2_2017("dar2_2017", "dar2_2017", 1.5, 0, 20);
		RooRealVar dar3_2017("dar3_2017", "dar3_2017", 2.0, 0, 20);
		RooRealVar dar4_2017("dar4_2017", "dar4_2017", 2.0, 0, 20);
		RooRealVar dar5_2017("dar5_2017", "dar5_2017", 2.0, 0, 20);
		RooRealVar dar6_2017("dar6_2017", "dar6_2017", 2.0, 0, 20);
		RooArgList dlist_2017(dar1_2017, dar2_2017, dar3_2017, dar4_2017, dar5_2017, dar6_2017);
		RooBernstein bkg_model_bern6_2017("bkg_model_bern6_2017", "bkg_model_bern6_2017", *m2mu, dlist_2017);
		bkg_model_bern6_2017.fitTo(data_obs);		
		
		RooRealVar dar1_2018("dar1_2018", "dar1_2018", 0.2, 0, 20);
		RooRealVar dar2_2018("dar2_2018", "dar2_2018", 1.5, 0, 20);
		RooRealVar dar3_2018("dar3_2018", "dar3_2018", 2.0, 0, 20);
		RooRealVar dar4_2018("dar4_2018", "dar4_2018", 2.0, 0, 20);
		RooRealVar dar5_2018("dar5_2018", "dar5_2018", 2.0, 0, 20);
		RooRealVar dar6_2018("dar6_2018", "dar6_2018", 2.0, 0, 20);
		RooArgList dlist_2018(dar1_2018, dar2_2018, dar3_2018, dar4_2018, dar5_2018, dar6_2018);
		RooBernstein bkg_model_bern6_2018("bkg_model_bern6_2018", "bkg_model_bern6_2018", *m2mu, dlist_2018);
		bkg_model_bern6_2018.fitTo(data_obs);		
		*/

		
		
		RooRealVar bar1_2017("bar1_2017", "bar1_2017", -0.5, -20, 20);                            
		RooRealVar bf1_2017("bf1_2017","bf1_2017",0.2,0.0,1.0);   
		RooExponential exp1_2017("exp1_2017", "exp1_2017", *m2mu, bar1_2017);
		RooRealVar bar2_2017("bar2_2017", "bar2_2017", -0.5, -20, 20);                                                                                          
		RooRealVar bf2_2017("bf2_2017","bf2_2017",0.2,0.0,1.0);                              
		RooExponential exp2_2017("exp2_2017", "exp2_2017", *m2mu, bar2_2017);
		RooRealVar bar3_2017("bar3_2017", "bar3_2017", -0.5, -20, 20);
		RooRealVar bf3_2017("bf3_2017","bf3_2017",0.2,0.0,1.0);
		RooExponential exp3_2017("exp3_2017", "exp3_2017", *m2mu, bar3_2017);
		RooRealVar bar4_2017("bar4_2017", "bar4_2017", -0.5, -20, 20);
		RooRealVar bf4_2017("bf4_2017","bf4_2017",0.2,0.0,1.0);
		RooExponential exp4_2017("exp4_2017", "exp4_2017", *m2mu, bar4_2017);
		RooRealVar bar5_2017("bar5_2017", "bar5_2017", -0.5, -20, 20);
		RooRealVar bf5_2017("bf5_2017","bf5_2017",0.2,0.0,1.0);
		RooExponential exp5_2017("exp5_2017", "exp5_2017", *m2mu, bar5_2017);
                RooRealVar bar6_2017("bar6_2017", "bar6_2017", -0.5, -20, 20);
                RooRealVar bf6_2017("bf6_2017","bf6_2017",0.2,0.0,1.0);
                RooExponential exp6_2017("exp6_2017", "exp6_2017", *m2mu, bar6_2017);
                RooRealVar bar7_2017("bar7_2017", "bar7_2017", -0.5, -20, 20);
                RooExponential exp7_2017("exp7_2017", "exp7_2017", *m2mu, bar7_2017);

		RooArgList explist_2017(exp1_2017,exp2_2017,exp3_2017,exp4_2017,exp5_2017,exp6_2017,exp7_2017);
		RooArgList expclist_2017(bf1_2017,bf2_2017,bf3_2017,bf4_2017,bf5_2017,bf6_2017);
		RooAddPdf bkg_model_exp7_2017("bkg_model_exp7_2017","bkg_model_exp7_2017",explist_2017,expclist_2017,true);
		bkg_model_exp7_2017.fitTo(data_obs);


		RooRealVar bar1_2018("bar1_2018", "bar1_2018", -0.5, -20, 20);                            
		RooRealVar bf1_2018("bf1_2018","bf1_2018",0.2,0.0,1.0);   
		RooExponential exp1_2018("exp1_2018", "exp1_2018", *m2mu, bar1_2018);
		RooRealVar bar2_2018("bar2_2018", "bar2_2018", -0.5, -20, 20);                                                                                           
		RooRealVar bf2_2018("bf2_2018","bf2_2018",0.2,0.0,1.0);                                                                                                  
		RooExponential exp2_2018("exp2_2018", "exp2_2018", *m2mu, bar2_2018);
		RooRealVar bar3_2018("bar3_2018", "bar3_2018", -0.5, -20, 20);
		RooRealVar bf3_2018("bf3_2018","bf3_2018",0.2,0.0,1.0);
		RooExponential exp3_2018("exp3_2018", "exp3_2018", *m2mu, bar3_2018);
		RooRealVar bar4_2018("bar4_2018", "bar4_2018", -0.5, -20, 20);
		RooRealVar bf4_2018("bf4_2018","bf4_2018",0.2,0.0,1.0);
		RooExponential exp4_2018("exp4_2018", "exp4_2018", *m2mu, bar4_2018);
		RooRealVar bar5_2018("bar5_2018", "bar5_2018", -0.5, -20, 20);
		RooRealVar bf5_2018("bf5_2018","bf5_2018",0.2,0.0,1.0);
		RooExponential exp5_2018("exp5_2018", "exp5_2018", *m2mu, bar5_2018);
                RooRealVar bar6_2018("bar6_2018", "bar6_2018", -0.5, -20, 20);
                RooRealVar bf6_2018("bf6_2018","bf6_2018",0.2,0.0,1.0);
                RooExponential exp6_2018("exp6_2018", "exp6_2018", *m2mu, bar6_2018);
                RooRealVar bar7_2018("bar7_2018", "bar7_2018", -0.5, -20, 20);
                RooExponential exp7_2018("exp7_2018", "exp7_2018", *m2mu, bar7_2018);

		RooArgList explist_2018(exp1_2018,exp2_2018,exp3_2018,exp4_2018,exp5_2018,exp6_2018,exp7_2018);
		RooArgList expclist_2018(bf1_2018,bf2_2018,bf3_2018,bf4_2018,bf5_2018,bf6_2018);
		RooAddPdf bkg_model_exp7_2018("bkg_model_exp7_2018","bkg_model_exp7_2018",explist_2018,expclist_2018,true);
		bkg_model_exp7_2018.fitTo(data_obs);



                RooRealVar pow_1_2017("pow_1_2017","exponent of power law",0,-20,20);
                RooRealVar pf1_2017("pf1_2017","frac of power law",0.2,0.0,1.0);
                RooGenericPdf plaw1_2017("plaw1_2017","TMath::Power(@0,@1)",RooArgList(*m2mu,pow_1_2017));
                RooRealVar qar1_2017("qar1_2017", "qar1_2017", 0.2, 0, 20);
                RooRealVar qar2_2017("qar2_2017", "qar2_2017", 1.5, 0, 20);
                RooRealVar qar3_2017("qar3_2017", "qar3_2017", 2.0, 0, 20);
                RooArgList qlist_2017(qar1_2017, qar2_2017, qar3_2017);
                RooRealVar bfp1_2017("bfp1_2017","frac of bernstein",0.2,0.0,1.0);
                RooBernstein bern4_2017("bkg_model_bern3_2017", "bkg_model_bern3_2017", *m2mu, qlist_2017);
                RooArgList plawlist1_2017(plaw1_2017, bern4_2017);
                RooArgList plawclist1_2017(pf1_2017, bfp1_2017);
                RooAddPdf bkg_model_bern4p1_2017("bkg_model_bern4p1_2017","bkg_model_bern4p1_2017",plawlist1_2017,plawclist1_2017,true);
                bkg_model_bern4p1_2017.fitTo(data_obs);



                RooRealVar pow_1_2018("pow_1_2018","exponent of power law",0,-20,20);
                RooRealVar pf1_2018("pf1_2018","frac of power law",0.2,0.0,1.0);
                RooGenericPdf plaw1_2018("plaw1_2018","TMath::Power(@0,@1)",RooArgList(*m2mu,pow_1_2018));
                RooRealVar qar1_2018("qar1_2018", "qar1_2018", 0.2, 0, 20);
                RooRealVar qar2_2018("qar2_2018", "qar2_2018", 1.5, 0, 20);
                RooRealVar qar3_2018("qar3_2018", "qar3_2018", 2.0, 0, 20);
                RooArgList qlist_2018(qar1_2018, qar2_2018, qar3_2018);
                RooRealVar bfp1_2018("bfp1_2018","frac of bernstein",0.2,0.0,1.0);
                RooBernstein bern4_2018("bkg_model_bern3_2018", "bkg_model_bern3_2018", *m2mu, qlist_2018);
                RooArgList plawlist1_2018(plaw1_2018, bern4_2018);
                RooArgList plawclist1_2018(pf1_2018, bfp1_2018);
                RooAddPdf bkg_model_bern4p1_2018("bkg_model_bern4p1_2018","bkg_model_bern4p1_2018",plawlist1_2018,plawclist1_2018,true);
                bkg_model_bern4p1_2018.fitTo(data_obs);





                RooCategory pdf_index_2017("pdf_index_2017","Index of the background PDF which is active");
		RooArgList bkg_pdf_list_2017;
                bkg_pdf_list_2017.add(bkg_model_bern4_2017);
                bkg_pdf_list_2017.add(bkg_model_pol4xexp_2017);
                bkg_pdf_list_2017.add(bkg_model_exp7_2017);
                bkg_pdf_list_2017.add(bkg_model_bern4p1_2017);
                RooMultiPdf bkg_model_2017("bkg_model_2017", "All Pdfs", pdf_index_2017, bkg_pdf_list_2017);
		bkg_model_2017.setCorrectionFactor(0.0);

                RooCategory  pdf_index_2018("pdf_index_2018","Index of the background PDF which is active");
                RooArgList bkg_pdf_list_2018;		
                bkg_pdf_list_2018.add(bkg_model_bern4_2018);
                bkg_pdf_list_2018.add(bkg_model_pol4xexp_2018);
                bkg_pdf_list_2018.add(bkg_model_exp7_2018);
                bkg_pdf_list_2018.add(bkg_model_bern4p1_2018);
                RooMultiPdf bkg_model_2018("bkg_model_2018", "All Pdfs", pdf_index_2018, bkg_pdf_list_2018);	       
                bkg_model_2018.setCorrectionFactor(0.0);




                gStyle->SetOptTitle(0);
                gStyle->SetOptStat(0);

                gStyle->SetPadTopMargin(0.07);
                gStyle->SetPadBottomMargin(0.3);
                gStyle->SetPadLeftMargin(0.12);
                gStyle->SetPadRightMargin(0.07);

                gStyle->SetNdivisions(508, "X");
                gStyle->SetNdivisions(508, "Y");

                TCanvas c_all("c_all", "c_all", 800, 1000);
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


		//save into ROO workspace
		RooWorkspace dpworkspace("dpworkspace", "");
		dpworkspace.import(data_obs);
		dpworkspace.import(*signalModel);
		if (year[y] == "2017"){
		  dpworkspace.import(bkg_model_2017);
		}else if (year[y] == "2018"){
		  dpworkspace.import(bkg_model_2018); 
		}
		dpworkspace.writeToFile(Form("output/dpWorkspace"+year[y]+suff+"_%d.root",i));

		//write the datacard
		char inputShape[200];
		sprintf(inputShape,"output/dpCard_"+year[y]+suff+"_m%.3f_%d.txt",mass,i);
		ofstream newcardShape;
		newcardShape.open(inputShape);
		newcardShape << Form("imax * number of channels\n");
		newcardShape << Form("jmax * number of background\n");
		newcardShape << Form("kmax * number of nuisance parameters\n");
		newcardShape << Form("shapes data_obs	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:data_obs\n",i);
		newcardShape << Form("shapes bkg_mass	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:bkg_model_"+year[y]+"\n",i);
		newcardShape << Form("shapes signalModel_generic	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:signalModel_generic\n",i);
		newcardShape << Form("bin		CatAB\n");
		newcardShape << Form("observation 	-1.0\n");
		newcardShape << Form("bin     		CatAB		CatAB		\n");
		newcardShape << Form("process 		signalModel_generic  	bkg_mass	\n");
		newcardShape << Form("process 		0    		1	   	\n");
		newcardShape << Form("rate    		%f  		%f		\n",
				     effcuts*effgraph->Eval(mass,0,"S")*luminosity, catA->Integral());
		//newcardShape << Form("lumi13TeV_2017 lnN 	1.023 	-\n");
		newcardShape << Form("lumi13TeV_2018 lnN 	1.026 	-\n");
		newcardShape << Form("id_eff_mva_2018 lnN	%f 	-\n", selSys);
		newcardShape << Form("eff_trig_2018 lnN         %f        -\n", triggSys);
		newcardShape << Form("pdf_index_"+year[y]+" discrete \n");
		//newcardShape << Form("sig_shape_2018 lnN        1.10 	-\n");
		//newcardShape << Form("eff_mu_13TeV_2017 lnN	1.015 	-\n");
		if (year[y] == "2017"){
		  newcardShape << Form("bkg_norm_2017 rateParam CatAB bkg_mass %f\n",1.0);
		}
		if (year[y] == "2018"){
		  newcardShape << Form("bkg_norm_2018 rateParam CatAB bkg_mass %f\n",1.0);
		}
		newcardShape << Form("res_unc param 0.0 1.0 \n");
		newcardShape << Form("");

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



                TH1 *testhisto = data_obs.createHistogram("testhisto", *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                TH1 *testpdf1;
                TH1 *testpdf2;
                TH1 *testpdf3;
                TH1 *testpdf4;
                TH1 *testpdf5;
                TH1 *testpdf6;
                if (year[y] == "2017"){
                  testpdf1   = bkg_model_bern4_2017.createHistogram("testpdf1"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf2   = bkg_model_pol4xexp_2017.createHistogram("testpdf2"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf3   = bkg_model_exp7_2017.createHistogram("testpdf3"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf4   = bkg_model_bern4p1_2017.createHistogram("testpdf4"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  //testpdf5   = bkg_model_bern5_2017.createHistogram("testpdf5"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  //testpdf6   = bkg_model_bern6_2017.createHistogram("testpdf6"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                }else if (year[y] == "2018"){
                  testpdf1   = bkg_model_bern4_2018.createHistogram("testpdf1"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf2   = bkg_model_pol4xexp_2018.createHistogram("testpdf2"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf3   = bkg_model_exp7_2018.createHistogram("testpdf3"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  testpdf4   = bkg_model_bern4p1_2018.createHistogram("testpdf4"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  //testpdf5   = bkg_model_bern5_2018.createHistogram("testpdf5"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                  //testpdf6   = bkg_model_bern6_2018.createHistogram("testpdf6"  , *m2mu, RooFit::Binning(nbins, massLow, massHigh));
                }


                testpdf1->Scale(testhisto->Integral(1,-1)/testpdf1->Integral(1,-1));
                testpdf2->Scale(testhisto->Integral(1,-1)/testpdf2->Integral(1,-1));
                testpdf3->Scale(testhisto->Integral(1,-1)/testpdf3->Integral(1,-1));
                testpdf4->Scale(testhisto->Integral(1,-1)/testpdf4->Integral(1,-1));
                //testpdf5->Scale(testhisto->Integral(1,-1)/testpdf5->Integral(1,-1));
                //testpdf6->Scale(testhisto->Integral(1,-1)/testpdf6->Integral(1,-1));

                int noofparm1 = 4;
                int noofparm2 = 5;
                int noofparm3 = 13;
                int noofparm4 = 6;
                int noofparm5 = 5;
                int noofparm6 = 6;
                Double_t mychi21 = 0.;
                Double_t mychi22 = 0.;
                Double_t mychi23 = 0.;
                Double_t mychi24 = 0.;
                Double_t mychi25 = 0.;
                Double_t mychi26 = 0.;

                for(int i=1; i<=testhisto->GetNbinsX(); i++){
                  Double_t tmp_mychi21 = ((testpdf1->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf1->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf1->GetBinContent(i);
                  Double_t tmp_mychi22 = ((testpdf2->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf2->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf2->GetBinContent(i);
                  Double_t tmp_mychi23 = ((testpdf3->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf3->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf3->GetBinContent(i);
                  Double_t tmp_mychi24 = ((testpdf4->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf4->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf4->GetBinContent(i);
                  //Double_t tmp_mychi25 = ((testpdf5->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf5->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf5->GetBinContent(i);
                  //Double_t tmp_mychi26 = ((testpdf6->GetBinContent(i) - testhisto->GetBinContent(i)) * (testpdf6->GetBinContent(i) - testhisto->GetBinContent(i))) / testpdf6->GetBinContent(i);
                  //cout << "testpdf1 increment " << testpdf1->GetBinContent(i) << "\n";
                  mychi21 += tmp_mychi21;
                  mychi22 += tmp_mychi22;
                  mychi23 += tmp_mychi23;
                  mychi24 += tmp_mychi24;
                  //mychi25 += tmp_mychi25;
                  //mychi26 += tmp_mychi26;

                }
                Double_t mychi2_final1 = mychi21/(nbins - (noofparm1 + 1));
                Double_t mychi2_final2 = mychi22/(nbins - (noofparm2 + 1));
                Double_t mychi2_final3 = mychi23/(nbins - (noofparm3 + 1));
                Double_t mychi2_final4 = mychi24/(nbins - (noofparm4 + 1));
                //Double_t mychi2_final5 = mychi25/(nbins - (noofparm5 + 1));
		//Double_t mychi2_final6 = mychi26/(nbins - (noofparm6 + 1));
		
		//if ( mychi2_final1 < mychi2_final2  &&  mychi2_final1 < mychi2_final3  &&  mychi2_final1 < mychi2_final4 )
		  chi1Values->SetPoint(i, mass, mychi2_final1);
		  //if ( mychi2_final2 < mychi2_final1  &&  mychi2_final2 < mychi2_final3  &&  mychi2_final2 < mychi2_final4 )
		   chi2Values->SetPoint(i, mass, mychi2_final2);
		   //if ( mychi2_final3 < mychi2_final1  &&  mychi2_final3 < mychi2_final2  &&  mychi2_final3 < mychi2_final4 )
		  chi3Values->SetPoint(i, mass, mychi2_final3);
		  //if ( mychi2_final4 < mychi2_final1  &&  mychi2_final4 < mychi2_final2  &&  mychi2_final4 < mychi2_final3 )
		  chi4Values->SetPoint(i, mass, mychi2_final4);
		 /*if ( mychi2_final5 < mychi2_final1  &&  mychi2_final5 < mychi2_final2  &&  mychi2_final5 < mychi2_final3 &&  mychi2_final5 < mychi2_final4 &&  mychi2_final5 < mychi2_final6)
		  chi5Values->SetPoint(i, mass, mychi2_final5);
		 if ( mychi2_final6 < mychi2_final1  &&  mychi2_final6 < mychi2_final2  &&  mychi2_final6 < mychi2_final3 &&  mychi2_final6 < mychi2_final4 &&  mychi2_final6 < mychi2_final5)
		  chi6Values->SetPoint(i, mass, mychi2_final6);
		 */

		  Double_t RooPlot_chi21 = frame->chiSquare("new_pdf", "new_hist", noofparm1 + 1);
		  Double_t RooPlot_chi22 = frame->chiSquare("new_pdf", "new_hist", noofparm2 + 1);
		  Double_t RooPlot_chi23 = frame->chiSquare("new_pdf", "new_hist", noofparm3 + 1);
		  Double_t RooPlot_chi24 = frame->chiSquare("new_pdf", "new_hist", noofparm4 + 1);

		  TLatex latex1;
		  latex1.SetNDC();
		  latex1.SetTextSize(0.04);
		  latex1.SetTextAlign(33);
		  TString h_string1 = "Chi2/ndf = " + std::to_string(mychi2_final1);
		  h_string1 = "#chi^{2}/(N_{bins}-N_{d.o.f}+1) =" + std::to_string(mychi2_final1);
		  latex1.DrawLatex(0.92,0.980, h_string1);

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
		  line1->SetLineColor(kBlack);
		  line1->SetLineWidth(1);
		  line1->Draw("same");

		  TH1F * Ratio1 = ((TH1F*)testhisto->Clone("Ratio1"));
		  for(int i = 1; i<testhisto->GetNbinsX(); i++){
		    Ratio1->SetBinContent(i, (testhisto->GetBinContent(i) - testpdf1->GetBinContent(i))/testhisto->GetBinError(i));
		  }
		  Ratio1->SetFillColor(3);
		  Ratio1->SetLineColor(3);
		  Ratio1->Draw("histsame");

		  TH1F * Ratio2 = ((TH1F*)testhisto->Clone("Ratio2"));
		  for(int i = 1; i<testhisto->GetNbinsX(); i++){
		    Ratio2->SetBinContent(i, (testhisto->GetBinContent(i) - testpdf2->GetBinContent(i))/testhisto->GetBinError(i));
		  }
		  Ratio2->SetFillStyle(3);
		  Ratio2->Draw("histsame");
		  c_all.SaveAs(Form("output/catA_slimWindow_pol3VSpol4_%d_"+year[y]+".png",i));

	}
	f_ws->Close();

	
	TCanvas c_fVal("c_fVal", "c_fVal", 950, 1020);
	chi1Values->GetXaxis()->SetRangeUser(0.8, 9);
	//TCanvas c_fVal("c_fVal", "c_fVal", 950, 820);
	//chi1Values->GetXaxis()->SetRangeUser(1.1, 2.2);

	chi1Values->GetXaxis()->SetTitle("m_{#mu#mu}");
	chi1Values->SetTitle("#chi^{2} / N_{d.o.f.}");
	chi1Values->SetMarkerColor(2);
	chi1Values->SetMarkerStyle(20);
	chi1Values->SetMarkerSize(1);
        chi1Values->Draw("ap");
	chi2Values->SetMarkerColor(3);
	chi2Values->SetMarkerStyle(20);
	chi2Values->SetMarkerSize(1);
	chi2Values->Draw("sp");
	chi3Values->SetMarkerColor(4);
	chi3Values->SetMarkerStyle(20);
	chi3Values->SetMarkerSize(1);
	chi3Values->Draw("sp");
	chi4Values->SetMarkerColor(5);
	chi4Values->SetMarkerStyle(20);
	chi4Values->SetMarkerSize(1);
	chi4Values->Draw("sp");
	/*
	chi5Values->SetMarkerColor(6);
	chi5Values->SetMarkerStyle(20);
	chi5Values->SetMarkerSize(1);
	chi5Values->Draw("sp");
	chi6Values->SetMarkerColor(7);
	chi6Values->SetMarkerStyle(20);
	chi6Values->SetMarkerSize(1);
	chi6Values->Draw("sp");
	*/
	auto legend = new TLegend(0.6,0.3,0.9,0.4);
	legend->AddEntry(chi1Values,"4th order Bern.","p");
	legend->AddEntry(chi2Values,"4th order Cheb.","p");
	legend->AddEntry(chi3Values,"Sum of 7 Exps","p");
	legend->AddEntry(chi4Values,"3rd order Bern + Power law","p");
	//legend->AddEntry(chi5Values,"5th order Bern.","p");
	//legend->AddEntry(chi6Values,"6th order Bern.","p");
        legend->SetBorderSize(0);
	legend->Draw();
	//c_fVal.SetLogy();
        c_fVal.SaveAs("Chi2Dist"+year[y]+".png");
	


  }

}
