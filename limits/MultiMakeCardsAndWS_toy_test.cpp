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
#include "../../../HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

using namespace std;

void MultiMakeCardsAndWS_toy_test(){

  TString year[2] = {"2017","2018"};
  for(int y = 1; y < 2; y++){ //year
  

  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
        TFile* file = NULL;  // Above 3 GeV
        if (year[y] == "2017"){
          file=TFile::Open("/home/tier3/wangzqe/research/darkphoton/CWR_toy_study/output_toyhist_2017.root");
        }
        else if (year[y] == "2018"){
          file=TFile::Open("/home/tier3/wangzqe/research/darkphoton/CWR_toy_study/output_toyhist_2018.root");
        }
	//X-SECTION GRAPH
	double m[11] 		= {2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   9.0,   10.0,  12.5,  20.0};
	double xSec[11] 	= {8541, 7514, 3323, 2055, 1422, 1043, 793.6, 621.1, 484.3, 292.8, 98.95};//[pb], model-dependent
	//TGraph* xsecgraph 	= new TGraph(17,m,xSec);
	TGraph* xsecgraph 	= new TGraph(11,m,xSec);

	//SELECTION UNCERTAINTY
	TFile* uncFile = TFile::Open("~/dphist/unc/idunc_mass.root");
	TH1F* selUncH      =     NULL;
	if (year[y] == "2017"){
	  selUncH = (TH1F*)uncFile->Get("unc2017");
	} else if (year[y] == "2018"){
	  selUncH = (TH1F*)uncFile->Get("unc2017");
	}
	
	//LUMINOSITY
	double luminosity = 0; //4000.;//pb-1
	if (year[y] == "2017") luminosity = 35300;

	else if (year[y] == "2018") luminosity = 61300;

	//EFFICIENCY
	//get acceptance from hist
	//TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_"+year[y]+"_mll_dR_wieghted.root");
	TFile* eff_file = TFile::Open("llptCutEfficiencies/modifiedMllEffD"+year[y]+"high.root");
        TFile* eff_file1 = TFile::Open("llptCutEfficiencies/modifiedMllEffD"+year[y]+"low.root");
	TEfficiency *teff, *teff1;
	teff = ((TEfficiency*)eff_file->Get("honemllD_clone"));
	teff1 = ((TEfficiency*)eff_file1->Get("honemllD_clone"));
        teff->Draw();
        teff->Paint("");
        teff1->Draw();
        teff1->Paint("");

	TGraphAsymmErrors* effgraph = teff->GetPaintedGraph();
        TGraphAsymmErrors* effgraph1 = teff1->GetPaintedGraph();

	//EFFICIENCY SYSTEMATIC
	//get estimated efficiency from list
	TFile* trigEffSystFile = TFile::Open("SystEst"+year[y]+".root");
	TH1F *tsys = ((TH1F*)trigEffSystFile->Get("triggerSys"));
	int nbins_tsys=tsys->GetNbinsX();

	//scale
	double eps2scale = 1;//0.01;//1;//0.1;//0.002; //this scales eps 0.02 -> 

	double unfittable_regions[8][2] = {{0,0.22}, {0.53,0.575}, {0.74,0.85}, {0.97,1.12}, {2.8,3.85}, {9.0,11}};

	//ID EFFICIENCY
	TFile* IDfileMVA = NULL;
	IDfileMVA=TFile::Open("~/dphist/lowDY/lowDY_mva.root");	
	TFile* IDfileMVA2 = NULL;
	IDfileMVA2=TFile::Open("~/dphist/lowDY/lowDY_mvajpsi.root");	
	TFile* IDfileNO = NULL;
	IDfileNO=TFile::Open("~/dphist/lowDY/lowDY_noid.root");	


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

	w->Print();

        // 58,352 pase 0.65%
        // 115,708 pase 0.325%
	for(int i=355; i<738; i++){
                double pvd_scale=1;
                if(i%2==0) continue;
                //if(i!=225) continue;
	  	//get the histograms
	        TH1D* catA;
	        TH1D* catB;
                catA=(TH1D*)file->Get(Form("massforLimit_CatA%d",i));

		//repeat for both the histograms w/ and w/o the MVA ID
	  	//get the histograms
		TH1D* catAMVA;
                TH1D* catBMVA;
                //int imap=round(139.3214+i*0.6511); //pase 0.65%
                int imap=round(i*0.5-0.5);   //pase 0.325%
		  catAMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatA%d",imap));
                  catBMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatB%d",imap));
	  	catAMVA->Add(catBMVA);
		catAMVA->Rebin(2);
	  	delete catBMVA;
		double countMVA = catAMVA->Integral();


	  	TH1D* catANO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatA%d",imap));
	  	TH1D* catBNO=(TH1D*)IDfileNO->Get(Form("massforLimit_CatB%d",imap));
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

	  	mass = 0.5*(massLow+massHigh);
		if (mass < 1.171) continue;
                if(mass>7.5) continue;
		if ((mass > 2.59) && (mass <= 4.25)) continue;

		double effcuts = countMVA / countNO;
		if (mass < 2.0) effcuts = 0.05*mass+0.68;
                if(mass>4.0) effcuts = effcuts*0.84;  //high mass ptcut scale
                if(mass<3.0) effcuts = effcuts*(1.21-0.08*mass);  //low mass ptcut scale

		cout << "The ID efficiency is " << effcuts << " at mass " << mass ; 
		cout << ".  The numerator is " << countMVA << " and the denominator is " << countNO << "\n"; 
                double triggereff=effgraph->Eval(mass,0,"S");
                if(i<590) triggereff=effgraph1->Eval(mass,0,"S");

		//Calculate log normal uncertainty for trigger and selection efficiency
		double triggSysVal = tsys->GetBinContent(tsys->FindBin(mass));
		double triggSys = 1.00 + abs(triggSysVal); 
		cout << ".  The value of trigger syst is " <<  triggSys << "\n"; 

		double selSysVal = selUncH->GetBinContent(selUncH->FindBin(mass));
                // mass   1   2   5   8 GeV
                // 2017  6.5 6.5 5.6 4.4 %
                // 2018  4.1 3.9 4.4 3.9 %
               
               if (year[y] == "2017"){
                   if(i<590) selSysVal=0.065;
                   else selSysVal=0.06;
               }
              if (year[y] == "2018"){
                   selSysVal=0.05;
               }
               //ID unc for jpsi at high mass 
               double selSys = 1.00 + abs(selSysVal);
                cout << ".  The value of sel syst is " <<  selSys << "\n";

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
		RooDataHist data_obs("data_obs", "", RooArgList(*m2mu), catA);
		RooRealVar bkg_norm("bkg_norm", "",catA->Integral());

                RooHistPdf bkg_pdf("bkg_pdf", "Smooth PDF 2", RooArgList(*m2mu), data_obs, 3);                
 
		RooRealVar par1_2017("par1_2017", "par1_2017", 0.2, 0, 10);
		RooRealVar par2_2017("par2_2017", "par2_2017", 1.5, 0, 10);
		RooRealVar par3_2017("par3_2017", "par3_2017", 2.0, 0, 10);
		RooRealVar par4_2017("par4_2017", "par4_2017", 2.0, 0, 10);
		RooArgList alist_2017(par1_2017, par2_2017, par3_2017, par4_2017);
		RooBernstein bkg_model_bern4_2017("bkg_model_bern4_2017", "bkg_model_bern4_2017", *m2mu, alist_2017);
		bkg_model_bern4_2017.fitTo(data_obs);		
		

		RooRealVar par1_2018("par1_2018", "par1_2018", 0.2, 0, 10);
		RooRealVar par2_2018("par2_2018", "par2_2018", 1.5, 0, 10);
		RooRealVar par3_2018("par3_2018", "par3_2018", 2.0, 0, 10);
		RooRealVar par4_2018("par4_2018", "par4_2018", 2.0, 0, 10);
		RooArgList alist_2018(par1_2018, par2_2018, par3_2018, par4_2018);
		RooBernstein bkg_model_bern4_2018("bkg_model_bern4_2018", "bkg_model_bern4_2018", *m2mu, alist_2018);
		bkg_model_bern4_2018.fitTo(data_obs);
	
                	
		//save into ROO workspace
		RooWorkspace dpworkspace("dpworkspace", "");
		dpworkspace.import(data_obs);
		dpworkspace.import(*signalModel);
                dpworkspace.import(bkg_pdf);
		if (year[y] == "2017"){
		  dpworkspace.import(bkg_model_bern4_2017);
		}else if (year[y] == "2018"){
		  dpworkspace.import(bkg_model_bern4_2018); 
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
		newcardShape << Form("shapes bkg_mass	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:bkg_model_bern4_"+year[y]+"\n",i);
		newcardShape << Form("shapes signalModel_generic	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:signalModel_generic\n",i);
		newcardShape << Form("bin		CatAB\n");
		newcardShape << Form("observation 	-1.0\n");
		newcardShape << Form("bin     		CatAB		CatAB		\n");
		newcardShape << Form("process 		signalModel_generic  	bkg_mass	\n");
		newcardShape << Form("process 		0    		1	   	\n");
		newcardShape << Form("rate    		%f  		%f		\n",
				     effcuts*triggereff*luminosity*pvd_scale, catA->Integral());
                if (year[y] == "2017"){
                  newcardShape << Form("lumi13TeV_2017 lnN      1.020   -\n");
                  newcardShape << Form("lumi13TeV_corr lnN      1.011   -\n");
                }
                if (year[y] == "2018"){
                  newcardShape << Form("lumi13TeV_2018 lnN      1.015   -\n");
                  newcardShape << Form("lumi13TeV_corr lnN      1.020   -\n");
                }
		newcardShape << Form("id_eff_mva_2018 lnN	%f 	-\n", selSys);
		newcardShape << Form("eff_trig_2018 lnN         %f        -\n", triggSys);
                newcardShape << Form("pdf_index_"+year[y]+" discrete \n");
		//newcardShape << Form("sig_shape_2018 lnN        1.10 	-\n");
		//newcardShape << Form("eff_mu_13TeV_2017 lnN	1.015 	-\n");
                double eff_cut_unc=1.05; // ID eff entanglement
                eff_cut_unc=1.08;
                newcardShape << Form("eff_cut lnN         %f        -\n",eff_cut_unc );

		if (year[y] == "2017"){
		  newcardShape << Form("bkg_norm_2017 rateParam CatAB bkg_mass 1.0\n");
		}
		if (year[y] == "2018"){
		  newcardShape << Form("bkg_norm_2018 rateParam CatAB bkg_mass 1.0\n");
		}
                newcardShape << Form("res_rel_generic param 0.013 0.0026\n");
		newcardShape << Form("");

		newcardShape.close();

	}
	f_ws->Close();

  }
	
}
