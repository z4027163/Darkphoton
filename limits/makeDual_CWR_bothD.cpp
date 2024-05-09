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

void makeDual_CWR_bothD(){

  TString year[2] = {"2017","2018"};
/*
  float norm_dual[2] = {4815.7,14754.5};    //dkk SR
  float norm2_dual[2]={5022.8,15410.5};     //dkpi SR
  float norm03_dual[2]={22182.1,15528.3};   //dkk CR
  float norm3_dual[2]={23191.1,16248.1};    //dkpi CR
*/

  float norm_dual[2] = {5650.4,16883};    //dtokk SR
  float norm2_dual[2]={5103,15246};      //dtokk CR
  float norm03_dual[2]={26879,18694};    //dtokpi SR
  float norm3_dual[2]={29771,20600};     //dtokpi CR


  for(int y = 0; y < 2; y++){ //year
  

  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  	TFile* file = NULL;  // Above 3 GeV
  	TFile* file2 = NULL; // Below 3 GeV
	//if (year == "2017") file=TFile::Open("/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2017/ScoutingRunD/mergedHistos_v1.root");
        if (year[y] == "2017"){
          file=TFile::Open("~/dphist/CWR/pase0p5/mergedHistos_mva_2017.root"); //38.7
          file2=TFile::Open("~/dphist/CWR/pase0p5/mergedHistos_jpsi0p015_2017.root");
        }
        else if (year[y] == "2018"){
          file=TFile::Open("~/dphist/CWR/pase0p5/mergedHistos_mva_2018.root"); //61.3 fb -1
          file2=TFile::Open("~/dphist/CWR/pase0p5/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
        }
  //PREPARE EXPECTED NUMBER OF SIGNAL EVENTS PER CATEGORY
	//X-SECTION GRAPH
	//double m[17] 		= {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 14.0, 16.0, 18.0, 20.0};
	//double xSec[17] 	= {10,  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  10,  10,  10,  10,  10};//[pb]
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
	//if (year[y] == "2017") luminosity = 8380;
        //if (year[y] == "2017") luminosity = 12630;

	//else if (year[y] == "2018") luminosity = 13980;
	else if (year[y] == "2018") luminosity = 61300;
	//else if (year[y] == "2018") luminosity = 21040;

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
        //reverse
//        { int N=effgraph->GetN(); int i=0; while(i<N) {effgraph->SetPointY(i,effgraph->GetPintY(i)*0.05)} }

	//EFFICIENCY SYSTEMATIC
	//get estimated efficiency from list
	TFile* trigEffSystFile = TFile::Open("SystEst"+year[y]+".root");
	TH1F *tsys = ((TH1F*)trigEffSystFile->Get("triggerSys"));
	int nbins_tsys=tsys->GetNbinsX();


	effgraph->SaveAs("output_dual/effgraph.root");
	accgraph->SaveAs("output_dual/accF.root");
	xsecgraph->SaveAs("output_dual/xsecgraph.root");

	
	//scale
	double eps2scale = 1;//0.01;//1;//0.1;//0.002; //this scales eps 0.02 -> 

	//*****----->>>>> nSignal = xsecgraph->Eval(mass)*eff*luminosity*acceptances[i]

	//define unfittable ranges
	//float unfittable_mins[9] = {0,    0.28, 0.5,   0.7,   0.95, 2.75,  3.55,  8.9};

	double unfittable_regions[8][2] = {{0,0.22}, {0.53,0.575}, {0.74,0.85}, {0.97,1.12}, {2.8,3.85}, {9.0,11}};



	//ID EFFICIENCY
	TFile* IDfileMVA = NULL;
	IDfileMVA=TFile::Open("~/dphist/lowDY/lowDY_mva.root");	
	TFile* IDfileMVA2 = NULL;
	IDfileMVA2=TFile::Open("~/dphist/lowDY/lowDY_mvajpsi.root");	
//        IDfileMVA2=TFile::Open("~/dphist/reverse/lowDY/lowDY_mvajpsi.root");
	TFile* IDfileNO = NULL;
	IDfileNO=TFile::Open("~/dphist/lowDY/lowDY_noid.root");	



   //LOOP OVER MASS INDICES AND MAKE THE CARDS/WORKSPACES
        TFile* dkk_ws = TFile::Open("dkk_ws.root", "READ");
        RooWorkspace *w_kk = (RooWorkspace*)dkk_ws->Get("dpworkspace");
        w_kk->var("dkkmean")->setConstant(true);//all this should actually be automatic. Check!!
        w_kk->var("dkkwidth")->setConstant(true);
        w_kk->var("dkksigma")->setConstant(true);

        TFile* dkpi_ws = TFile::Open("dkpi_ws.root", "READ");
        RooWorkspace *w_kpi = (RooWorkspace*)dkpi_ws->Get("dpworkspace");
        w_kpi->var("dkpimean")->setConstant(true);//all this should actually be automatic. Check!!
        w_kpi->var("sigL1")->setConstant(true);
        w_kpi->var("sigR1")->setConstant(true);

        TFile* dpp_ws = TFile::Open("dpp_ws.root", "READ");
        RooWorkspace *w_pp = (RooWorkspace*)dpp_ws->Get("dpworkspace");
        w_pp->var("dppmean")->setConstant(true);//all this should actually be automatic. Check!!
        w_pp->var("dppwidth")->setConstant(true);
        w_pp->var("dppsigma")->setConstant(true);


	double mass = -1.;
	double rel_reso=0.013;//temporary
        double rel_reso_kk = 0.0153;
	TFile* f_ws = TFile::Open(("../mass_calibration/pdfs"+(string)year[y]+".root").c_str(), "READ");
	RooWorkspace *w = (RooWorkspace*)f_ws->Get("dpworkspace");
	w->loadSnapshot("calibrated");
	w->var("alpha1")->setConstant(true);//all this should actually be automatic. Check!!
	w->var("alpha2")->setConstant(true);
	w->var("n1")->setConstant(true);
	w->var("frac_gau")->setConstant(true);
	w->var("gau_reso_scale")->setConstant(true);

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

        double rescale[6]={0.05,0.142,0.771,0.018,0.896,1.0}; //[0] is reverse id
        TString trg[6]={"","_trg1","_trg2","_trg3","_trg4",""};
        int kk=5;
       
	for(int i=416; i<430; i++){
	  	//get the histograms
	        TH1D* catA;
	        TH1D* catB;
	        if (i > 580){ 
		  catA=(TH1D*)file->Get(Form("massforLimit_CatA%d",i));
		  catB=(TH1D*)file->Get(Form("massforLimit_CatB%d",i));
	        }
		else{
		  catA=(TH1D*)file2->Get(Form("massforLimit%s_CatA%d",trg[kk].Data(),i));
		  catB=(TH1D*)file2->Get(Form("massforLimit%s_CatB%d",trg[kk].Data(),i));
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
                int imap=round(0.5*i-0.5);
		if (i > 580){
		  catAMVA=(TH1D*)IDfileMVA->Get(Form("massforLimit_CatA%d",imap));
		  catBMVA=(TH1D*)IDfileMVA->Get(Form("massforLimit_CatB%d",imap));
		} 
		else{
		  catAMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatA%d",imap));
                  catBMVA=(TH1D*)IDfileMVA2->Get(Form("massforLimit_CatB%d",imap));
		}
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
                //reverse
//                if (mass < 3.0) effcuts = 0.03;
		cout << "The ID efficiency is " << effcuts << " at mass " << mass ; 
		cout << ".  The numerator is " << countMVA << " and the denominator is " << countNO << "\n"; 

                effValues->SetPoint(i, mass, effcuts);
                accValues->SetPoint(i, mass, accgraph->Eval(mass,0,""));
                plotValues->SetPoint(i, mass, effcuts*accgraph->Eval(mass,0,"S")*effgraph->Eval(mass,0,"S")*rescale[kk]);

		//Calculate log normal uncertainty for trigger and selection efficiency
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


                RooRealVar lar1_2017("lar1_2017", "lar1_2017", 0.0, -10.0, 10.0);
                RooRealVar lar2_2017("lar2_2017", "lar2_2017", 0.0, -10.0, 10.0);
                RooRealVar lar3_2017("lar3_2017", "lar3_2017", 0.0, -10.0, 10.0);
                RooRealVar lar4_2017("lar4_2017", "lar4_2017", 0.0, -10.0, 10.0);
                RooArgList llist_2017(lar1_2017, lar2_2017, lar3_2017, lar4_2017);
                RooPolynomial bkg_model_line4_2017("bkg_model_line_2017", "bkg_model_line_2017", *m2mu, llist_2017, 1);
                //Exponentials
                RooRealVar car1_2017("car1_2017", "car1_2017", -0.5, -10, 10);
                RooExponential bkg_model_exp4_2017("bkg_model_exp4_2017", "bkg_model_exp4_2017", *m2mu, car1_2017);
                //Product of the two
                RooProdPdf bkg_model_pol4xexp_2017("bkg_model_pol4xexp_2017", "bkg_model_pol4xexp_2017", bkg_model_line4_2017, bkg_model_exp4_2017);
                bkg_model_pol4xexp_2017.chi2FitTo(data_obs);

                RooRealVar lar1_2018("lar1_2018", "lar1_2018", 0.0, -10.0, 10.0);
                RooRealVar lar2_2018("lar2_2018", "lar2_2018", 0.0, -10.0, 10.0);
                RooRealVar lar3_2018("lar3_2018", "lar3_2018", 0.0, -10.0, 10.0);
                RooRealVar lar4_2018("lar4_2018", "lar4_2018", 0.0, -10.0, 10.0);
                RooArgList llist_2018(lar1_2018, lar2_2018, lar3_2018, lar4_2018);
                RooPolynomial bkg_model_line4_2018("bkg_model_line_2018", "bkg_model_line_2018", *m2mu, llist_2018, 1);
                //Exponentials
                RooRealVar car1_2018("car1_2018", "car1_2018", -0.5, -10, 10);
                RooExponential bkg_model_exp4_2018("bkg_model_exp4_2018", "bkg_model_exp4_2018", *m2mu, car1_2018);
                //Product of the two
                RooProdPdf bkg_model_pol4xexp_2018("bkg_model_pol4xexp_2018", "bkg_model_pol4xexp_2018", bkg_model_line4_2018, bkg_model_exp4_2018);
                bkg_model_pol4xexp_2018.chi2FitTo(data_obs);



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
		
		
		RooRealVar bar1_2017("bar1_2017", "bar1_2017", -0.5, -10, 10);                            
		RooRealVar bf1_2017("bf1_2017","bf1_2017",0.2,0.0,1.0);   
		RooExponential exp1_2017("exp1_2017", "exp1_2017", *m2mu, bar1_2017);
		RooRealVar bar2_2017("bar2_2017", "bar2_2017", -0.5, -10, 10);                                                                                          
		RooRealVar bf2_2017("bf2_2017","bf2_2017",0.2,0.0,1.0);                              
		RooExponential exp2_2017("exp2_2017", "exp2_2017", *m2mu, bar2_2017);
		RooRealVar bar3_2017("bar3_2017", "bar3_2017", -0.5, -10, 10);
		RooRealVar bf3_2017("bf3_2017","bf3_2017",0.2,0.0,1.0);
		RooExponential exp3_2017("exp3_2017", "exp3_2017", *m2mu, bar3_2017);
		RooRealVar bar4_2017("bar4_2017", "bar4_2017", -0.5, -10, 10);
		RooRealVar bf4_2017("bf4_2017","bf4_2017",0.2,0.0,1.0);
		RooExponential exp4_2017("exp4_2017", "exp4_2017", *m2mu, bar4_2017);
		RooRealVar bar5_2017("bar5_2017", "bar5_2017", -0.5, -10, 10);
		RooRealVar bf5_2017("bf5_2017","bf5_2017",0.2,0.0,1.0);
		RooExponential exp5_2017("exp5_2017", "exp5_2017", *m2mu, bar5_2017);
                RooRealVar bar6_2017("bar6_2017", "bar6_2017", -0.5, -10, 10);
                RooExponential exp6_2017("exp6_2017", "exp6_2017", *m2mu, bar6_2017);

		RooArgList explist_2017(exp1_2017,exp2_2017,exp3_2017,exp4_2017,exp5_2017,exp6_2017);
		RooArgList expclist_2017(bf1_2017,bf2_2017,bf3_2017,bf4_2017,bf5_2017);
		RooAddPdf bkg_model_exp7_2017("bkg_model_exp7_2017","bkg_model_exp7_2017",explist_2017,expclist_2017,true);
		bkg_model_exp7_2017.fitTo(data_obs);


		RooRealVar bar1_2018("bar1_2018", "bar1_2018", -0.5, -10, 10);                            
		RooRealVar bf1_2018("bf1_2018","bf1_2018",0.2,0.0,1.0);   
		RooExponential exp1_2018("exp1_2018", "exp1_2018", *m2mu, bar1_2018);
		RooRealVar bar2_2018("bar2_2018", "bar2_2018", -0.5, -10, 10);                                                                                           
		RooRealVar bf2_2018("bf2_2018","bf2_2018",0.2,0.0,1.0);                                                                                                  
		RooExponential exp2_2018("exp2_2018", "exp2_2018", *m2mu, bar2_2018);
		RooRealVar bar3_2018("bar3_2018", "bar3_2018", -0.5, -10, 10);
		RooRealVar bf3_2018("bf3_2018","bf3_2018",0.2,0.0,1.0);
		RooExponential exp3_2018("exp3_2018", "exp3_2018", *m2mu, bar3_2018);
		RooRealVar bar4_2018("bar4_2018", "bar4_2018", -0.5, -10, 10);
		RooRealVar bf4_2018("bf4_2018","bf4_2018",0.2,0.0,1.0);
		RooExponential exp4_2018("exp4_2018", "exp4_2018", *m2mu, bar4_2018);
		RooRealVar bar5_2018("bar5_2018", "bar5_2018", -0.5, -10, 10);
		RooRealVar bf5_2018("bf5_2018","bf5_2018",0.2,0.0,1.0);
		RooExponential exp5_2018("exp5_2018", "exp5_2018", *m2mu, bar5_2018);
                RooRealVar bar6_2018("bar6_2018", "bar6_2018", -0.5, -10, 10);
		RooRealVar bf6_2018("bf6_2018","bf6_2018",0.2,0.0,1.0);
                RooExponential exp6_2018("exp6_2018", "exp6_2018", *m2mu, bar6_2018);
                RooRealVar bar7_2018("bar7_2018", "bar7_2018", -0.5, -10, 10);
                RooExponential exp7_2018("exp7_2018", "exp7_2018", *m2mu, bar7_2018);

		RooArgList explist_2018(exp1_2018,exp2_2018,exp3_2018,exp4_2018,exp5_2018,exp6_2018);
		RooArgList expclist_2018(bf1_2018,bf2_2018,bf3_2018,bf4_2018,bf5_2018);
		RooAddPdf bkg_model_exp7_2018("bkg_model_exp7_2018","bkg_model_exp7_2018",explist_2018,expclist_2018,true);
		bkg_model_exp7_2018.fitTo(data_obs);

                RooRealVar pow_1_2017("pow_1_2017","exponent of power law",0,-10,10);
                RooRealVar pf1_2017("pf1_2017","frac of power law",0.2,0.0,1.0);
                RooGenericPdf plaw1_2017("plaw1_2017","TMath::Power(@0,@1)",RooArgList(*m2mu,pow_1_2017));
                RooRealVar qar1_2017("qar1_2017", "qar1_2017", 0.2, 0, 10);
                RooRealVar qar2_2017("qar2_2017", "qar2_2017", 1.5, 0, 10);
                RooRealVar qar3_2017("qar3_2017", "qar3_2017", 2.0, 0, 10);
                RooArgList qlist_2017(qar1_2017, qar2_2017, qar3_2017);
                RooRealVar bfp1_2017("bfp1_2017","frac of bernstein",0.2,0.0,1.0);
                RooBernstein bern3_2017("bkg_model_bern3_2017", "bkg_model_bern3_2017", *m2mu, qlist_2017);
                RooArgList plawlist1_2017(plaw1_2017, bern3_2017);
                RooArgList plawclist1_2017(pf1_2017, bfp1_2017);
                RooAddPdf bkg_model_bern3p1_2017("bkg_model_bern3p1_2017","bkg_model_bern3p1_2017",plawlist1_2017,plawclist1_2017,true);
                bkg_model_bern3p1_2017.fitTo(data_obs);

                RooRealVar pow_1_2018("pow_1_2018","exponent of power law",0,-10,10);
                RooRealVar pf1_2018("pf1_2018","frac of power law",0.2,0.0,1.0);
                RooGenericPdf plaw1_2018("plaw1_2018","TMath::Power(@0,@1)",RooArgList(*m2mu,pow_1_2018));
                RooRealVar qar1_2018("qar1_2018", "qar1_2018", 0.2, 0, 10);
                RooRealVar qar2_2018("qar2_2018", "qar2_2018", 1.5, 0, 10);
                RooRealVar qar3_2018("qar3_2018", "qar3_2018", 2.0, 0, 10);
                RooArgList qlist_2018(qar1_2018, qar2_2018, qar3_2018);
                RooRealVar bfp1_2018("bfp1_2018","frac of bernstein",0.2,0.0,1.0);
                RooBernstein bern3_2018("bkg_model_bern3_2018", "bkg_model_bern3_2018", *m2mu, qlist_2018);
                RooArgList plawlist1_2018(plaw1_2018, bern3_2018);
                RooArgList plawclist1_2018(pf1_2018, bfp1_2018);
                RooAddPdf bkg_model_bern3p1_2018("bkg_model_bern3p1_2018","bkg_model_bern3p1_2018",plawlist1_2018,plawclist1_2018,true);
                bkg_model_bern3p1_2018.fitTo(data_obs);

                RooCategory pdf_index_2017("pdf_index_2017","Index of the background PDF which is active");
		RooArgList bkg_pdf_list_2017;
                bkg_pdf_list_2017.add(bkg_model_bern4_2017);

                bkg_pdf_list_2017.add(bkg_model_pol4xexp_2017);
                bkg_pdf_list_2017.add(bkg_model_exp7_2017);
                bkg_pdf_list_2017.add(bkg_model_bern3p1_2017);

                RooMultiPdf bkg_model_2017("bkg_model_2017", "All Pdfs", pdf_index_2017, bkg_pdf_list_2017);
		bkg_model_2017.setCorrectionFactor(0.5);

                RooCategory  pdf_index_2018("pdf_index_2018","Index of the background PDF which is active");
                RooArgList bkg_pdf_list_2018;		
                bkg_pdf_list_2018.add(bkg_model_bern4_2018);

                bkg_pdf_list_2018.add(bkg_model_pol4xexp_2018);
                bkg_pdf_list_2018.add(bkg_model_exp7_2018);
                bkg_pdf_list_2018.add(bkg_model_bern3p1_2018);

                RooMultiPdf bkg_model_2018("bkg_model_2018", "All Pdfs", pdf_index_2018, bkg_pdf_list_2018);	       
                bkg_model_2018.setCorrectionFactor(0.5);

                //save dtokk peak
                RooVoigtian* dkkModel = (RooVoigtian*)w_kk->pdf("voig");
                 
		//save into ROO workspace
		RooWorkspace dpworkspace("dpworkspace", "");
		dpworkspace.import(data_obs);
		dpworkspace.import(*signalModel);
                dpworkspace.import(*dkkModel);


                RooVoigtian* dkpiModel = (RooVoigtian*)w_kpi->pdf("Gb1");
                dpworkspace.import(*dkpiModel);

                RooVoigtian* dppModel = (RooVoigtian*)w_pp->pdf("voig_pipi");
                dpworkspace.import(*dppModel);

		if (year[y] == "2017"){
		  dpworkspace.import(bkg_model_2017);
		}else if (year[y] == "2018"){
		  dpworkspace.import(bkg_model_2018); 
		}
		dpworkspace.writeToFile(Form("output_dual/dpWorkspace"+year[y]+suff+"_%d.root",i));

		//write the datacard
		char inputShape[200];
		sprintf(inputShape,"output_dual/dpCard_"+year[y]+suff+"_m%.3f_%d.txt",mass,i);
		ofstream newcardShape;
		newcardShape.open(inputShape);
		newcardShape << Form("imax * number of channels\n");
		newcardShape << Form("jmax * number of background\n");
		newcardShape << Form("kmax * number of nuisance parameters\n");
		newcardShape << Form("shapes data_obs	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:data_obs\n",i);
		newcardShape << Form("shapes bkg_mass	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:bkg_model_"+year[y]+"\n",i);
		newcardShape << Form("shapes signalModel_generic	CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:signalModel_generic\n",i);
                newcardShape << Form("shapes dkk_mass CatAB dpWorkspace"+year[y]+suff+"_%d.root dpworkspace:voig\n",i);
                newcardShape << Form("shapes dkpi_mass CatAB dpWorkspace"+year[y]+"IterV3_%d.root dpworkspace:Gb1\n",i);

		newcardShape << Form("bin		CatAB    \n");
		newcardShape << Form("observation 	-1.0      \n");
		newcardShape << Form("bin     		CatAB		CatAB		CatAB	   CatAB  \n");
		newcardShape << Form("process 		signalModel_generic     bkg_mass      dkk_mass      dkpi_mass  \n");
		newcardShape << Form("process 		0    		1           2           3    \n");
		newcardShape << Form("rate    		%f  		%f	  %f      %f  \n",
				     effcuts*effgraph->Eval(mass,0,"S")*luminosity*rescale[kk], catA->Integral(),norm_dual[y],norm03_dual[y]);
                if (year[y] == "2017"){
                  newcardShape << Form("lumi13TeV_2017 lnN      1.020   -  -  -\n");
                  newcardShape << Form("lumi13TeV_corr lnN      1.011   -  -  -\n");
                  newcardShape << Form("sig_shape_2017 lnN      1.013  -  -  -\n");
                }
                if (year[y] == "2018"){
                  newcardShape << Form("lumi13TeV_2018 lnN      1.015   -  -  -\n");
                  newcardShape << Form("lumi13TeV_corr lnN      1.020   -  -  -\n");
                  newcardShape << Form("sig_shape_2018 lnN      1.027  -  -  -\n");
                }
		newcardShape << Form("id_eff_mva_"+year[y]+" lnN	%f   -    -    -  \n", selSys);
		newcardShape << Form("eff_trig_"+year[y]+" lnN         %f   -    -   -  \n", triggSys);
                double eff_cut_unc=1.08; // ID eff entanglement
                newcardShape << Form("eff_cut lnN         %f      -    -    -\n",eff_cut_unc ); 
                newcardShape << Form("pdf_index_"+year[y]+" discrete \n");
                newcardShape << Form("dkk_ratio lnN            -    -   1.20  - \n");
                newcardShape << Form("dkpi_ratio lnN            -    -   -    1.25  \n");
                newcardShape << Form("sig_ratio lnN            -    -    -     -  \n"); 
		newcardShape << Form("bkg_norm_"+year[y]+" rateParam CatAB bkg_mass 1.0\n");
                newcardShape << Form("res_rel_generic param 0.013 0.0026\n");

		newcardShape.close();

	}
	f_ws->Close();


	for (int j=1; j<=nbins_tsys; j++){
	  double val = 1.00 + abs(tsys->GetBinContent(j));
	  cout <<  "Tris sys " << val << "\n"; 
	}

  }
	
}
