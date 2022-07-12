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

void fullDist(TString year="2017"){
  

  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  	TFile* file = NULL;  // Above 3 GeV
  	TFile* file2 = NULL; // Below 3 GeV
  	TFile* file3 = NULL; // Below 3 GeV
  	TFile* file4 = NULL; // Below 3 GeV
	//if (year == "2017") file=TFile::Open("/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2017/ScoutingRunD/mergedHistos_v1.root");
	if (year == "2017"){
	  file=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_mva_2017.root"); //38.7
	  file2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_jpsi0p015_2017.root"); //38.7
               
	  file3=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_mva_2018.root"); //61.3 fb -1
	  file4=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
	}
  //PREPARE EXPECTED NUMBER OF SIGNAL EVENTS PER CATEGORY
	//X-SECTION GRAPH
	//double m[17] 		= {0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 14.0, 16.0, 18.0, 20.0};
	//double xSec[17] 	= {10,  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  10,  10,  10,  10,  10};//[pb]
	double m[11] 		= {2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,   9.0,   10.0,  12.5,  20.0};
	double xSec[11] 	= {8541, 7514, 3323, 2055, 1422, 1043, 793.6, 621.1, 484.3, 292.8, 98.95};//[pb], model-dependent
	//TGraph* xsecgraph 	= new TGraph(17,m,xSec);
	TGraph* xsecgraph 	= new TGraph(11,m,xSec);

	//ACCEPTANCE
	/*TFile* acc_file = TFile::Open("acceptances_dy.root");
	TEfficiency* acc_teff = (TEfficiency*)acc_file->Get("cmsacc");
	int nbins_acc=acc_teff->GetPassedHistogram()->GetNbinsX();
	double acceptances[nbins_acc];
	double m_acceptances[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
		acceptances[j-1] = acc_teff->GetEfficiency(j);
		m_acceptances[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
		}*/
	TFile* acc_file = TFile::Open("acc_dyturbo.root");
	TH1D* acc_teff = (TH1D*)acc_file->Get("s_m");
	int nbins_acc=acc_teff->GetNbinsX();
	double acceptances[nbins_acc];
	double m_acceptances[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
		acceptances[j-1] = acc_teff->GetBinContent(j);
		cout<<"The acceptence is \n"<<acceptances[j-1]<<endl;		
		m_acceptances[j-1] = acc_teff->GetBinCenter(j);
		}
	TGraph* accgraph 	= new TGraph(nbins_acc,m_acceptances,acceptances);
	
	//TF1* accF = (TF1*)acc_file->Get("fit_func");
	//LUMINOSITY
	double luminosity = 0; //4000.;//pb-1
	if (year == "2017") luminosity = 35300;
	//if (year == "2017") luminosity = 4000;

	//else if (year == "2018") luminosity = 6600;
	else if (year == "2018") luminosity = 61300;
	//else if (year == "2018") luminosity = 1183;

	//EFFICIENCY
	//get acceptance from hist
	//TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_"+year+"_mll_dR_wieghted.root");
	TFile* eff_file = TFile::Open("modifiedMllEff"+year+".root");
	//TFile* eff_file = TFile::Open("l1_corrCuts_eff_Data_newAllTrigLowMass_2018_mll_dR_wieghted.root");
	TEfficiency *teff;
	if (year == "2017"){
	  teff = ((TEfficiency*)eff_file->Get("honemllD_clone"));
	}
	else if (year == "2018"){
	  teff = ((TEfficiency*)eff_file->Get("honemll_clone"));
	}

	//cout<<teff<<endl;
	teff->Draw();
	teff->Paint("");
	TGraphAsymmErrors* effgraph = teff->GetPaintedGraph();
	//EFFICIENCY SYSTEMATIC
	//get estimated efficiency from list
	TFile* trigEffSystFile = TFile::Open("SystEst"+year+".root");
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

	        TH1D* catA;
	        TH1D* catB;
	        TH1D* catC;
	        TH1D* catD;
		catA=(TH1D*)file->Get(Form("massforLimitFull"));
		catB=(TH1D*)file2->Get(Form("massforLimitFull"));
		catC=(TH1D*)file3->Get(Form("massforLimitFull"));
		catD=(TH1D*)file4->Get(Form("massforLimitFull"));
	  	//we're using only one category, so we sum all histos
	  	catA->Add(catC);
	  	catB->Add(catD);
	  	delete catC;
	  	delete catD;

		//repeat for both the histograms w/ and w/o the MVA ID
	  	//get the histograms
		


	  	double massLow  =  catA->GetXaxis()->GetXmin();
		double massHigh =  catA->GetXaxis()->GetXmax();
		double massBinWidth = massHigh-massLow;
	  

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


		//cout<<"Spline: "<<effAgraph->Eval(mass,0,"S")<<endl;
		//cout<<"Graph : "<<effAgraph->Eval(mass)<<endl;

		//define the signal model
		// la sua variabile massa osservabile come massa qui, semplicemente cambiandogli il range.

		//RooDataHist data_obs("data_obs", "", RooArgList(*m2mu), catA);


		//RooPlot *frame = m2mu->frame();
		//data_obs.plotOn(frame);
		TCanvas c_all("c_all", "c_all", 1200, 750);
		catA->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
		catA->GetYaxis()->SetTitle("Events / (0.002 GeV)");
		catA->SetTitle("");
		catA->SetStats(1==0);
		catA->Draw("");
		catB->SetLineColor(2);
		catB->Draw("same");

		auto legend = new TLegend(0.55,0.75,0.87,0.85);
		legend->AddEntry(catA,"Upsilon Trained Selection","l");
		legend->AddEntry(catB,"J/#psi Trained Selection","l");
		legend->SetBorderSize(0);
		legend->Draw();
		auto cmsTag = new TLatex(0.13,0.917,"#scale[1.1]{CMS}");
		cmsTag->SetNDC();
		cmsTag->SetTextAlign(11);
		cmsTag->Draw();
		auto cmsTag2 = new TLatex(0.2,0.917,"#scale[0.825]{#bf{#it{Preliminary}}}");
		cmsTag2->SetNDC();
		cmsTag2->SetTextAlign(11);
		cmsTag2->Draw();
		auto tag = new TLatex(0.52,0.917,"#scale[0.8]{Scouting Triggers, 96.6 fb^{-1} (13 TeV)}");
		tag->SetNDC();
		tag->SetTextAlign(11);
		tag->Draw();

		c_all.SetLogy();
		c_all.SaveAs(Form("fullMassDist.pdf"));



}
