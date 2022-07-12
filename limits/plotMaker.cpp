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

void plotMaker(TString year="2017"){
  

  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  	TFile* file = NULL;  // Above 3 GeV
  	TFile* file2 = NULL; // Below 3 GeV
	//if (year == "2017") file=TFile::Open("/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2017/ScoutingRunD/mergedHistos_v1.root");
	if (year == "2017"){
	  file=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_mva_2017.root"); //38.7
	  file2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_jpsi0p015_2017.root"); //38.7
        }
        else if (year == "2018"){
	  file=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_mva_2018.root"); //61.3 fb -1
	  file2=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
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


	double massGG[8]          = {1., 2., 4., 5., 6., 7., 8., 20.};
	double acc_amumuGG[8]     = {0.0774, 0.0798, 0.0847, 0.1017, 0.1061, 0.1132, 0.1355, 0.3518};
	TGraph* accgraphGG        = new TGraph(8,massGG,acc_amumuGG);
	accgraphGG->SaveAs("testGG.root");

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
	TGraph* effValues = new TGraph(190);
	TGraph* accValues = new TGraph(190);
	TGraph* plotValues = new TGraph(190);

	double rel_reso=0.013;//temporary

	for(int i=0; i<400; i++){
	  	//get the histograms
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
                mass = 0.5*(massLow+massHigh);
                if (mass < 1.0) continue;	  

		//compute mass point and define ROOFit variables
	  	bool dontfit=false;


		double effcuts = countMVA / countNO;
		//if (mass < 2.0) effcuts = 0.383615;
		if (mass < 2.0) effcuts = 0.05*mass+0.68;
		//cout << "The ID efficiency is " << effcuts << " at mass " << mass ; 
		//cout << ".  The numerator is " << countMVA << " and the denominator is " << countNO << "\n"; 


                effValues->SetPoint(i, mass, effcuts);
                accValues->SetPoint(i, mass, accgraph->Eval(mass,0,"S"));
                plotValues->SetPoint(i, mass, effcuts*effgraph->Eval(mass,0,"S"));
		//Calculate log normal uncertainty for trigger efficiency
		double triggSysVal = tsys->GetBinContent(tsys->FindBin(mass));
		double triggSys = 1.00 + abs(triggSysVal); 
		//cout << ".  The value of trigger syst is " <<  triggSys << "\n"; 
		
	}
	double product[nbins_acc];
	double m_product[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
	  double massVal = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	  product[j-1] = acc_teff->GetEfficiency(j)*plotValues->Eval(massVal,0);
	  cout<<"The product is \n"<< product[j-1]<<endl;		
	  m_product[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	} TGraph* prodgraph 	= new TGraph(nbins_acc,m_product,product);

	double acceptenceGG[nbins_acc];
	double m_acceptenceGG[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
	  double massVal = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	  acceptenceGG[j-1] = accgraphGG->Eval(massVal,0);
	  cout<<"The product is \n"<< acceptenceGG[j-1]<<endl;		
	  m_acceptenceGG[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	} TGraph* accgraphGGPlot = new TGraph(nbins_acc,m_acceptenceGG,acceptenceGG);

	double productGG[nbins_acc];
	double m_productGG[nbins_acc];
	for (int j=1; j<=nbins_acc; j++){
	  double massVal = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	  productGG[j-1] = accgraphGG->Eval(massVal,0)*plotValues->Eval(massVal,0);
	  cout<<"The product is \n"<< productGG[j-1]<<endl;		
	  m_productGG[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
	} TGraph* productGGPlot = new TGraph(nbins_acc,m_productGG,productGG);


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

	TCanvas c_fVal("c_fVal", "c_fVal", 1200, 980);
        accgraph->GetYaxis()->SetRangeUser(0.00001, 1);
        accgraph->GetXaxis()->SetRangeUser(0.8, 9);
        accgraph->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
        accgraph->GetYaxis()->SetTitle("Arbitrary Units");
	accgraph->SetMarkerColor(4);
	accgraph->SetMarkerStyle(20);
	accgraph->SetMarkerSize(2);
	c_fVal.cd(4);
	accgraph->SetLineColor(4);
	accgraph->SetLineWidth(2);
        accgraph->Draw("apl");
	accgraph->SetTitle("");

	prodgraph->SetMarkerStyle(21);
	prodgraph->SetMarkerColor(4);
	prodgraph->SetMarkerSize(2);
        prodgraph->SetLineWidth(2);
        prodgraph->SetLineColor(4);
	prodgraph->Draw("SPl");

	accgraphGGPlot->SetMarkerStyle(20);
	accgraphGGPlot->SetMarkerColor(2);
	accgraphGGPlot->SetMarkerSize(2);
	accgraphGGPlot->SetLineColor(2);
        accgraphGGPlot->SetLineWidth(2);
	accgraphGGPlot->Draw("SPl");

	productGGPlot->SetMarkerStyle(21);
	productGGPlot->SetMarkerColor(2);
	productGGPlot->SetMarkerSize(2);
	productGGPlot->SetLineColor(2);
	productGGPlot->SetLineWidth(2);
	productGGPlot->Draw("SPl");

        TLine *l = new TLine(3.6,0,3.6,1);
        //l->SetLineColor(0);
        l->Draw("s");

        auto legend = new TLegend(0.5,0.2,0.87,0.4);
        legend->AddEntry(accgraph,"Acceptance DY","p");
        legend->AddEntry(prodgraph,"Efficiency (ID + Trigger) #times Acceptance DY","p");
        legend->AddEntry(accgraphGGPlot,"Acceptance Gluon Fusion","p");
        legend->AddEntry(productGGPlot,"Efficiency (ID + Trigger) #times Acceptance ggF","p");
	legend->SetBorderSize(0);
        legend->Draw();
	auto cmsTag = new TLatex(0.13,0.917,"#scale[1.1]{CMS}");
	cmsTag->SetNDC();
	cmsTag->SetTextAlign(11);
	cmsTag->Draw();
	auto JpsiTag = new TLatex(0.16,0.817,"#scale[0.5]{J/#psi Trained MVA}");
	JpsiTag->SetNDC();
	JpsiTag->SetTextAlign(11);
	JpsiTag->Draw();
	auto upTag = new TLatex(0.4,0.817,"#scale[0.5]{#Upsilon Trained MVA}");
	upTag->SetNDC();
	upTag->SetTextAlign(11);
	upTag->Draw();
	auto cmsTag2 = new TLatex(0.22,0.917,"#scale[0.825]{#bf{#it{Preliminary}}}");
	cmsTag2->SetNDC();
	cmsTag2->SetTextAlign(11);
	cmsTag2->Draw();
        c_fVal.SetLogy();
        c_fVal.SaveAs("Multiplot.pdf");


	for (int j=1; j<=nbins_tsys; j++){
	  double val = 1.00 + abs(tsys->GetBinContent(j));
	  cout <<  "Tris sys " << val << "\n"; 
	}



}
