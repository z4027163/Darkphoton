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

void YearRatioTest(){


  //WHICH YEAR
	TString suff="IterV3";
  //INPUT FILE WITH HISTOGRAMS TO FIT BACKGROUND
  	TFile* file_2017 = NULL;  // Above 3 GeV
  	TFile* file2_2017 = NULL; // Below 3 GeV
  	TFile* file_2018 = NULL;  // Above 3 GeV
  	TFile* file2_2018 = NULL; // Below 3 GeV
	//if (year == "2017") file=TFile::Open("/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2017/ScoutingRunD/mergedHistos_v1.root");
	file_2017=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2017.root"); //38.7
	file2_2017=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2017.root"); //38.7
        file_2018=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_mva_2018.root"); //61.3 fb -1
	file2_2018=TFile::Open("/afs/cern.ch/work/w/wangz/public/darkphoton/sigma_p013/mergedHistos_jpsi0p015_2018.root"); //61.3 fb -1
	




	double unfittable_regions[8][2] = {{0,0.22}, {0.53,0.575}, {0.74,0.85}, {0.97,1.12}, {2.8,3.85}, {9.0,11}};



	for(int i=0; i<400; i++){
	  	//get the histograms
                if (i<177) continue;
		if (i>368) continue;
	        TH1D* catA_2017;
	        TH1D* catB_2017;
	        TH1D* catA_2018;
	        TH1D* catB_2018;
	        if (i > 290){ 
		  catA_2017=(TH1D*)file_2017->Get(Form("massforLimit_CatA%d",i));
		  catB_2017=(TH1D*)file_2017->Get(Form("massforLimit_CatB%d",i));
		  catA_2018=(TH1D*)file_2018->Get(Form("massforLimit_CatA%d",i));
		  catB_2018=(TH1D*)file_2018->Get(Form("massforLimit_CatB%d",i));
	        }
		else{
		  catA_2017=(TH1D*)file2_2017->Get(Form("massforLimit_CatA%d",i));
		  catB_2017=(TH1D*)file2_2017->Get(Form("massforLimit_CatB%d",i));
		  catA_2018=(TH1D*)file2_2018->Get(Form("massforLimit_CatA%d",i));
		  catB_2018=(TH1D*)file2_2018->Get(Form("massforLimit_CatB%d",i));
	        }
	  	//TH1D* catC=(TH1D*)file->Get(Form("massforLimit_CatC%d",i));

	  	//we're using only one category, so we sum all histos
	  	catA_2017->Add(catB_2017);
	  	catA_2018->Add(catB_2018);
	  	//catA->Add(catC);
	  	delete catB_2017;
	  	delete catB_2018;


		Double_t factor = 1.;
		catA_2017->Scale(factor/catA_2017->Integral());
		catA_2018->Scale(factor/catA_2018->Integral());

		//Reduce to the difference
		catA_2017->Add(catA_2018, -1);

                gStyle->SetOptTitle(0);
                gStyle->SetOptStat(0);

                gStyle->SetPadTopMargin(0.07);
                gStyle->SetPadBottomMargin(0.3);
                gStyle->SetPadLeftMargin(0.12);
                gStyle->SetPadRightMargin(0.07);

                gStyle->SetNdivisions(508, "X");
                gStyle->SetNdivisions(508, "Y");

                TCanvas c_all("c_all", "c_all", 950, 1020);

                catA_2017->SetFillStyle(3);
                catA_2017->SetLineColor(2);
                catA_2017->Draw("hist");

                c_all.SaveAs(Form("yearCheck/yearDifference%d.png",i));

	}


	

	
}
