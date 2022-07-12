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

void discreteProfileCards(TString year="2017"){
  

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
	TFile* f_ws = TFile::Open(("../mass_calibration/pdfs"+(string)year+".root").c_str(), "READ");
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
                plotValues->SetPoint(i, mass, effcuts*accgraph->Eval(mass,0,"S")*effgraph->Eval(mass,0,"S"));

		//Calculate log normal uncertainty for trigger efficiency
		double triggSysVal = tsys->GetBinContent(tsys->FindBin(mass));
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

	       
		RooRealVar par1("par1", "par1", 0.2, 0, 10);
		RooRealVar par2("par2", "par2", 1.5, 0, 10);
		RooRealVar par3("par3", "par3", 2.0, 0, 10);
		RooRealVar par4("par4", "par4", 2.0, 0, 10);
		RooArgList alist(par1, par2, par3, par4);
		RooBernstein bkg_model_bern4("bkg_model_bern4", "bkg_model_bern4", *m2mu, alist);
		bkg_model_bern4.fitTo(data_obs);			
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

		
                //4st order polynomial
                RooRealVar lar1("lar1", "lar1", 0.0, -5.0, 5.0);
                RooRealVar lar2("lar2", "lar2", 0.0, -5.0, 5.0);
                RooRealVar lar3("lar3", "lar3", 0.0, -5.0, 5.0);
                RooRealVar lar4("lar4", "lar4", 0.0, -5.0, 5.0);
                RooArgList llist(lar1, lar2, lar3, lar4  );
                RooPolynomial bkg_model_line3("bkg_model_line3", "bkg_model_line3", *m2mu, llist, 1);
                //Exponentials
                RooRealVar car1("car1", "car1", -0.5, -7, 7);
                RooExponential bkg_model_exp3("bkg_model_exp3", "bkg_model_exp3", *m2mu, car1);
                //Product of the two
                RooProdPdf bkg_model_expxline4("bkg_model_expxline4", "bkg_model_expxline4", bkg_model_line3, bkg_model_exp3);
                bkg_model_expxline4.chi2FitTo(data_obs);

                RooCategory bkg_pdf_cat("pdf_index","Index of the background PDF which is active");
                RooArgList bkg_pdf_list;
                bkg_pdf_list.add(bkg_model_bern4);
                bkg_pdf_list.add(bkg_model_expxline4);

                RooMultiPdf bkg_model("bkg_model", "All Pdfs", bkg_pdf_cat, bkg_pdf_list);

		
		

		RooPlot *frame = m2mu->frame();
		data_obs.plotOn(frame);
		bkg_model.plotOn(frame);
		TCanvas c_all("c_all", "c_all", 800, 500);
		frame->Draw("goff");
		c_all.SaveAs(Form("output/catA_%d_"+year+".png",i));

		//save into ROO workspace
		RooWorkspace dpworkspace("dpworkspace", "");
		dpworkspace.import(data_obs);
		dpworkspace.import(*signalModel);
		dpworkspace.import(bkg_model);
		dpworkspace.writeToFile(Form("output/dpWorkspace"+year+suff+"_%d.root",i));

		//write the datacard
		char inputShape[200];
		sprintf(inputShape,"output/dpCard_"+year+suff+"_m%.3f_%d.txt",mass,i);
		ofstream newcardShape;
		newcardShape.open(inputShape);
		newcardShape << Form("imax * number of channels\n");
		newcardShape << Form("jmax * number of background\n");
		newcardShape << Form("kmax * number of nuisance parameters\n");
		newcardShape << Form("shapes data_obs	CatAB dpWorkspace"+year+suff+"_%d.root dpworkspace:data_obs\n",i);
		newcardShape << Form("shapes bkg_mass	CatAB dpWorkspace"+year+suff+"_%d.root dpworkspace:bkg_model\n",i);
		newcardShape << Form("shapes signalModel_generic	CatAB dpWorkspace"+year+suff+"_%d.root dpworkspace:signalModel_generic\n",i);
		newcardShape << Form("bin		CatAB\n");
		newcardShape << Form("observation 	-1.0\n");
		newcardShape << Form("bin     		CatAB		CatAB		\n");
		newcardShape << Form("process 		signalModel_generic  	bkg_mass	\n");
		newcardShape << Form("process 		0    		1	   	\n");
		newcardShape << Form("rate    		%f  		%f		\n",
				     effcuts*effgraph->Eval(mass,0,"S")*luminosity, catA->Integral());
		//newcardShape << Form("lumi13TeV_2017 lnN 	1.023 	-\n");
		newcardShape << Form("lumi13TeV_2018 lnN 	1.026 	-\n");
		newcardShape << Form("id_eff_mva_2018 lnN	1.10 	-\n");
		newcardShape << Form("eff_trig_2018 lnN         %f        -\n", triggSys);
		//newcardShape << Form("sig_shape_2018 lnN        1.10 	-\n");
		//newcardShape << Form("eff_mu_13TeV_2017 lnN	1.015 	-\n");
		//newcardShape << Form("bkg_norm rateParam CatA bkg_mass %f\n",catA->Integral());
		//newcardShape << Form("resA param %f %f\n",resA.getValV(),resA.getValV()*0.1);
		newcardShape.close();
		
	}
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

	TCanvas c_fVal("c_fVal", "c_fVal", 1250, 1020);
        //effValues->GetYaxis()->SetRangeUser(0.00001, 1);
        accValues->GetXaxis()->SetRangeUser(0.8, 9);
        accValues->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
        accValues->Draw("a*");
	accValues->SetTitle("");
	//plotValues->Draw("aSAME");
        auto legend = new TLegend(0.6,0.1,0.9,0.4);
        legend->AddEntry(accValues,"Acceptance","p");
        legend->AddEntry(plotValues,"Efficiency (ID + Trigger) #times Acceptance","p");
        legend->Draw();
	auto cmsTag= new TLatex(0.13,0.917,"#scale[1.1]{CMS}");
	cmsTag->SetNDC();
	cmsTag->SetTextAlign(11);
	cmsTag->Draw();
	auto cmsTag2 = new TLatex(0.215,0.917,"#scale[0.825]{#bf{#it{Preliminary}}}");
	cmsTag2->SetNDC();
	cmsTag2->SetTextAlign(11);
	cmsTag2->Draw();
        c_fVal.SetLogy();
        c_fVal.SaveAs("Multiplot.png");


	for (int j=1; j<=nbins_tsys; j++){
	  double val = 1.00 + abs(tsys->GetBinContent(j));
	  cout <<  "Tris sys " << val << "\n"; 
	}



}
