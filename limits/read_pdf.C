#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;
 
void read_pdf(){

        TFile* f_ws = TFile::Open("output_dual/dpWorkspace2018IterV3_202.root", "READ");
        RooWorkspace *w = (RooWorkspace*)f_ws->Get("dpworkspace");
        RooRealVar *x = (RooRealVar*)w->var("m2mu");

        RooVoigtian* dkkModel = (RooVoigtian*)w->pdf("bkg_model_2018");
        RooPlot *frame = x->frame(Title("peak fit"));
        dkkModel->plotOn(frame);

        TCanvas *c = new TCanvas("c","c");
        c->cd();
        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->Draw();
}

