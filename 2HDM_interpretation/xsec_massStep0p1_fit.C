#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


void xsec_massStep0p1_fit() {

  gStyle->SetOptFit(1);
  const int nBins = 21;
  const float fbin = 8;
  Double_t mass_xsec[nBins] = {1.1, 1.2, 1.3, 1.5, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.3, 3.5, 3.8, 4., 5., 6., 7., 8.};
  Double_t xsec_higlu[nBins] = {3037930, 1006950, 476158, 189754, 110832, 97620, 80089, 70219, 61263, 54932, 49797, 48612, 51237, 49231, 47662, 44123, 41760, 31963, 24612, 19678, 16289};

  TGraph* g_fit_xsec = new TGraph(nBins, mass_xsec, xsec_higlu);
  TCanvas *c_fit_xsec = new TCanvas ("c_fit_xsec", "c_fit_xsec", 700, 500);

  g_fit_xsec -> GetXaxis() -> SetTitle("Mass [GeV]");
  g_fit_xsec -> GetYaxis() -> SetTitle("Cross section [pb]");
 
  g_fit_xsec -> SetLineColor(kRed);
  g_fit_xsec -> SetLineWidth(2);
  g_fit_xsec -> SetMarkerColor(kBlue);
  //c1->SetFrameFillColor(41);
  //c_fit_xsec -> SetGrid();
  //c_fit_xsec -> SetLogx();
  //c_fit_xsec -> SetLogy();

  Double_t xsec_line;
  Double_t xsec_spline;
  std::vector<float> xsec_l;
  std::vector<float> xsec_s;
  
  g_fit_xsec -> Draw("*AL");

  for (float nbin = 1.1; nbin < fbin; (nbin = nbin+0.1))
    {
      //cout << "nbin here is:\t" << nbin << endl;
      xsec_line = g_fit_xsec -> Eval(nbin);
      xsec_spline = g_fit_xsec -> Eval(nbin, 0, "S");
      cout << nbin << '\t' << "xsec_line\t" << xsec_line << '\t' << "xsec_spline\t" << xsec_spline << endl;
      xsec_l.push_back(xsec_line);
      xsec_s.push_back(xsec_spline);
    }
  //TSpline *s = (TSpline*)gROOT->FindObject("stemp");
  //xsec_spline = s -> Eval(6.6);

  g_fit_xsec -> SetTitle("HIGLU cross section of a SM-like scalar particle");
  //c_fit_xsec->Print("fit_xsec_noGrid.pdf");
  //c_fit_xsec->Print("fit_xsec_noGrid.png");

}
