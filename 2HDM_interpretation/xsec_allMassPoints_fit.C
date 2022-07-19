#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


void xsec_allMassPoints_fit() {

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

  Double_t mass_allLimits[144] = {1.176, 1.187, 1.199, 1.211, 1.223, 1.235, 1.248, 1.26, 1.273, 1.286, 1.299, 1.312, 1.325, 1.338, 1.351, 1.365, 1.378, 1.392, 1.406, 1.42, 1.434, 1.449, 1.463,
				  1.478, 1.493, 1.508, 1.523, 1.538, 1.553, 1.569, 1.584, 1.6, 1.616, 1.632, 1.649, 1.665, 1.682, 1.699, 1.716, 1.733, 1.75, 1.768, 1.785, 1.803, 1.821, 1.839,
				  1.858, 1.876, 1.895, 1.914, 1.933, 1.953, 1.972, 1.992, 2.012, 2.032, 2.052, 2.073, 2.094, 2.114, 2.136, 2.157, 2.179, 2.2, 2.222, 2.245, 2.267, 2.29, 2.313,
				  2.336, 2.359, 2.383, 2.406, 2.43, 2.455, 2.479, 2.504, 2.529, 2.554, 2.58, 4.201, 4.243, 4.286, 4.328, 4.372, 4.415, 4.46, 4.504, 4.549, 4.595, 4.641, 4.687,
				  4.734, 4.781, 4.829, 4.877, 4.926, 4.975, 5.025, 5.075, 5.126, 5.177, 5.229, 5.282, 5.334, 5.388, 5.442, 5.496, 5.551, 5.606, 5.663, 5.719, 5.776, 5.834, 5.892,
				  5.951, 6.011, 6.071, 6.132, 6.193, 6.255, 6.318, 6.381, 6.445, 6.509, 6.574, 6.64, 6.706, 6.773, 6.841, 6.909, 6.978, 7.048, 7.119, 7.19, 7.262, 7.334, 7.408,
				  7.482, 7.557, 7.632, 7.709, 7.786, 7.864};

  for (int nbin = 0; nbin < 144; nbin++)
    {
      //cout << "mass her is :\t" << mass_allLimits[nbin] << endl;
      xsec_line = g_fit_xsec -> Eval(mass_allLimits[nbin]);
      xsec_spline = g_fit_xsec -> Eval(mass_allLimits[nbin], 0, "S");
      //cout << nbin << '\t' << "xsec_line\t" << xsec_line << '\t' << "xsec_spline\t" << xsec_spline << endl;
      cout << mass_allLimits[nbin] << '\t' << xsec_line << '\t' << xsec_spline << endl;
      xsec_l.push_back(xsec_line);
      xsec_s.push_back(xsec_spline);
    }
  //TSpline *s = (TSpline*)gROOT->FindObject("stemp");
  //xsec_spline = s -> Eval(6.6);

  g_fit_xsec -> SetTitle("HIGLU cross section of a SM-like scalar particle");
  //c_fit_xsec->Print("fit_xsec_noGrid.pdf");
  //c_fit_xsec->Print("fit_xsec_noGrid.png");

}
