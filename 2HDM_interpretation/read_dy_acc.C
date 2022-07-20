#include <iostream>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>

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

void read_dy_acc(){

   TFile *acc_file = new TFile("acceptances_dy.root");
   TEfficiency* acc_teff = (TEfficiency*)acc_file->Get("cmsacc");
   int nbins_acc=acc_teff->GetPassedHistogram()->GetNbinsX();
   double acceptances[nbins_acc];
   double m_acceptances[nbins_acc];
   for (int j=1; j<=nbins_acc; j++){
  //      if(j>9) continue;
        acceptances[j-1] = acc_teff->GetEfficiency(j);
        m_acceptances[j-1] = acc_teff->GetPassedHistogram()->GetBinCenter(j);
    //    cout << m_acceptances[j-1] << ",";
        cout << acceptances[j-1] << ",";
    //    cout << "mass=" << m_acceptances[j-1] << " acc=" << acceptances[j-1] << endl;
   }

}
