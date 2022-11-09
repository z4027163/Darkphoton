#include <iostream>
#include <TH1.h>
#include <TTree.h>
#include <math.h>

using namespace std;

void skim(){

  TString dir="/eos/user/w/wangz/darkphoton/tnp/os_pt12/";
  TString name1 = "2018_tree.root";
  TString tname = "fitter_tree";
  TFile *oldfile= new TFile(dir+name1);
  TTree *tree = (TTree*)oldfile->Get("tpTree/fitter_tree");
  int entries = tree->GetEntries();
  
  float mass, chi, dR, trkiso, drovermass,dxy,dz;
  float pt, eta, phi;
  int m2id;
  int ntklayers, nmhits, nphits;
  float vtxchi2;
  float PVd, IP;

  tree->SetBranchAddress("mass",&mass);
  tree->SetBranchAddress("ntklayers",&ntklayers);
  tree->SetBranchAddress("chi",&chi);
  tree->SetBranchAddress("dR",&dR);
  tree->SetBranchAddress("nmhits",&nmhits);
  tree->SetBranchAddress("nphits",&nphits);
  tree->SetBranchAddress("trkiso",&trkiso);
  tree->SetBranchAddress("m2id",&m2id);
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("eta",&eta);
  tree->SetBranchAddress("phi",&phi);
  tree->SetBranchAddress("dxy",&dxy);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("vtxchi2",&vtxchi2);
  tree->SetBranchAddress("PVd",&PVd);
  tree->SetBranchAddress("IP",&IP);

  float mass_c, chi_c, dR_c,trkiso_c,dxy_c, dz_c;
  int m2id_c;
  int ntklayers_c, nmhits_c, nphits_c;
  float pt_c, eta_c, phi_c;
  float vtxchi2_c;
  float PVd_c, IP_c;
  TFile *newfile = new TFile(dir+"skim_"+name1,"RECREATE");
  TDirectory *tof = newfile->mkdir("tpTree");
  TTree *newtree = new TTree(tname,tname);
  newtree->Branch("mass",&mass_c);
  newtree->Branch("ntklayers",&ntklayers_c);
  newtree->Branch("chi",&chi_c);
  newtree->Branch("dr_n",&dR_c);
  newtree->Branch("nmhits",&nmhits_c);
  newtree->Branch("nphits",&nphits_c);
  newtree->Branch("trkiso",&trkiso_c);
  newtree->Branch("m2id",&m2id_c);
  newtree->Branch("pt",&pt_c);
  newtree->Branch("eta",&eta_c);
  newtree->Branch("phi",&phi_c);
  newtree->Branch("dxy",&dxy_c);
  newtree->Branch("dz",&dz_c);
  newtree->Branch("vtxchi2",&vtxchi2_c);
  newtree->Branch("PVd",&PVd_c);
  newtree->Branch("IP_c",&IP_c);

/*
  double mean=0;
  int num=0;
  for(int i=0; i<entries; i++){
    if(i%10!=0) continue;
    tree->GetEntry(i);
    if(mass>3&&mass<3.2){
      mean+=dR/mass;
      num++;
    }
  }

  mean=mean/double(num);
*/
  for(int i=0; i<entries; i++){    
    if(i%10!=0) continue;
    tree->GetEntry(i);
//    if(mass<1 || mass>12) continue;
    if(mass<3 || mass>3.2) continue;
    if(PVd>0.05) continue;
    mass_c=mass;
    ntklayers_c=ntklayers;
    chi_c=chi;
//    dR_c=fabs(dR/mass/mean-1);
    nmhits_c=nmhits;
    nphits_c=nphits;
    trkiso_c=trkiso;
    m2id_c=m2id;
    pt_c=pt;
    eta_c=eta;
    phi_c=phi;
    vtxchi2_c=vtxchi2;
    PVd_c=PVd;
    IP_c=IP;
//    if(mass<1 || (mass>2.5 && mass<4.2) || mass>8.5) continue;
//    if(mass<60 || mass>120) continue;
    newtree->Fill();
  }
  tof->cd();
  newtree->Write(tname,TObject::kOverwrite);
  newfile->Close();
}
