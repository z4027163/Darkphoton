#include <iostream>
#include <TH1.h>
#include <TTree.h>
#include <math.h>

#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
const double mPI = 3.141592654;
float iso_c,dr_n,chi_c,dxy_c,nlayer_c,nmhits_c,nphits_c,IP_c,PVd_c, vtxchi2_c;

double DELTAPHI( double phi1, double phi2 ){
        if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
          return -999;
        }
        float dp=fabs(phi1-phi2);
        if (dp>mPI) dp-=float(2*mPI);
        return dp;
}

using namespace std;

//1.j/psi ID  2. j/psi+IPsig ID  3.upsilon0p2 ID  4. upsilon 0p2 +IP ID  5.upsilon 0p2 +PVd ID  6.upsilon 0p05 +PVd ID 7. upsilon 0p03 +PVd 8.mc ID  9.mc0p2+IPsig ID
void MVAread(TMVA::Reader *reader, int type=1){

   reader->AddVariable( "ntklayers", &nlayer_c );
   reader->AddVariable( "chi", &chi_c );
   reader->AddVariable( "nmhits", &nmhits_c );
   reader->AddVariable( "nphits", &nphits_c );
   reader->AddVariable( "trkiso", &iso_c );
   reader->AddVariable( "vtxchi2",&vtxchi2_c);
   if(type==2 || type==4 || type==9) reader->AddVariable( "IP",&IP_c);
   if(type==5 || type==6 || type==7 || type==10) reader->AddVariable("PVd",&PVd_c);

   if(type<=0 || type>10){ cout << "INPUT TYPE BUG" << endl; return;}
   TString fold[10]={"jpsi0p2","jpsi0p2_IP","upsilon0p2","upsilon0p2_IP","upsilon0p2_PVd","upsilon0p05_PVd","upsilon0p03_PVd","mc0p2_noIP","mc0p2_IP","upsilon0p015_PVd"};
   TString dir="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/MVA/darkphoton/pt4/"+fold[type-1]+"/";

   TString methodName = "BDT method";
   TString weightfile = dir + TString("TMVAClassification_BDT.weights.xml");
   reader->BookMVA( methodName, weightfile );
}

void generate_mva_9(TString treepath = "2018A/tree_1.root", TString outfilename = "./scoutData4MuonsM70.root"){

  TString dir="/eos/user/w/wangz/darkphoton/";
  TString name1 = treepath;
  TString tname = "tree";
  TFile *oldfile= new TFile(dir+name1);
  TTree *tree = (TTree*)oldfile->Get("tpTree/fitter_tree");
  int entries = tree->GetEntries();

    float charge[4];
    float  mupt[4];
    float  mueta[4];
    float  muphi[4];
    int   muid[4];
    float  much[4];
    int  muphits[4];
    int  mumhits[4];
    int  munlayer[4];
    float  muchi[4];
    float  mutrkiso[4];
    float Dxy[4], Dz[4];
    int nmuon;
    float IP, vtxchi2;
    float PVd;

    tree->SetBranchAddress("pt"  , mupt);
    tree->SetBranchAddress("eta" , mueta);
    tree->SetBranchAddress("phi" , muphi);
    tree->SetBranchAddress("trkiso",mutrkiso);
    tree->SetBranchAddress("nphits",muphits);
    tree->SetBranchAddress("nmhits",mumhits);
    tree->SetBranchAddress("muid",muid);
    tree->SetBranchAddress("chi",muchi);
    tree->SetBranchAddress("charge",charge);
    tree->SetBranchAddress("ntklayers",munlayer);
    tree->SetBranchAddress("nmuon", &nmuon);
    tree->SetBranchAddress("dxy",Dxy);
    tree->SetBranchAddress("dz",Dz);
    tree->SetBranchAddress("vtxchi2",&vtxchi2);
    tree->SetBranchAddress("IP",&IP);
    tree->SetBranchAddress("PVd",&PVd);

  float mass,m1pt,m2pt,m1eta,m2eta,m1phi,m2phi,m1id,m2id;
  float m1mva1,m1mva2,m1mva3,m1mva4,m1mva5,m1mva6,m1mva7,m1mva8,m1mva9,m1mva10;
  float m2mva1,m2mva2,m2mva3,m2mva4,m2mva5,m2mva6,m2mva7,m2mva8,m2mva9,m2mva10;
  int pair;

  TFile *newfile = new TFile(outfilename,"RECREATE");
  TTree *newtree = new TTree(tname,tname);
//  newtree->AutoSave();
//  newfile->SaveSelf();
  newtree->Branch("mass",&mass);
  newtree->Branch("m1id",&m1id);
  newtree->Branch("m1pt",&m1pt);
  newtree->Branch("m1eta",&m1eta);
  newtree->Branch("m1phi",&m1phi);
  newtree->Branch("m1mva1",&m1mva1);
  newtree->Branch("m1mva2",&m1mva2);
  newtree->Branch("m1mva3",&m1mva3);
  newtree->Branch("m1mva4",&m1mva4);
  newtree->Branch("m1mva5",&m1mva5);
  newtree->Branch("m1mva6",&m1mva6);
  newtree->Branch("m1mva7",&m1mva7);
  newtree->Branch("m1mva8",&m1mva8);
  newtree->Branch("m1mva9",&m1mva9);
  newtree->Branch("m2mva1",&m2mva1);
  newtree->Branch("m2mva2",&m2mva2);
  newtree->Branch("m2mva3",&m2mva3);
  newtree->Branch("m2mva4",&m2mva4);
  newtree->Branch("m2mva5",&m2mva5);
  newtree->Branch("m2mva6",&m2mva6);
  newtree->Branch("m2mva7",&m2mva7);
  newtree->Branch("m2mva8",&m2mva8);
  newtree->Branch("m2mva9",&m2mva9);
  newtree->Branch("m1mva10",&m1mva10);
  newtree->Branch("m2mva10",&m2mva10);
  newtree->Branch("m2id",&m2id);
  newtree->Branch("m2pt",&m2pt);
  newtree->Branch("m2eta",&m2eta);
  newtree->Branch("m2phi",&m2phi);
  newtree->Branch("pair",&pair);
  newtree->Branch("PVd",&PVd_c);
  //import mva weight


   TMVA::Reader *reader1 = new TMVA::Reader( "!Color:!Silent" );   
   MVAread(reader1, 1); 

   TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader2, 2);

   TMVA::Reader *reader3 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader3, 3);

   TMVA::Reader *reader4 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader4, 4);

   TMVA::Reader *reader5 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader5, 5);

   TMVA::Reader *reader6 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader6, 6);

   TMVA::Reader *reader7 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader7, 7);

   TMVA::Reader *reader8 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader8, 8);

   TMVA::Reader *reader9 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader9, 9);

   TMVA::Reader *reader10 = new TMVA::Reader( "!Color:!Silent" );
   MVAread(reader10, 10);
  entries=entries/10;
  cout << "entries=" << entries << endl;

 
  for(int i=0; i<entries; i++){    
    tree->GetEntry(i);
    int indx1, indx2;
    float tempt=0;
    if(nmuon<2) continue;
    //pick first 2 muons
    for(int j=0; j<nmuon; j++){
       if(tempt<mupt[j]){
          indx1=j;
          tempt=mupt[j];
       }
    }
    float tempt2=0;
    for(int j=0; j<nmuon; j++){
       if(j==indx1) continue;
       if(tempt2<mupt[j]){
          indx2=j;
          tempt2=mupt[j];
       }
    }

   m1pt=mupt[indx1];
   m2pt=mupt[indx2];
   m1eta=mueta[indx1];
   m2eta=mueta[indx2];
   m1phi=muphi[indx1];
   m2phi=muphi[indx2];
   
   pair=0;
   if(charge[indx1]*charge[indx2]>0) pair=1;
   else if(charge[indx1]*charge[indx2]<0) pair=-1;

      
   TLorentzVector mm;
	
   TLorentzVector m1;
   TLorentzVector m2;
	
   m1.SetPtEtaPhiM(m1pt, m1eta, m1phi, 0.1057);
   m2.SetPtEtaPhiM(m2pt, m2eta, m2phi, 0.1057); 
   
   mm=m1+m2;
   mass=mm.M();
   float temdr=sqrt((mueta[indx1]-mueta[indx2])*(mueta[indx1]-mueta[indx2])+DELTAPHI(muphi[indx1],muphi[indx2])*DELTAPHI(muphi[indx1],muphi[indx2]));

   float slidePt1 = mass/3.;
   float slidePt2 = mass/4.;
   if(m1pt<slidePt1 || m2pt<slidePt2) continue;
   
//indx1
   iso_c=mutrkiso[indx1];
   chi_c=muchi[indx1];
   dxy_c=Dxy[indx1];
   nlayer_c=munlayer[indx1];
   nmhits_c=mumhits[indx1];
   nphits_c=muphits[indx1];
   vtxchi2_c=vtxchi2;
   IP_c=IP; 
   PVd_c=PVd;
//   if(charge[indx1]*charge[indx2]>0) dr_n=fabs(temdr/mass/mean_ss-1);
//   if(charge[indx1]*charge[indx2]<0) dr_n=fabs(temdr/mass/mean-1); 
    
   m1mva1 = reader1->EvaluateMVA( "BDT method");
   m1mva2 = reader2->EvaluateMVA( "BDT method");
   m1mva3 = reader3->EvaluateMVA( "BDT method");
   m1mva4 = reader4->EvaluateMVA( "BDT method");
   if(PVd<0.2) m1mva5 = reader5->EvaluateMVA( "BDT method");
   else m1mva5=-1;
   if(PVd<0.05) m1mva6 = reader6->EvaluateMVA( "BDT method");
   else m1mva6=-1;
   if(PVd<0.03) m1mva7 = reader7->EvaluateMVA( "BDT method");
   else m1mva7=-1;
   m1mva8 = reader8->EvaluateMVA( "BDT method");
   m1mva9 = reader9->EvaluateMVA( "BDT method");
   if(PVd<0.015) m1mva10=reader10->EvaluateMVA( "BDT method"); 
   else m1mva10=-1;

   m1id=0;
   if(iso_c<0.15 && (nphits_c > 0) && nlayer_c > 5 && chi_c < 10.) m1id=1;       
//indx2
   iso_c=mutrkiso[indx2];
   chi_c=muchi[indx2];
   dxy_c=Dxy[indx2];
   nlayer_c=munlayer[indx2];
   nmhits_c=mumhits[indx2];
   nphits_c=muphits[indx2];

   m2mva1 = reader1->EvaluateMVA( "BDT method");
   m2mva2 = reader2->EvaluateMVA( "BDT method");
   m2mva3 = reader3->EvaluateMVA( "BDT method");
   m2mva4 = reader4->EvaluateMVA( "BDT method");
   if(PVd<0.2) m2mva5 = reader5->EvaluateMVA( "BDT method");
   else m2mva5=-1;
   if(PVd<0.05) m2mva6 = reader6->EvaluateMVA( "BDT method");
   else m2mva6=-1;
   if(PVd<0.03) m2mva7 = reader7->EvaluateMVA( "BDT method");
   else m2mva7=-1;
   m2mva8 = reader8->EvaluateMVA( "BDT method");
   m2mva9 = reader9->EvaluateMVA( "BDT method");
   if(PVd<0.015) m2mva10=reader10->EvaluateMVA( "BDT method");
   else m2mva10=-1;
 
   m2id=0;
   if(iso_c<0.15 && (nphits_c > 0) && nlayer_c > 5 && chi_c < 10.) m2id=1;

   newtree->Fill();
  }

//  newfile->cd();
//  newtree->Write(tname,TObject::kOverwrite);
  newfile->Write();
  delete reader1;
  delete reader2;
  delete reader3;
  delete reader4;
  delete reader5;
  delete reader6;
  delete reader7;
  delete reader8;
  delete reader9;
  delete reader10;
  newfile->Close();
}
