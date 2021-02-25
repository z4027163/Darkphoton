//#include "sumwgt.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include <TTreeReader.h>
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include <TTreeReaderValue.h>
#include "TLorentzVector.h"
#define mPI 3.141593
using namespace std;

float iso_c,chi_c,nlayer_c,nmhits_c,nphits_c,PVd_c, vtxchi2_c;
double DELTAPHI( double phi1, double phi2 ){
        if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
          return -999;
        }
        float dp=fabs(phi1-phi2);
        if (dp>mPI) dp-=float(2*mPI);
        return dp;
}


void PVdstudy(string treepath = "tree_2.root", const char* outfilename = "./PVdstudy.root", bool isMC=false) {

   bool debug=false; 
   vector<string> filesToRun;
   string dirIn="/eos/user/w/wangz/darkphoton/";
   filesToRun.push_back(dirIn.append(treepath));
   TFile* outfile = new TFile(outfilename, "RECREATE");
   TTree* outtree = new TTree("tree", "tree");
     



    TH1F* hpt_up = new TH1F("dimuonpt_up","upsilon dimuon pt",100,0,100);
    TH1F* hdr_up = new TH1F("dR_up","dR_up",50,0,5);
    TH1F* hip_up= new TH1F("IP_up","IP_up",80,0,8);
    TH1F* hpvd_up= new TH1F("PVd_up","PVd_up",80,0,0.2);
    TH1F* hpvd_dr_up[5];
    TH1F* hip_dr_up[5];
    TH1F* hpvd_pt_up[5];

    TH1F* hpt_jpsi = new TH1F("dimuonpt_jpsi","jpsi dimuon pt",100,0,100);
    TH1F* hdr_jpsi = new TH1F("dR_jpsi","dR_jpsi",50,0,5);
    TH1F* hpvd_jpsi= new TH1F("PVd_jpsi","PVd_jpsi",80,0,0.2);
    TH1F* hip_jpsi = new TH1F("IP_jpsi","IP_jpsi",80,0,8);
    TH1F* hpvd_dr_jpsi[5];
    TH1F* hip_dr_jpsi[5];
    TH1F* hpvd_pt_jpsi[5];

    TH1F* heta = new TH1F("eta","eta",38,-1.9,1.9);
    TH1F* heta_mass[5];

     
    double binpt[6]={0,5,10,15,20,200};
    double bindr[6]={0,0.3,0.5,0.7,0.9,5};
    double binmass[6]={1,2,4,6,8,12};
    for(int j=0;j<5;j++){
     hpvd_dr_up[j]= new TH1F(Form("up_drbin%d",j),Form("up_drbin%d",j),80,0,0.2);
     hpvd_pt_up[j]= new TH1F(Form("up_ptbin%d",j),Form("up_ptbin%d",j),80,0,0.2);
     hpvd_dr_jpsi[j]= new TH1F(Form("jpsi_drbin%d",j),Form("jpsi_drbin%d",j),80,0,0.2);
     hpvd_pt_jpsi[j]= new TH1F(Form("jpsi_ptbin%d",j),Form("jpsi_ptbin%d",j),80,0,0.2);
     hip_dr_up[j]= new TH1F(Form("ipup_drbin%d",j),Form("ipup_drbin%d",j),80,0,8);
     hip_dr_jpsi[j]= new TH1F(Form("ipjpsi_drbin%d",j),Form("ipjpsi_drbin%d",j),80,0,8);
     heta_mass[j]= new TH1F(Form("eta_massbin%d",j),Form("eta_massbin%d",j),38,-1.9,1.9);
    }

   for(unsigned int f=0; f<filesToRun.size(); f++){
    cout<<filesToRun[f]<<endl;
    TChain* chain = new TChain("fitter_tree");
    chain->Add((TString)filesToRun[f]);

    TTreeReader reader(chain);


    TTreeReaderArray<float>          charge  (reader, "charge"    );
    TTreeReaderArray<float>          mupt  (reader, "pt"    );
    TTreeReaderArray<float>          mueta  (reader, "eta"    );
    TTreeReaderArray<float>          muphi  (reader, "phi"    );
    TTreeReaderArray<int>            muid  (reader, "muid.id"    );
    TTreeReaderArray<int>         muphits  (reader, "nphits.m2phits"    );
    TTreeReaderArray<int>          mumhits  (reader, "nmhits.m2mhits"    );
    TTreeReaderArray<int>          munlayer  (reader, "ntklayers"    );
    TTreeReaderArray<float>          mutrkiso  (reader, "trkiso"    );
    TTreeReaderArray<float>          muchi  (reader, "chi"    );
    TTreeReaderValue<int>                          nmuon  (reader, "nmuon"    );
    TTreeReaderValue<float>                        vtxchi2  (reader, "vtxchi2.vtxchie2"    );
    TTreeReaderValue<float>                        PVd  (reader, "PVd"    );
    TTreeReaderValue<float>                        IP  (reader, "IP"    );
    TTreeReaderArray<float>          mva  (reader, "mva"    );

    int count[4]={0};
    while(reader.Next()) {
      count[0]++;
      if(debug) cout << "PVd=" << *PVd << endl;
//      if(*PVd>0.015) continue;
      int indx1, indx2;
      float tempt=0;
      if(*nmuon<2) continue;
      count[1]++;
      PVd_c=*PVd; 
      vtxchi2_c=*vtxchi2;
      std::vector<unsigned> goodmuons;
      for (int i = 0; i < *nmuon; i++) {
          if(mupt[i]<4 || abs(mueta[i])>1.9) continue;
//          if(mva[i]<0) continue;
//          if(muid[i]<0.5) continue;
          goodmuons.push_back(i);
      }

        if (goodmuons.size() < 2) continue;
        unsigned idx1 = goodmuons[0];
        unsigned idx2 = goodmuons[1];
        count[2]++;           
        if (mupt[idx1] < mupt[idx2]) {
            idx1 = goodmuons[1];
            idx2 = goodmuons[0];
	}
	
        TLorentzVector mm;
	
        TLorentzVector m1;
        TLorentzVector m2;

        m1.SetPtEtaPhiM(mupt[idx1], mueta[idx1], muphi[idx1], 0.1057);
        m2.SetPtEtaPhiM(mupt[idx2], mueta[idx2], muphi[idx2], 0.1057);

        mm =  m1+m2;

        float mass,m1pt,m2pt,m1eta,m2eta,m1phi,m2phi,m1id,m2id,m1ch,m2ch,m1iso,m2iso;	
        m1pt   = m1.Pt();
        m1eta  = m1.Eta();
        m1phi  = m1.Phi();
        m1iso  = mutrkiso[idx1];
        m1id   = muid[idx1];
        m1ch   = charge[idx1];

        m2pt   = m2.Pt();
        m2eta  = m2.Eta();
        m2phi  = m2.Phi();
        m2iso  = mutrkiso[idx2]; 
        m2id   = muid[idx2];
	m2ch   = charge[idx2];
	

        double weight=1;
	//2018

	//'L1_DoubleMu4_SQ_OS_dR_Max1p2', [12]
	//'L1_DoubleMu4p5_SQ_OS_dR_Max1p2', [16]
	//
        mass   = mm.M();
	//	mass4  = mmmm.M();
          
	float dimupt=mm.Pt();
        double dr= sqrt((m1eta-m2eta)*(m1eta-m2eta)+DELTAPHI(m1phi,m2phi)*DELTAPHI(m1phi,m2phi));
  
	float slidePt1 = mm.M()/3.;
	if(slidePt1<4.) slidePt1 = 4.;
	float slidePt2 = mm.M()/4.;
	if(slidePt2<4.) slidePt2 = 4.;

	float maxEta=TMath::Max(abs(m1eta),abs(m2eta));
      

	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){ 
        count[3]++;
         if(mass>9.2 && mass<9.7){
           hpvd_up->Fill(*PVd);
           hip_up->Fill(*IP); 
           hdr_up->Fill(dr);
           hpt_up->Fill(dimupt);
           for(int k=0;k<5;k++){
             if(dimupt>=binpt[k] && dimupt<binpt[k+1]) hpvd_pt_up[k]->Fill(*PVd);
             if(dr>=bindr[k] && dr<bindr[k+1]) {hpvd_dr_up[k]->Fill(*PVd); hip_dr_up[k]->Fill(*IP);}
           }
         }
        if(mass>3 && mass<3.2){
//          if(mass<3.2){
          hpvd_jpsi->Fill(*PVd);
           hdr_jpsi->Fill(dr);
           hpt_jpsi->Fill(dimupt);
           hip_jpsi->Fill(*IP);
           for(int k=0;k<5;k++){
             if(dimupt>=binpt[k] && dimupt<binpt[k+1]) hpvd_pt_jpsi[k]->Fill(*PVd);
             if(dr>=bindr[k] && dr<bindr[k+1]){hpvd_dr_jpsi[k]->Fill(*PVd); hip_dr_jpsi[k]->Fill(*IP);}
           }
        }
        for(int k=0;k<5;k++){
          if(mass>=binmass[k] && mass<binmass[k+1]) {heta_mass[k]->Fill(m1eta); heta_mass[k]->Fill(m2eta);}
        }
        if(mass>1&&mass<10) {heta->Fill(m1eta); heta->Fill(m2eta);}
	}
      goodmuons.clear();
      }
   }
   outfile->cd();
   hpvd_up->Write();
   hdr_up->Write();
   hpt_up->Write();
   hip_up->Write();
   hip_jpsi->Write();
   hpvd_jpsi->Write();
   hdr_jpsi->Write();
   hpt_jpsi->Write();
   heta->Write();
   for(int j=0;j<5;j++){
     hpvd_dr_up[j]->Write();
     hpvd_pt_up[j]->Write();
     hpvd_dr_jpsi[j]->Write();
     hpvd_pt_jpsi[j]->Write();
     hip_dr_jpsi[j]->Write();
     hip_dr_up[j]->Write();
     heta_mass[j]->Write();
   }   

    outfile->Close();

}
