//#include "sumwgt.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include <TTreeReader.h>
#include "TH1D.h"
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


void trimscoutFromfulltree_bin(string treepath = "full_tem.root", const char* outfilename = "./scoutDataMVA_v2.root", bool isMC=false) {

   bool debug=false; 
   vector<string> filesToRun;
//   string dirIn="/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/";
   string dirIn="/eos/user/w/wangz/darkphoton/";
   filesToRun.push_back(dirIn.append(treepath));
   TFile* outfile = new TFile(outfilename, "RECREATE");
   TTree* outtree = new TTree("tree", "tree");
     

    TH1D* massforLimit_CatA[400];
    TH1D* massforLimit_CatB[400];
    TH1D* massforLimit_CatC[400];
    TH1D* mmpt[400];
    float m=0.2;

    for(int j=0; j<400; j++){
      m = m+(m*0.01); 
      massforLimit_CatA[j] = new TH1D(Form("massforLimit_CatA%d",j),Form("massforLimit_CatA%d",j),100,m-(m*0.013*5.),m+(m*0.013*5.));  massforLimit_CatA[j]->Sumw2();
      massforLimit_CatB[j] = new TH1D(Form("massforLimit_CatB%d",j),Form("massforLimit_CatB%d",j),100,m-(m*0.013*5.),m+(m*0.013*5.));  massforLimit_CatB[j]->Sumw2();
      massforLimit_CatC[j] = new TH1D(Form("massforLimit_CatC%d",j),Form("massforLimit_CatC%d",j),100,m-(m*0.013*5.),m+(m*0.013*5.));  massforLimit_CatC[j]->Sumw2();
      mmpt[j] = new TH1D(Form("mmpt%d",j),Form("mmpt%d",j),100,0,100);  mmpt[j]->Sumw2();
      //massforLimitBlind[j] = new TH1D(Form("massforLimitBlind%d",j),Form("massforLimitBlind%d",j),48,m-(m*0.01*6.),m+(m*0.01*6.));  massforLimitBlind[j]->Sumw2();
      
    }
    
    

    
    TH2F* BS = new TH2F("BS","BS",50,-1.,1.,50,-1.,1.);

    TH1F* forLimitMassZ = new TH1F("forLimitMassZ","forLimitMassZ",1200,0.2,120.);
    TH1F* forLimitMassInclZ = new TH1F("forLimitMassInclZ","forLimitMassInclZ",1200,0.2,120.);
    TH1F* forLimitMassNonInclZ = new TH1F("forLimitMassNonInclZ","forLimitMassNonInclZ",1200,0.2,120.);
    TH1F* forResolutionAMassZ = new TH1F("forResolutionAMassZ","forResolutionAMassZ",240,60.,120.);
    TH1F* forResolutionBMassZ = new TH1F("forResolutionBMassZ","forResolutionBMassZ",240,60.,120.);
    TH1F* forResolutionCMassZ = new TH1F("forResolutionCMassZ","forResolutionCMassZ",240,60.,120.);
    TH1F* forResolutionAMassJPsi = new TH1F("forResolutionAMassJPsi","forResolutionAMassJPsi",45,2.6,3.5);
    TH1F* forResolutionBMassJPsi = new TH1F("forResolutionBMassJPsi","forResolutionBMassJPsi",45,2.6,3.5);
    TH1F* forResolutionCMassJPsi = new TH1F("forResolutionCMassJPsi","forResolutionCMassJPsi",45,2.6,3.5);


    TH1F* massforLimitFull = new TH1F("massforLimitFull","massforLimitFull",5500,0.,11.);
    TH1F* massforLimitEta = new TH1F("massforLimitEta","massforLimitEta",260,0.490,0.620);
    TH1F* massforLimitOmega = new TH1F("massforLimitOmega","massforLimitOmega",400,0.68,0.88);
    TH1F* massforLimitPhi = new TH1F("massforLimitPhi","massforLimitPhi",480,0.90,1.14);
    TH1F* massforLimitJPsiPsi = new TH1F("massforLimitJPsiPsi","massforLimitJPsiPsi",3000,2.7,4.2);
    TH1F* massforLimitUpsilon = new TH1F("massforLimitUpsilon","massforLimitUpsilon",5000,8.5,11.);

    TH1F* massforLimitFullTight = new TH1F("massforLimitFullTight","massforLimitFullTight",5000,0.,10.);
    TH1F* massforLimitEtaTight = new TH1F("massforLimitEtaTight","massforLimitEtaTight",260,0.490,0.620);
    TH1F* massforLimitOmegaTight = new TH1F("massforLimitOmegaTight","massforLimitOmegaTight",400,0.68,0.88);
    TH1F* massforLimitPhiTight = new TH1F("massforLimitPhiTight","massforLimitPhiTight",480,0.90,1.14);
    TH1F* massforLimitJPsiPsiTight = new TH1F("massforLimitJPsiPsiTight","massforLimitJPsiPsiTight",3000,2.7,4.2);
    TH1F* massforLimitUpsilonTight = new TH1F("massforLimitUpsilonTight","massforLimitUpsilonTight",5000,8.5,11.);

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
    TTreeReaderArray<float>          mva  (reader, "mva"    );  //mva:MVA7, mva2:MVAlow
    TTreeReaderValue<float>          IP   (reader, "IP");

    int count[4]={0};
    while(reader.Next()) {
      count[0]++;
      if(debug) cout << "PVd=" << *PVd << endl;
      if(*PVd>0.015) continue;     //MVA7 selection
//      if(*IP>3.5) continue;      //MVAlow selection
      int indx1, indx2;
      float tempt=0;
      if(*nmuon<2) continue;
      count[1]++;
      PVd_c=*PVd; 
      vtxchi2_c=*vtxchi2;
      std::vector<unsigned> goodmuons;
      for (int i = 0; i < *nmuon; i++) {
          if(mupt[i]<4 || abs(mueta[i])>1.9) continue;
//          if(mva[i]<-0.1) continue;     //MVAlow
          if(mva[i]<-0.0) continue;       //MVA7
//          if(muid[i]<0.5) continue;
          goodmuons.push_back(i);
      }

        if (goodmuons.size() < 2) continue;
        unsigned idx1 = goodmuons[0];
        unsigned idx2 = goodmuons[1];
        if(debug) cout << "idx1,2=" << idx1 << "," << idx2 << endl; 
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
        double ptll=mm.Pt(); 
	
	float slidePt1 = mm.M()/3.;
	if(slidePt1<4.) slidePt1 = 4.;
	float slidePt2 = mm.M()/4.;
	if(slidePt2<4.) slidePt2 = 4.;

	float maxEta=TMath::Max(abs(m1eta),abs(m2eta));



	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2) forLimitMassInclZ->Fill(mass,weight);
//	if(!( l1Result->at(2) || l1Result->at(6) || l1Result->at(7) || l1Result->at(10) || l1Result->at(11) || l1Result->at(13))) continue;

	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 ) forLimitMassNonInclZ->Fill(mass,weight);


	if(m1ch*m2ch<0. &&  m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9) forLimitMassZ->Fill(mass,weight);

	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){ 
        count[3]++;

	massforLimitFull->Fill(mass,weight); 
	massforLimitEta->Fill(mass,weight);
	massforLimitOmega->Fill(mass,weight);
	massforLimitPhi->Fill(mass,weight);
	massforLimitJPsiPsi->Fill(mass,weight);
	massforLimitUpsilon->Fill(mass,weight);


	massforLimitFullTight->Fill(mass,weight); 
	massforLimitEtaTight->Fill(mass,weight);
	massforLimitOmegaTight->Fill(mass,weight);
	massforLimitPhiTight->Fill(mass,weight);
	massforLimitJPsiPsiTight->Fill(mass,weight);
	massforLimitUpsilonTight->Fill(mass,weight);

	
	}


	
	if(mass>60. && mass<120.){
			if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9){ forResolutionAMassZ->Fill(mass,weight); }
			if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.2 ){ forResolutionBMassZ->Fill(mass,weight); }
			if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 ){ forResolutionCMassZ->Fill(mass,weight); }

	}


	if(mass>2.6 && mass<3.5){
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9 ){ forResolutionAMassJPsi->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.2){ forResolutionBMassJPsi->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4){ forResolutionCMassJPsi->Fill(mass,weight); }

	}


	float ma=0.2;
	for(int j=0; j<400.; j++){
      		ma = ma+(ma*0.01); 

     		if(mass > ma-(ma*0.013*5.) && mass < ma+(ma*0.013*5.)) {

		  	if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9){ massforLimit_CatA[j]->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.9 ){ massforLimit_CatB[j]->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 ){ massforLimit_CatC[j]->Fill(mass,weight); }
                        if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){mmpt[j]->Fill(ptll,weight);}
		}
      
   	 }	
	
      goodmuons.clear();
      }
      cout << "final showdown " << "count1=" << count[0] << " count2=" << count[1] << " count3=" << count[2] << " count4=" << count[3] << endl;
   }
   outfile->cd();
   massforLimitFull->Write();
   massforLimitEta->Write();
   massforLimitOmega->Write();
   massforLimitPhi->Write();
   massforLimitJPsiPsi->Write();
   massforLimitUpsilon->Write();

   massforLimitFullTight->Write();
   massforLimitEtaTight->Write();
   massforLimitOmegaTight->Write();
   massforLimitPhiTight->Write();
   massforLimitJPsiPsiTight->Write();
   massforLimitUpsilonTight->Write();


   for(int j=0; j<400.;j++){
	massforLimit_CatA[j]->Write();
     	massforLimit_CatB[j]->Write();
     	massforLimit_CatC[j]->Write();
        mmpt[j]->Write();
     }
   forLimitMassZ->Write();
   forLimitMassInclZ->Write();
   forLimitMassNonInclZ->Write();
   forResolutionAMassZ->Write();
   forResolutionBMassZ->Write();
   forResolutionCMassZ->Write();
   forResolutionAMassJPsi->Write();
   forResolutionBMassJPsi->Write();
   forResolutionCMassJPsi->Write();
    outfile->Close();

}
