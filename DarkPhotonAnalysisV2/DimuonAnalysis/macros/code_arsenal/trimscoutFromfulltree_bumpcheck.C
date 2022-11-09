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
double DELTAPHI( float phi1, double phi2 ){
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
    TH1D* massforLimit_trg1_CatA[400];
    TH1D* massforLimit_trg1_CatB[400];
    TH1D* massforLimit_trg2_CatA[400];
    TH1D* massforLimit_trg2_CatB[400];
    TH1D* massforLimit_trg3_CatA[400];
    TH1D* massforLimit_trg3_CatB[400];
    TH1D* massforLimit_trg4_CatA[400];
    TH1D* massforLimit_trg4_CatB[400];
    TH1D* mmpt[400];
    float m=0.2;

    for(int j=0; j<400; j++){
      m = m+(m*0.01); 
      massforLimit_CatA[j] = new TH1D(Form("massforLimit_CatA%d",j),Form("massforLimit_CatA%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_CatA[j]->Sumw2();
      massforLimit_CatB[j] = new TH1D(Form("massforLimit_CatB%d",j),Form("massforLimit_CatB%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_CatB[j]->Sumw2();
      massforLimit_CatC[j] = new TH1D(Form("massforLimit_CatC%d",j),Form("massforLimit_CatC%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_CatC[j]->Sumw2();
      mmpt[j] = new TH1D(Form("mmpt%d",j),Form("mmpt%d",j),100,0,100);  mmpt[j]->Sumw2();
      //massforLimitBlind[j] = new TH1D(Form("massforLimitBlind%d",j),Form("massforLimitBlind%d",j),48,m-(m*0.01*6.),m+(m*0.01*6.));  massforLimitBlind[j]->Sumw2();
      massforLimit_trg1_CatA[j] = new TH1D(Form("massforLimit_trg1_CatA%d",j),Form("massforLimit_trg1_CatA%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg1_CatA[j]->Sumw2();
      massforLimit_trg1_CatB[j] = new TH1D(Form("massforLimit_trg1_CatB%d",j),Form("massforLimit_trg1_CatB%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg1_CatB[j]->Sumw2();
      massforLimit_trg2_CatA[j] = new TH1D(Form("massforLimit_trg2_CatA%d",j),Form("massforLimit_trg2_CatA%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg2_CatA[j]->Sumw2();
      massforLimit_trg2_CatB[j] = new TH1D(Form("massforLimit_trg2_CatB%d",j),Form("massforLimit_trg2_CatB%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg2_CatB[j]->Sumw2();
      massforLimit_trg3_CatA[j] = new TH1D(Form("massforLimit_trg3_CatA%d",j),Form("massforLimit_trg3_CatA%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg3_CatA[j]->Sumw2();
      massforLimit_trg3_CatB[j] = new TH1D(Form("massforLimit_trg3_CatB%d",j),Form("massforLimit_trg3_CatB%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg3_CatB[j]->Sumw2();
      massforLimit_trg4_CatA[j] = new TH1D(Form("massforLimit_trg4_CatA%d",j),Form("massforLimit_trg4_CatA%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg4_CatA[j]->Sumw2();
      massforLimit_trg4_CatB[j] = new TH1D(Form("massforLimit_trg4_CatB%d",j),Form("massforLimit_trg4_CatB%d",j),100,m-(m*0.013*8.),m+(m*0.013*8.));  massforLimit_trg4_CatB[j]->Sumw2();
      
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
    TH1F* massforLimitFull_trg1 = new TH1F("massforLimitFull_trg1","massforLimitFull_trg1",5500,0.,11.);
    TH1F* massforLimitFull_trg2 = new TH1F("massforLimitFull_trg2","massforLimitFull_trg2",5500,0.,11.);
    TH1F* massforLimitFull_trg3 = new TH1F("massforLimitFull_trg3","massforLimitFull_trg3",5500,0.,11.);
    TH1F* massforLimitFull_trg4 = new TH1F("massforLimitFull_trg4","massforLimitFull_trg4",5500,0.,11.);
    
//for study the bump
   TH1F* dz_bump = new TH1F("dz_bump","dz_bump",100,-25,25);
   TH1F* dz_side1 = new TH1F("dz_side1","dz_side1",100,-25,25);
   TH1F* dz_side2 = new TH1F("dz_side2","dz_side2",100,-25,25);
   TH1F* mmdz_bump = new TH1F("mmdz_bump","mmdz_bump",100,0,0.2);
   TH1F* mmdz_side1 = new TH1F("mmdz_side1","mmdz_side1",100,0,0.2);
   TH1F* mmdz_side2 = new TH1F("mmdz_side2","mmdz_side2",100,0,0.2);
   TH1F* mupt_bump = new TH1F("mupt_bump","mupt_bump",100,0,50);
   TH1F* mupt_side1 = new TH1F("mupt_side1","mupt_side1",100,0,50);
   TH1F* mupt_side2 = new TH1F("mupt_side2","mupt_side2",100,0,50);
   TH1F* mueta_bump = new TH1F("mueta_bump","mueta_bump",76,-1.9,1.9);
   TH1F* mueta_side1 = new TH1F("mueta_side1","mueta_side1",76,-1.9,1.9);
   TH1F* mueta_side2 = new TH1F("mueta_side2","mueta_side2",76,-1.9,1.9);
   TH1F* dxy_bump = new TH1F("dxy_bump","dxy_bump",400,-0.2,0.2);
   TH1F* dxy_side1 = new TH1F("dxy_side1","dxy_side1",400,-0.2,0.2);
   TH1F* dxy_side2 = new TH1F("dxy_side2","dxy_side2",400,-0.2,0.2);
   TH1F* mmpt_bump = new TH1F("mmpt_bump","mmpt_bump",100,0,50);
   TH1F* mmpt_side1 = new TH1F("mmpt_side1","mmpt_side1",100,0,50);
   TH1F* mmpt_side2 = new TH1F("mmpt_side2","mmpt_side2",100,0,50);
   TH1F* dr_bump = new TH1F("dr_bump","dr_bump",100,0,0.5);
   TH1F* dr_side1 = new TH1F("dr_side1","dr_side1",100,0,0.5);
   TH1F* dr_side2 = new TH1F("dr_side2","dr_side2",100,0,0.5);
   TH1F* drm_bump = new TH1F("drm_bump","drm_bump",100,0,0.4);
   TH1F* drm_side1 = new TH1F("drm_side1","drm_side1",100,0,0.4);
   TH1F* drm_side2 = new TH1F("drm_side2","drm_side2",100,0,0.4);
   TH1F* ip_bump = new TH1F("ip_bump","ip_bump",35,0,3.5);
   TH1F* ip_side1 = new TH1F("ip_side1","ip_side1",35,0,3.5);
   TH1F* ip_side2 = new TH1F("ip_side2","ip_side2",35,0,3.5);

 
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
    TTreeReaderArray<float>          dxy  (reader, "dxy");
    TTreeReaderArray<float>          dz  (reader, "dz");
    TTreeReaderValue<unsigned int>            l1trg (reader, "l1trg");
    TTreeReaderValue<float>                        PVd  (reader, "PVd"    );
    TTreeReaderArray<float>          mva  (reader, "mva2"    );  //mva:MVA7, mva2:MVAlow
    TTreeReaderValue<float>          IP   (reader, "IP");
    int count[4]={0};
    while(reader.Next()) {
      count[0]++;
      if(debug) cout << "PVd=" << *PVd << endl;
      if(*PVd>0.2) continue;
//      if(*PVd>0.015) continue;     //MVA7 selection
      if(*IP>3.5) continue;      //MVAlow selection
//      if(*IP<3.5 || *IP>11) continue;
      int indx1, indx2;
      float tempt=0;
      if(*nmuon<2) continue;
      //if(*nmuon!=2) continue;
      count[1]++;
      PVd_c=*PVd; 
      vtxchi2_c=*vtxchi2;
      std::vector<unsigned> goodmuons;
      for (int i = 0; i < *nmuon; i++) {
          if(mupt[i]<4 || abs(mueta[i])>1.9) continue;
          if(mva[i]<-0.1) continue;     //MVAlow
//          if(fabs(dxy[i])>0.015) continue; 
//          if(fabs(dz[i])>5) continue;
//          if(mva[i]<-0.0) continue;       //MVA7
//          if(muid[i]<0.5) continue;
//reverse id
//          if(mva[i]>-0.1) continue;     //MVAlow
//          if(mva[i]<-0.0) continue;       //MVA7
          
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
        float m1dz,m2dz,m1dxy,m2dxy;	
        m1pt   = m1.Pt();
        m1eta  = m1.Eta();
        m1phi  = m1.Phi();
        m1iso  = mutrkiso[idx1];
        m1id   = muid[idx1];
        m1ch   = charge[idx1];
        m1dz   = dz[idx1];
        m1dxy  = dxy[idx1];
        

        m2pt   = m2.Pt();
        m2eta  = m2.Eta();
        m2phi  = m2.Phi();
        m2iso  = mutrkiso[idx2]; 
        m2id   = muid[idx2];
	m2ch   = charge[idx2];
	m2dz   = dz[idx2];
        m2dxy  = dxy[idx2];

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
        float dr= sqrt((m1eta-m2eta)*(m1eta-m2eta)+DELTAPHI(m1phi,m2phi)*DELTAPHI(m1phi,m2phi));
        float drm = dr/mass;

        if((*l1trg&4)==4 && mass<3) continue;
	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2) forLimitMassInclZ->Fill(mass,weight);
//	if(!( l1Result->at(2) || l1Result->at(6) || l1Result->at(7) || l1Result->at(10) || l1Result->at(11) || l1Result->at(13))) continue;

	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 ) forLimitMassNonInclZ->Fill(mass,weight);


	if(m1ch*m2ch<0. &&  m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9) forLimitMassZ->Fill(mass,weight);

	if(m1ch*m2ch<0. && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){ 
        count[3]++;

	massforLimitFull->Fill(mass,weight); 

        if((*l1trg&1)==1) massforLimitFull_trg1->Fill(mass,weight);
        if((*l1trg&2)==2) massforLimitFull_trg2->Fill(mass,weight);
        if((*l1trg&4)==4) massforLimitFull_trg3->Fill(mass,weight);
        if((*l1trg&8)==8) massforLimitFull_trg4->Fill(mass,weight);
      //TEST
       // if((*l1trg&4)==4 && (*l1trg&2)!=2 && (*l1trg&8)!=8) massforLimitFull_trg1->Fill(mass,weight);
       // if((*l1trg&4)==4 && ((*l1trg&2)==2 || (*l1trg&8)==8)) massforLimitFull_trg3->Fill(mass,weight);
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

        //for bump check
        if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){
          if(fabs(mass-1.584)<0.03){
              dz_bump->Fill(m1dz,weight);
              dz_bump->Fill(m2dz,weight);
              mmdz_bump->Fill(fabs(m1dz-m2dz),weight);
              mupt_bump->Fill(m1pt,weight);
              mupt_bump->Fill(m2pt,weight);
              mueta_bump->Fill(m1eta,weight);
              mueta_bump->Fill(m2eta,weight);
              dxy_bump->Fill(m1dxy,weight);
              dxy_bump->Fill(m2dxy,weight);
              mmpt_bump->Fill(ptll,weight);
              dr_bump->Fill(dr,weight);
              drm_bump->Fill(drm,weight);
              ip_bump->Fill(*IP,weight);

          }
          else if(mass>1.524 && mass<1.554){
              dz_side1->Fill(m1dz,weight);
              dz_side1->Fill(m2dz,weight);
              mmdz_side1->Fill(fabs(m1dz-m2dz),weight);
              mupt_side1->Fill(m1pt,weight);
              mupt_side1->Fill(m2pt,weight);
              mueta_side1->Fill(m1eta,weight);
              mueta_side1->Fill(m2eta,weight);
              dxy_side1->Fill(m1dxy,weight);
              dxy_side1->Fill(m2dxy,weight);
              mmpt_side1->Fill(ptll,weight);
              dr_side1->Fill(dr,weight);
              drm_side1->Fill(drm,weight);
              ip_side1->Fill(*IP,weight);
          }
          else if(mass>1.614 && mass<1.644){
              dz_side2->Fill(m1dz,weight);
              dz_side2->Fill(m2dz,weight);
              mmdz_side2->Fill(fabs(m1dz-m2dz),weight);
              mupt_side2->Fill(m1pt,weight);
              mupt_side2->Fill(m2pt,weight);
              mueta_side2->Fill(m1eta,weight);
              mueta_side2->Fill(m2eta,weight);
              dxy_side2->Fill(m1dxy,weight);
              dxy_side2->Fill(m2dxy,weight);
              mmpt_side2->Fill(ptll,weight);
              dr_side2->Fill(dr,weight);
              drm_side2->Fill(drm,weight);
              ip_side2->Fill(*IP,weight);
          }
        }
	float ma=0.2;
	for(int j=0; j<400.; j++){
      		ma = ma+(ma*0.01); 

     		if(mass > ma-(ma*0.013*8.) && mass < ma+(ma*0.013*8.)) {

		  	if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9){ massforLimit_CatA[j]->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.9 ){ massforLimit_CatB[j]->Fill(mass,weight); }
			if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 ){ massforLimit_CatC[j]->Fill(mass,weight); }
                        if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9){mmpt[j]->Fill(ptll,weight);}

                        //trg separation
                        if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9){
                            if((*l1trg&1)==1) massforLimit_trg1_CatA[j]->Fill(mass,weight);
                            if((*l1trg&2)==2) massforLimit_trg2_CatA[j]->Fill(mass,weight);
                            if((*l1trg&4)==4) massforLimit_trg3_CatA[j]->Fill(mass,weight);
                            if((*l1trg&8)==8) massforLimit_trg4_CatA[j]->Fill(mass,weight);
                        }
                        if(m1ch*m2ch<0.  && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.9 ){
                            if((*l1trg&1)==1) massforLimit_trg1_CatB[j]->Fill(mass,weight);
                            if((*l1trg&2)==2) massforLimit_trg2_CatB[j]->Fill(mass,weight);
                            if((*l1trg&4)==4) massforLimit_trg3_CatB[j]->Fill(mass,weight);
                            if((*l1trg&8)==8) massforLimit_trg4_CatB[j]->Fill(mass,weight);
                        }
		}
      
   	 }	
	
      goodmuons.clear();
      }
      cout << "final showdown " << "count1=" << count[0] << " count2=" << count[1] << " count3=" << count[2] << " count4=" << count[3] << endl;
   }
   outfile->cd();
   massforLimitFull->Write();
   massforLimitFull_trg1->Write();
   massforLimitFull_trg2->Write();
   massforLimitFull_trg3->Write();
   massforLimitFull_trg4->Write();
 
   dz_bump->Write();
   mupt_bump->Write();
   mueta_bump->Write();
   dz_side1->Write();
   mupt_side1->Write();
   mueta_side1->Write();
   dz_side2->Write();
   mupt_side2->Write();
   mueta_side2->Write(); 
   dxy_bump->Write();
   dxy_side1->Write();
   dxy_side2->Write();
   mmpt_bump->Write();
   mmpt_side1->Write();
   mmpt_side2->Write();
   mmdz_bump->Write();
   mmdz_side1->Write();
   mmdz_side2->Write();
   dr_bump->Write();
   dr_side1->Write();
   dr_side2->Write();
   drm_bump->Write();
   drm_side1->Write();
   drm_side2->Write();
   ip_bump->Write();
   ip_side1->Write();
   ip_side2->Write();
   for(int j=0; j<400.;j++){
	massforLimit_CatA[j]->Write();
     	massforLimit_CatB[j]->Write();
     	massforLimit_CatC[j]->Write();
      //  mmpt[j]->Write();
/*
        massforLimit_trg1_CatA[j]->Write();
        massforLimit_trg1_CatB[j]->Write();
        massforLimit_trg2_CatA[j]->Write();
        massforLimit_trg2_CatB[j]->Write();
        massforLimit_trg3_CatA[j]->Write();
        massforLimit_trg3_CatB[j]->Write();
        massforLimit_trg4_CatA[j]->Write();
        massforLimit_trg4_CatB[j]->Write();
*/
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
