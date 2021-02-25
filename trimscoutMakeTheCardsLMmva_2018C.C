//#include "sumwgt.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

float iso_c,chi_c,nlayer_c,nmhits_c,nphits_c,PVd_c, vtxchi2_c;
void MVAread(TMVA::Reader *reader, int type=1){

   reader->AddVariable( "ntklayers", &nlayer_c );
   reader->AddVariable( "chi", &chi_c );
   reader->AddVariable( "nmhits", &nmhits_c );
   reader->AddVariable( "nphits", &nphits_c );
   reader->AddVariable( "trkiso", &iso_c );
   reader->AddVariable( "vtxchi2",&vtxchi2_c);
   if(type>2) reader->AddVariable("PVd",&PVd_c);
   if(type<=0 || type>4){ cout << "INPUT TYPE BUG" << endl; return;}
   TString fold[4]={"JPsi0p2_noIP","upsilon0p2","upsilon0p2_PVd","upsilon0p03_PVd"};
//   TString dir="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/MVA/darkphoton/pt4/"+fold[type-1]+"/";
   TString dir="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/MVA/darkphoton/pt4/"+fold[type-1]+"/";
   TString methodName = "BDT method";
   TString weightfile = dir + TString("TMVAClassification_BDT.weights.xml");
   reader->BookMVA( methodName, weightfile );
}

void trimscoutMakeTheCardsLMmva_2018C(string treepath = "scout_2.root", const char* outfilename = "./scoutData4MuonsM70.root", bool isMC=false) {

   bool debug=false; 
   vector<string> filesToRun;
   string dirIn="/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/";
   filesToRun.push_back(dirIn.append(treepath));
   TFile* outfile = new TFile(outfilename, "RECREATE");
//   TTree* outtree = new TTree("tree", "tree");
     

    TH1D* massforLimit_CatA[400];
    TH1D* massforLimit_CatB[400];
    TH1D* massforLimit_CatC[400];
    
    float m=0.2;

    for(int j=0; j<400; j++){
      m = m+(m*0.01); 
      massforLimit_CatA[j] = new TH1D(Form("massforLimit_CatA%d",j),Form("massforLimit_CatA%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatA[j]->Sumw2();
      massforLimit_CatB[j] = new TH1D(Form("massforLimit_CatB%d",j),Form("massforLimit_CatB%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatB[j]->Sumw2();
      massforLimit_CatC[j] = new TH1D(Form("massforLimit_CatC%d",j),Form("massforLimit_CatC%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatC[j]->Sumw2();

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
    cout << filesToRun[f].c_str() << endl;
    TFile *fin = new TFile(filesToRun[f].c_str());
    TTree *intree = (TTree*)fin->Get("mmtree/tree");
    int entries = intree->GetEntries();

    std::vector<float>          *mpt;    
    std::vector<float>          *meta;   
    std::vector<float>         *mphi;  
    std::vector<float>          *mcharge; 
    std::vector<char>          *mid;   
    std::vector<float>          *chi2;  
    std::vector<float>          *cpiso;  
    std::vector<float>          *chiso;    
    std::vector<float>          *phiso;   
    std::vector<float>          *nhiso;  
    std::vector<float>          *puiso;   
    std::vector<float>          *tkiso; 
    std::vector<bool>           *l1Result;
    std::vector<unsigned char>  *nmhits; 
    std::vector<unsigned char>  *nphits;  
    std::vector<unsigned char>  *ntklayers;
    UChar_t  hlt; 
    UInt_t  nverts;   
    std::vector<float>           *vtxX;
    std::vector<float>           *vtxY;   
    std::vector<float>           *vtxZ;
    std::vector<float>           *vtxchi2; 
    double                      rho;
    unsigned int               run;
    unsigned int                 lumSec;
  
    double *wgt;
    double *xsec; 

    TBranch *b_mpt=0;
    TBranch *b_meta=0;
    TBranch *b_mphi=0;
    TBranch *b_mcharge=0;
    TBranch *b_mid=0;
    TBranch *b_chi2=0;
    TBranch *b_tkiso=0;
    TBranch *b_l1Result=0;
    TBranch *b_nmhits=0;
    TBranch *b_nphits=0;
    TBranch *b_ntklayers=0;
    TBranch *b_hlt=0;
    TBranch *b_nverts=0;
    TBranch *b_vtxX=0;
    TBranch *b_vtxY=0;
    TBranch *b_vtxZ=0;
    TBranch *b_vtxchi2=0;
    TBranch *b_rho=0;
    TBranch *b_run=0;
    TBranch *b_lumSec=0;
    TBranch *b_wgt=0;
    TBranch *b_xsec=0;



    intree->SetBranchAddress("muonpt",&mpt,&b_mpt);
    intree->SetBranchAddress("muoneta",&meta, &b_meta);
    intree->SetBranchAddress("muonphi",&mphi,&b_mphi);
    intree->SetBranchAddress("muoncharge",&mcharge,&b_mcharge);
    intree->SetBranchAddress("muonid",&mid,&b_mid);
    intree->SetBranchAddress("chi2",&chi2,&b_chi2);
    intree->SetBranchAddress("tkiso",&tkiso,&b_tkiso);
    intree->SetBranchAddress("l1Result",&l1Result,&b_l1Result);
    intree->SetBranchAddress("nMuonHits",&nmhits,&b_nmhits);
    intree->SetBranchAddress("nPixelHits",&nphits,&b_nphits);
    intree->SetBranchAddress("nTkLayers",&ntklayers,&b_ntklayers);
    intree->SetBranchAddress("trig",&hlt,&b_hlt);
    intree->SetBranchAddress("vtxChi2",&vtxchi2,&b_vtxchi2);
    intree->SetBranchAddress("nvtx",&nverts,&b_nverts);
    intree->SetBranchAddress("vtxX",&vtxX,&b_vtxX);
    intree->SetBranchAddress("vtxY",&vtxY,&b_vtxY);
    intree->SetBranchAddress("rho",&rho,&b_rho);
    intree->SetBranchAddress("run",&run,&b_run);
    intree->SetBranchAddress("lumSec",&lumSec,&b_lumSec);
    
//   if (isMC) {
//         intree->SetBranchAddress("lumSec",&lumSec,&b_lumSec);	
//    }
	 
 

    double wgtsum  = 1.0;//isMC ? sumwgt(treepath) : 1.0;
    float  puwgt   = 1.0;
    float  mcwgt   = 1.0;
    float  m1pt    = 0.0;        
    float  m1eta   = 0.0;        
    float  m1phi   = 0.0;        
    float  m1iso   = 0.0;        
    float  m2pt    = 0.0;        
    float  m2eta   = 0.0;        
    float  m2phi   = 0.0;        
    float  m2iso   = 0.0;        
    float  mass    = 0.0;        
    float  mass4   = 0.0;        
    char   m1id    = 0;
    char   m2id    = 0;
    char   m3id    = 0;
    char   m4id    = 0;
    float  m1ch    = 0.; 
    float  m2ch    = 0.;
    float  m3ch    = 0.; 
    float  m4ch    = 0.; 
    vector<bool> l1Bools;
    unsigned nvtx  = 0;
    unsigned Run   = 0;
    unsigned LumSec   = 0;
/*
    if (isMC) {
    outtree->Branch("mcwgt" , &mcwgt , "mcwgt/F" );
    outtree->Branch("puwgt" , &puwgt , "puwgt/F" );
    }
    outtree->Branch("m1pt"  , &m1pt  , "m1pt/F"  );
    outtree->Branch("m1eta" , &m1eta , "m1eta/F" );
    outtree->Branch("m1phi" , &m1phi , "m1phi/F" );
    outtree->Branch("m1iso" , &m1iso , "m1iso/F" );
    outtree->Branch("m1id"  , &m1id  , "m1id/B"  );
    outtree->Branch("m2pt"  , &m2pt  , "m2pt/F"  );
    outtree->Branch("m2eta" , &m2eta , "m2eta/F" );
    outtree->Branch("m2phi" , &m2phi , "m2phi/F" );
    outtree->Branch("m2iso" , &m2iso , "m2iso/F" );
    outtree->Branch("m2id"  , &m2id  , "m2id/B"  );
    outtree->Branch("mass"  , &mass  , "mass/F"  );
    outtree->Branch("nvtx"  , &nvtx  , "nvtx/i"  );
    outtree->Branch("Run"   , &Run   , "Run/i"  );
    outtree->Branch("LumSec", &LumSec, "LumSec/i"  );

    //outtree->Branch("l1Bools"  , &l1Bools , "l1Bools" );
*/
    int p=0;
    TMVA::Reader *reader1 = new TMVA::Reader( "!Color:!Silent" );
    MVAread(reader1, 4);
    if(debug) cout << "total=" << entries << endl;
//    if(entries>100000) entries=100000;
//    entries=entries/50;
    for(int t=0; t<entries; t++){
      intree->GetEntry(t);
      if(t%10000==0) cout << "entry=" << t << endl;
      if ((int(hlt) & 2) == 0) continue;   
      if(!( l1Result->at(2) ||l1Result->at(7) || l1Result->at(11) || l1Result->at(6))) continue;
      bool passIso=false;
      bool passIsoLoose=false;

      double PVd=1;
      double vchi2=0;
      nvtx=nverts;
      float BSx = 0.092;
      float BSy = -0.06;
      if(treepath.find("_13TeV") != string::npos){ BSx = -0.029; BSy = 0.07;}
      if(isMC){BSx=0.0107; BSy=0.0417;}

      if(nvtx>0. ){
	PVd= sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) );
        vchi2=(*vtxchi2)[0];
      }
      if(PVd>0.015) continue;
      if(mpt->size()<2) continue;
      PVd_c=PVd; 
      vtxchi2_c=vchi2;
      double ea = (isMC ? 0.3 : 0.45);
      std::vector<int> goodmuons;
      if(debug) cout << "nmuon=" << mpt->size() << endl;
      for (int i = 0; i < mpt->size(); i++) {
          if((*mpt)[i]<4) continue;
          if(abs((*meta)[i])>1.9) continue;
	  iso_c=(*tkiso)[i];
          chi_c=(*chi2)[i];
          nlayer_c=int((*ntklayers)[i]);
          nmhits_c=int((*nmhits)[i]);
          nphits_c=int((*nphits)[i]);
          double mva=reader1->EvaluateMVA( "BDT method");
          if(debug) cout << "muon["<<i<<"] mva=" << mva << endl;
          if(mva<0) continue;
          goodmuons.push_back(i);
      }

        if (goodmuons.size() < 2) continue;
        int idx1 = goodmuons.at(0);
        int idx2 = goodmuons.at(1);
        if(debug) cout << "idx1,2=" << idx1 << "," << idx2 << endl; 
                
	if((*tkiso)[idx1]<0.02 && (*tkiso)[idx2]<0.02 ) passIso=true;
	if((*tkiso)[idx1]<0.15 && (*tkiso)[idx2]<0.15 ) passIsoLoose=true;
        passIsoLoose=true;

        if ((*mpt)[idx1] < (*mpt)[idx2]) {
            idx1 = goodmuons.at(1);
            idx2 = goodmuons.at(0);
	    }
	
        TLorentzVector mm;
	
        TLorentzVector m1;
        TLorentzVector m2;

        m1.SetPtEtaPhiM((*mpt)[idx1], (*meta)[idx1], (*mphi)[idx1], 0.1057);
        m2.SetPtEtaPhiM((*mpt)[idx2], (*meta)[idx2], (*mphi)[idx2], 0.1057);

        mm =  m1+m2;
	
        m1pt   = m1.Pt();
        m1eta  = m1.Eta();
        m1phi  = m1.Phi();
        m1iso  = (*tkiso)[idx1];
        m1id   = (*mid)[idx1];
        m1ch   = (*mcharge)[idx1];

        m2pt   = m2.Pt();
        m2eta  = m2.Eta();
        m2phi  = m2.Phi();
        m2iso  = (*tkiso)[idx2]; 
        m2id   = (*mid)[idx2];
	m2ch   = (*mcharge)[idx2];
	
	Run=run;
	LumSec=lumSec;

	//2018

	//'L1_DoubleMu4_SQ_OS_dR_Max1p2', [12]
	//'L1_DoubleMu4p5_SQ_OS_dR_Max1p2', [16]
	//
        mass   = mm.M();
	//	mass4  = mmmm.M();
        
	
        bool passPVconstraint = false;
	bool passPVconstraintTight = false;
	float slidePt1 = mm.M()/3.;
	if(slidePt1<4.) slidePt1 = 4.;
	float slidePt2 = mm.M()/4.;
	if(slidePt2<4.) slidePt2 = 4.;

	float maxEta=TMath::Max(abs(m1eta),abs(m2eta));

	
	

	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.2 ) passPVconstraint = true;
	  }

	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.1 ) passPVconstraintTight = true;
	  }

	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && passIsoLoose) forLimitMassInclZ->Fill(mass);
//	if(!( l1Result->at(2) || l1Result->at(6) || l1Result->at(7) || l1Result->at(10) || l1Result->at(11) || l1Result->at(13))) continue;

	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && passIsoLoose) forLimitMassNonInclZ->Fill(mass);


	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9 && passIsoLoose) forLimitMassZ->Fill(mass);

	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9 && passIsoLoose){ 

	massforLimitFull->Fill(mass); 
	massforLimitEta->Fill(mass);
	massforLimitOmega->Fill(mass);
	massforLimitPhi->Fill(mass);
	massforLimitJPsiPsi->Fill(mass);
	massforLimitUpsilon->Fill(mass);

	if(passPVconstraintTight){

	massforLimitFullTight->Fill(mass); 
	massforLimitEtaTight->Fill(mass);
	massforLimitOmegaTight->Fill(mass);
	massforLimitPhiTight->Fill(mass);
	massforLimitJPsiPsiTight->Fill(mass);
	massforLimitUpsilonTight->Fill(mass);

	}
	
	}


	
	if(mass>60. && mass<120.){
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9 && 		      passIsoLoose){ forResolutionAMassZ->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.2 && passIsoLoose){ forResolutionBMassZ->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 && passIsoLoose){ forResolutionCMassZ->Fill(mass); }

	}


	if(mass>2.6 && mass<3.5){
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9 && 		      passIsoLoose){ forResolutionAMassJPsi->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.2 && passIsoLoose){ forResolutionBMassJPsi->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 && passIsoLoose){ forResolutionCMassJPsi->Fill(mass); }

	}


	float ma=0.2;
	for(int j=0; j<400.; j++){
      		ma = ma+(ma*0.01); 

     		if(mass > ma-(ma*0.01*10.) && mass < ma+(ma*0.01*10.)) {

		  	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<0.9 && 		      passIsoLoose){ massforLimit_CatA[j]->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=0.9 && maxEta<1.9 && passIsoLoose){ massforLimit_CatB[j]->Fill(mass); }
			if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta>=1.2 && maxEta<2.4 && passIsoLoose){ massforLimit_CatC[j]->Fill(mass); }

		}
      
   	 }	
	
      goodmuons.clear();
      }
      delete reader1;
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
