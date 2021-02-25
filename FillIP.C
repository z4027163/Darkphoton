//#include "sumwgt.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include <TTreeReader.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TTreeReaderValue.h>
#include "TLorentzVector.h"

void FillIP(string treepath = "scout_1.root", const char* outfilename = "./IP.root", bool isMC=false) {

 
   vector<string> filesToRun;
   vector<bool> isData;
//   string dirIn="/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/2018/ScoutingRunC/ScoutingCaloMuon/crab_20200617_174734/200617_154740/0000/";
   string dirIn="";
   filesToRun.push_back(dirIn.append(treepath));isData.push_back(true);
   TFile* outfile = new TFile(outfilename, "RECREATE");
   TTree* outtree = new TTree("tree", "tree");
     


    
    TH2F* BS = new TH2F("BS","BS",50,-1.,1.,50,-1.,1.);

    TH1F* IPZ = new TH1F("IPZ","IPZ",50,0,10.);
    TH1F* IPup = new TH1F("IPup","IPup",50,0,10.);
    TH1F* PVdZ = new TH1F("PVdZ","PVdZ",50,0,0.2);
    TH1F* PVdup = new TH1F("PVdup","PVdup",50,0,0.2);
    TH1F* VtxerrZ = new TH1F("VtxErrZ","VtxErrZ",50,0,0.05);
    TH1F* Vtxerrup = new TH1F("VtxErrup","VtxErrup",50,0,0.05);
   for(unsigned int f=0; f<filesToRun.size(); f++){


     cout<<filesToRun[f]<<endl;
    TChain* chain = new TChain("mmtree/tree");
    chain->Add((TString)filesToRun[f]);
    
    TTreeReader reader(chain);
    
    TTreeReaderValue<double>*                      wgt;
    TTreeReaderValue<double>*                      xsec;
    
    
    TTreeReaderValue<std::vector<float> >          mpt      (reader, "muonpt"     );
    TTreeReaderValue<std::vector<float> >          meta     (reader, "muoneta"    );
    TTreeReaderValue<std::vector<float> >          mphi     (reader, "muonphi"    );
    TTreeReaderValue<std::vector<float> >          mcharge  (reader, "muoncharge" );
    TTreeReaderValue<std::vector<char>  >          mid      (reader, "muonid"     );
    TTreeReaderValue<std::vector<float> >          chi2     (reader, "chi2"       );
    TTreeReaderValue<std::vector<float> >          dxy      (reader, "dxy"        );
    TTreeReaderValue<std::vector<float> >          dz       (reader, "dz"         );
    TTreeReaderValue<std::vector<float> >          cpiso    (reader, "cpiso"      );
    TTreeReaderValue<std::vector<float> >          chiso    (reader, "chiso"      );
    TTreeReaderValue<std::vector<float> >          phiso    (reader, "phiso"      );
    TTreeReaderValue<std::vector<float> >          nhiso    (reader, "nhiso"      );
    TTreeReaderValue<std::vector<float> >          puiso    (reader, "puiso"      );
    TTreeReaderValue<std::vector<float> >          tkiso    (reader, "tkiso"      );
    TTreeReaderValue<std::vector<bool> >           l1Result (reader, "l1Result"   );
    TTreeReaderValue<std::vector<unsigned char> >  nmhits   (reader, "nMuonHits"  );
    TTreeReaderValue<std::vector<unsigned char> >  nphits   (reader, "nPixelHits" );
    TTreeReaderValue<std::vector<unsigned char> >  ntklayers(reader, "nTkLayers"  );
    TTreeReaderValue<unsigned char>                hlt      (reader, "trig"       );
    TTreeReaderValue<unsigned>                     nverts   (reader, "nvtx"       );
    cout<<"here4"<<endl;
    TTreeReaderValue<std::vector<float> >          vtxX     (reader, "vtxX"       );
    TTreeReaderValue<std::vector<float> >          vtxY     (reader, "vtxY"       );
    TTreeReaderValue<std::vector<float> >          vtxZ     (reader, "vtxZ"       );
    //TTreeReaderValue<std::vector<float> >          vtxchi2  (reader, "vtxchi2"    );
    TTreeReaderValue<std::vector<float> >          vtxXError    (reader, "vtxXError"       );
    TTreeReaderValue<std::vector<float> >          vtxYError    (reader, "vtxYError"       );
    TTreeReaderValue<std::vector<float> >          vtxZError    (reader, "vtxZError"       );
   
    TTreeReaderValue<double>                       rho      (reader, "rho"        );
    TTreeReaderValue<unsigned int>                 run      (reader, "run"        );
    TTreeReaderValue<unsigned int>                 lumSec   (reader, "lumSec"     );
  
  //TTreeReaderValue<std::vector<int> >		   pdgid    (reader, "pdgid" );
  //TTreeReaderValue<std::vector<int> > 	   motherid (reader, "motherid");
     
    
//   if (isMC || !isData[f]) {
//        wgt  = new TTreeReaderValue<double>(reader, "wgt" );
//        xsec = new TTreeReaderValue<double>(reader, "xsec");
	
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
    int p=0;
    int count[4]={0};
    int count0=0;
    while(reader.Next()) {
      count0++;
      //if(!isData[f]) {if(pdgid->size()<2 || motherid->size()<2) continue;}
      if (((*hlt) & 2) == 0) continue;   
      if(!( l1Result->at(2) ||l1Result->at(7) || l1Result->at(11) || l1Result->at(6))) continue;
      count[0]++;
      bool passIso=false;
      bool passIsoLoose=false;
        double ea = (isMC ? 0.3 : 0.45);
        std::vector<unsigned> goodmuons;
      double PVd=1;
        float BSx = 0.086;
        float BSy = -0.033;
      if(treepath.find("_13TeV") != string::npos){ BSx = -0.029; BSy = 0.07;}
      if(isMC){BSx=0.0107; BSy=0.0417;}
      nvtx = *nverts;
      if(nvtx>0. ){
        PVd= sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) );
      }
      if(PVd>0.2) continue;
      if(mpt->size()<2) continue;
      count[1]++;
        for (std::size_t i = 0; i < mpt->size(); i++) {
	  //if ((*nmhits)[i] == 0)     continue;
          if((*mpt)[i]<4) continue;
          if(abs((*meta)[i])>1.9) continue;
	  if ((*nphits)[i] == 0)     continue;
          if ((*ntklayers)[i] <= 5)  continue;
          if ((*chi2)[i] > 10.)      continue;
          if((*tkiso)[i]>0.15) continue;
            goodmuons.push_back(i);
        }

        if (goodmuons.size() < 2) continue;
        unsigned idx1 = goodmuons[0];
        unsigned idx2 = goodmuons[1];
    

	if((*tkiso)[goodmuons[0]]<0.02 && (*tkiso)[goodmuons[1]]<0.02 ) passIso=true;
	if((*tkiso)[goodmuons[0]]<0.15 && (*tkiso)[goodmuons[1]]<0.15 ) passIsoLoose=true;

        if(passIsoLoose) count[2]++;
        if ((*mpt)[goodmuons[0]] < (*mpt)[goodmuons[1]]) {
            idx1 = goodmuons[1];
            idx2 = goodmuons[0];
	    }
	
        TLorentzVector mm;
	
        TLorentzVector m1;
        TLorentzVector m2;
	
        
        m1.SetPtEtaPhiM((*mpt)[idx1], (*meta)[idx1], (*mphi)[idx1], 0.1057);
        m2.SetPtEtaPhiM((*mpt)[idx2], (*meta)[idx2], (*mphi)[idx2], 0.1057);
	//m3.SetPtEtaPhiM((*mpt)[idx3], (*meta)[idx3], (*mphi)[idx3], 0.1057);
        //m4.SetPtEtaPhiM((*mpt)[idx4], (*meta)[idx4], (*mphi)[idx4], 0.1057);

//        mm += m1;
//        mm += m2;
        mm=m1+m2;
        m1pt   = m1.Pt();
        m1eta  = m1.Eta();
        m1phi  = m1.Phi();
        m1iso  = (*cpiso)[idx1] + (*phiso)[idx1] + (*nhiso)[idx1] - ea*(*rho);
        m1id   = (*mid)[idx1];
        m1ch   = (*mcharge)[idx1];

        m2pt   = m2.Pt();
        m2eta  = m2.Eta();
        m2phi  = m2.Phi();
        m2iso  = (*cpiso)[idx2] + (*phiso)[idx2] + (*nhiso)[idx2] - ea*(*rho);
        m2id   = (*mid)[idx2];
	m2ch   = (*mcharge)[idx2];
	
	Run=*run;
	LumSec=*lumSec;

        double weight=1;
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

	
        double IP=-1;	
        double vtxErr=1;
	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.2 ) passPVconstraint = true;
          double temx=(*vtxX)[0] - BSx;
          double temy=(*vtxY)[0] - BSy;
          double temr=sqrt(temx*temx+temy*temy);
          vtxErr=sqrt(temx*temx*(*vtxXError)[0]*(*vtxXError)[0]+temy*temy*(*vtxYError)[0]*(*vtxYError)[0])/temr;
          IP=PVd/vtxErr;
	  }

	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.1 ) passPVconstraintTight = true;
	  }

	if(m1ch*m2ch<0. && passPVconstraint && m1pt>slidePt1 && m2pt>slidePt2 && maxEta<1.9 && passIsoLoose){ 
             if(mass>86 && mass<96) {IPZ->Fill(IP); PVdZ->Fill(PVd); VtxerrZ->Fill(vtxErr);}
             if(mass>9.2 && mass<9.7) {IPup->Fill(IP); PVdup->Fill(PVd);Vtxerrup->Fill(vtxErr);}
	/*
	0.49 - 062 : eta
	0.68 - 0.88 : omega
	0.90 - 1.14 : phi
	2.7 - 4.2 : jpsi + psi(2s)
	8.5 - 11 : Upsilons
	*/
	
	}


	
      }
   }

   IPZ->Write();
   IPup->Write();
   PVdZ->Write();
   PVdup->Write();
   VtxerrZ->Write();
   Vtxerrup->Write();
    outfile->Close();

}
