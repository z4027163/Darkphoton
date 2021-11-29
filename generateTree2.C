//#include "sumwgt.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include <TTreeReader.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TTreeReaderValue.h>
#include "TLorentzVector.h"
const double mPI = 3.141592654;

double DELTAPHI( double phi1, double phi2 ){
        if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
          return -999;
        }
        float dp=fabs(phi1-phi2);
        if (dp>mPI) dp-=float(2*mPI);
        return dp;
}

void generateTree2(string treepath = "scout_1.root", const char* outfilename = "./scoutData4MuonsM70.root", bool isMC=false) {
 
   vector<string> filesToRun;
   vector<bool> isData;
   string dirIn="/eos/user/w/wangz/ntuple_darkphoton/";
   //string dirIn="";
   filesToRun.push_back(dirIn.append(treepath));isData.push_back(true);
   TFile* outfile = new TFile(outfilename, "RECREATE");

    outfile->cd();
    TDirectory *tof = outfile->mkdir("tpTree");

   TTree* outtree = new TTree("fitter_tree", "fitter_tree");
     
   int count=0;
   int count2=0;
   int count3=0;
   int count4=0;
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
    TTreeReaderValue<std::vector<float> >          vtxchi2  (reader, "vtxChi2"    );
    TTreeReaderValue<std::vector<float> >          vtxXError    (reader, "vtxXError"       );
    TTreeReaderValue<std::vector<float> >          vtxYError    (reader, "vtxYError"       );
    TTreeReaderValue<std::vector<float> >          vtxZError    (reader, "vtxZError"       );
   
    TTreeReaderValue<double>                       rho      (reader, "rho"        );
    TTreeReaderValue<unsigned int>                 run      (reader, "run"        );
    TTreeReaderValue<unsigned int>                 lumSec   (reader, "lumSec"     );
    TTreeReaderValue<unsigned int>                 putrue   (reader, "putrue"     ); 
  //TTreeReaderValue<std::vector<int> >		   pdgid    (reader, "pdgid" );
  //TTreeReaderValue<std::vector<int> > 	   motherid (reader, "motherid");
     
/*    
   if (isMC || !isData[f]) {
        wgt  = new TTreeReaderValue<double>(reader, "wgt" );
        xsec = new TTreeReaderValue<double>(reader, "xsec");
	
    }
*/	 
 

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
    int   m1id    = 0;
    int   m2id    = 0;
    char   m3id    = 0;
    char   m4id    = 0;
    float  m1ch    = 0.; 
    float  m2ch    = 0.;
    int  m2phits    = 0.; 
    int  m2mhits = 0.;
    int  m2layer    = 0.; 
    float  m2chi = 0.;
    float  m2trkiso = 0;
    float dR=-1;
    vector<bool> l1Bools;
    unsigned nvtx  = 0;
    unsigned Run   = 0;
    unsigned LumSec   = 0;
    unsigned npu = 0;
    float Dxy=0;
    float Dz=0;
    float vchi2=0;
    float PVd=0;
    float IP=0;
    float vtxErr=0;

    if (isMC) {
    outtree->Branch("mcwgt" , &mcwgt , "mcwgt/F" );
    outtree->Branch("puwgt" , &puwgt , "puwgt/F" );
    }
    outtree->Branch("m1pt"  , &m1pt  , "m1pt/F"  );
    outtree->Branch("m1eta" , &m1eta , "m1eta/F" );
    outtree->Branch("m1phi" , &m1phi , "m1phi/F" );
    outtree->Branch("m1id"  , &m1id  , "m1id/I"  );
    outtree->Branch("pt"  , &m2pt  , "m2pt/F"  );
    outtree->Branch("eta" , &m2eta , "m2eta/F" );
    outtree->Branch("phi" , &m2phi , "m2phi/F" );
    outtree->Branch("trkiso",&m2trkiso,"trkiso/F");
    outtree->Branch("nphits",&m2phits,"m2phits/I");
    outtree->Branch("nmhits",&m2mhits,"m2mhits/I");
    outtree->Branch("chi",&m2chi,"m2chi/F");
    outtree->Branch("ntklayers",&m2layer,"ntklayers/I");
    outtree->Branch("m2id"  , &m2id  , "m2id/I"  );
    outtree->Branch("mass"  , &mass  , "mass/F"  );
    outtree->Branch("nvtx"  , &nvtx  , "nvtx/i"  );
    outtree->Branch("Run"   , &Run   , "Run/i"  );
    outtree->Branch("npu"   , &npu   , "npu/i"  );
    outtree->Branch("LumSec", &LumSec, "LumSec/i"  );
    outtree->Branch("dR", &dR,"dR/F");
    outtree->Branch("dxy",&Dxy,"dxy/F");
    outtree->Branch("dz",&Dz,"dz/F");
    outtree->Branch("vtxchi2",&vchi2,"vtxchi2/F");
    outtree->Branch("PVd",&PVd,"PVd/F");
    outtree->Branch("IP",&IP,"IP/F");
    outtree->Branch("vtxErr",&vtxErr,"vtxErr/F");
    //outtree->Branch("l1Bools"  , &l1Bools , "l1Bools" );
    int p=0;
    while(reader.Next()) {
      //if(!isData[f]) {if(pdgid->size()<2 || motherid->size()<2) continue;}
      if (((*hlt) & 2) == 0) continue;   
      bool passIso=false;
      bool passIsoLoose=false;
        double ea = (isMC ? 0.3 : 0.45);
        std::vector<unsigned> goodmuons;
/*
        for (std::size_t i = 0; i < mpt->size(); i++) {
	  if ((*nphits)[i] == 0)     continue;
          if ((*ntklayers)[i] <= 5)  continue;
          if ((*chi2)[i] > 10.)      continue;
            double iso = (*cpiso)[i] + (*nhiso)[i] + (*phiso)[i] - ea*(*rho);
	   

            goodmuons.push_back(i);
        }
        if (goodmuons.size() < 2) continue;
*/
        if(mpt->size()<2) continue;

        count2++;
      //only check first muon for tag
        if(((*mpt)[0])<12 || fabs((*meta)[0])>2.4) continue;
        if(((*nphits)[0] == 0) || (*ntklayers)[0] <= 5 || (*chi2)[0] > 10. || (*tkiso)[0]>0.15) continue; 

        count3++;
        unsigned idx1 = 0;
        unsigned idx2 = 1;
        int nprob=0;
        for(std::size_t i = 1; i < mpt->size(); i++) {
             if(((*mpt)[i])<3 || fabs((*meta)[i])>1.9) continue;
             if((*mcharge)[0]*(*mcharge)[i]>0) continue;
             nprob++;
             idx2=i;
        }       
        if(nprob!=1) continue;  
        TLorentzVector mm;
	
        TLorentzVector m1;
        TLorentzVector m2;
	
        count4++;
        m1.SetPtEtaPhiM((*mpt)[idx1], (*meta)[idx1], (*mphi)[idx1], 0.1057);
        m2.SetPtEtaPhiM((*mpt)[idx2], (*meta)[idx2], (*mphi)[idx2], 0.1057);
	//m3.SetPtEtaPhiM((*mpt)[idx3], (*meta)[idx3], (*mphi)[idx3], 0.1057);
        //m4.SetPtEtaPhiM((*mpt)[idx4], (*meta)[idx4], (*mphi)[idx4], 0.1057);
        
         
        mm += m1;
        mm += m2;
        m1pt   = m1.Pt();
        m1eta  = m1.Eta();
        m1phi  = m1.Phi();
        m1iso  = (*tkiso)[idx1];
        m1id=1;
        m1ch   = (*mcharge)[idx1];

        m2pt   = m2.Pt();
        m2eta  = m2.Eta();
        m2phi  = m2.Phi();
//        m2iso  = (*cpiso)[idx2] + (*phiso)[idx2] + (*nhiso)[idx2] - ea*(*rho);
//        cout << "cpiso=" << (*cpiso)[idx2] << " phiso=" << (*phiso)[idx2] << " nhiso=" << (*nhiso)[idx2] << " rho=" << *rho << endl;
        m2id=0;
        m2trkiso = (*tkiso)[idx2];
        m2phits=(*nphits)[idx2];
        m2mhits=(*nmhits)[idx2];
//        cout << "iso=" << (*tkiso)[idx2] << endl;
        m2layer=(*ntklayers)[idx2];
        m2chi=(*chi2)[idx2];
        if((*tkiso)[idx2]<0.15 && ((*nphits)[idx2] > 0) && (*ntklayers)[idx2] > 5 && (*chi2)[idx2] < 10.) m2id=1;
	m2ch   = (*mcharge)[idx2];
	
	Run=*run;
	LumSec=*lumSec;
        npu=*putrue;

        Dxy=(*dxy)[idx2];
        Dz=(*dz)[idx2]; 
        dR=sqrt((m1eta-m2eta)*(m1eta-m2eta)+DELTAPHI(m1phi,m2phi)*DELTAPHI(m1phi,m2phi));
//        if(m1ch*m2ch>0) continue;

	//'L1_DoubleMu4_SQ_OS_dR_Max1p2', [12]
	//'L1_DoubleMu4p5_SQ_OS_dR_Max1p2', [16]
	//
        mass   = mm.M();
	//	mass4  = mmmm.M();
        nvtx = *nverts;
        
	
        bool passPVconstraint = false;
	bool passPVconstraintTight = false;
	float BSx = 0.092;
	float BSy = -0.06;
//	if(treepath.find("_13TeV") != string::npos){ BSx = -0.029; BSy = 0.07;}
        if(isMC){ BSx=-0.02446; BSy=0.06924;} //upsilon sample
//        if(isMC){ BSx=-0.02465; BSy=0.06962;}   //jpsi sample
	float slidePt1 = mm.M()/3.;
	if(slidePt1<4.) slidePt1 = 4.;
	float slidePt2 = mm.M()/4.;
	if(slidePt2<4.) slidePt2 = 4.;

	float maxEta=TMath::Max(abs(m1eta),abs(m2eta));

	
	

	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.2 ) passPVconstraint = true;
          vchi2=(*vtxchi2)[0];
          PVd=sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) );
          // IP=r/(sigma_r*cos(theta))
          double temx=(*vtxX)[0] - BSx;
          double temy=(*vtxY)[0] - BSy;
          double temr=sqrt(temx*temx+temy*temy);
          vtxErr=sqrt(temx*temx*(*vtxXError)[0]*(*vtxXError)[0]+temy*temy*(*vtxYError)[0]*(*vtxYError)[0])/temr;
          IP=PVd/vtxErr;
	  }

	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.1 ) passPVconstraintTight = true;
	  }
//        if(!passPVconstraint) continue;
//        if(!( l1Result->at(2) || l1Result->at(6) || l1Result->at(7) || l1Result->at(10) || l1Result->at(11) || l1Result->at(13))) continue;
//        if(!( l1Result->at(15))) continue;
        if(nvtx<=0) continue;
 
	 count++;     	
         outtree->Fill();	
      }
   }
   cout << "2m=" << count2 << "\n 2m+tag=" << count3 << "\n 2m+tag+probe=" << count4 << endl;
   cout << "final count=" << count << endl;
   tof->cd();
   outtree->Write("fitter_tree",TObject::kOverwrite);
   tof->Close();
    outfile->Close();
}
