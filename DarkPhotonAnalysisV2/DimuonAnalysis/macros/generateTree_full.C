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

void generateTree_full(string treepath = "scout_2.root", const char* outfilename = "./scoutData4MuonsM70.root", int era=1,  bool isMC=false) {
 
   vector<string> filesToRun;
   vector<bool> isData;
   string dirIn="/eos/cms/store/group/phys_exotica/darkPhoton/jakob/newProd/";
//   string dirIn="";
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
  
  //TTreeReaderValue<std::vector<int> >		   pdgid    (reader, "pdgid" );
  //TTreeReaderValue<std::vector<int> > 	   motherid (reader, "motherid");
     
/*    
   if (isMC || !isData[f]) {
        wgt  = new TTreeReaderValue<double>(reader, "wgt" );
        xsec = new TTreeReaderValue<double>(reader, "xsec");
	
    }
*/	 
 

    double wgtsum  = 1.0;//isMC ? sumwgt(treepath) : 1.0;
    int nmuon =0;
    float  puwgt   = 1.0;
    float  mcwgt   = 1.0;
    float charge[4]={0};
    float  mupt[4]    = {0};        
    float  mueta[4]   = {0};        
    float  muphi[4]   = {0};        
    float  mass    = 0.0;        
    int   muid[4]    = {0};
    float  much[4]    = {0.}; 
    int  muphits[4]    = {0}; 
    int  mumhits[4] = {0};
    int  munlayer[4]    = {0}; 
    float  muchi[4] = {0.};
    float  mutrkiso[4] = {0};
    vector<bool> l1Bools;
    unsigned nvtx  = 0;
    unsigned Run   = 0;
    unsigned LumSec   = 0;
    float Dxy[4]={-100,-100,-100,-100};
    float Dz[4]={-100,-100,-100,-100};
    float vchi2 = 0;
    float PVd, IP,vtxErr;
    int l1trg;
    float vtx_z;
 
    if (isMC) {
    outtree->Branch("mcwgt" , &mcwgt , "mcwgt/F" );
    outtree->Branch("puwgt" , &puwgt , "puwgt/F" );
    }
    outtree->Branch("pt"  , mupt  , "pt[4]/F"  );
    outtree->Branch("eta" , mueta , "eta[4]/F" );
    outtree->Branch("phi" , muphi , "phi[4]/F" );
    outtree->Branch("trkiso",mutrkiso,"trkiso[4]/F");
    outtree->Branch("nphits",muphits,"m2phits[4]/I");
    outtree->Branch("nmhits",mumhits,"m2mhits[4]/I");
    outtree->Branch("muid",muid,"id[4]/I");
    outtree->Branch("chi",muchi,"chi[4]/F");
    outtree->Branch("charge",charge,"charge[4]/F");
    outtree->Branch("ntklayers",munlayer,"ntklayers[4]/I");
    outtree->Branch("nvtx"  , &nvtx  , "nvtx/i"  );
    outtree->Branch("Run"   , &Run   , "Run/i"  );
    outtree->Branch("LumSec", &LumSec, "LumSec/i"  );
    outtree->Branch("nmuon", &nmuon, "nmuon/I"  );
    outtree->Branch("dxy",Dxy,"dxy[4]/F");
    outtree->Branch("dz",Dz,"dz[4]/F");
    outtree->Branch("vtxchi2",&vchi2,"vtxchie2/F");
    outtree->Branch("PVd",&PVd,"PVd/F");
    outtree->Branch("IP",&IP,"IP/F");
    outtree->Branch("vtxErr",&vtxErr,"vtxErr/F");
    outtree->Branch("vtxZ",&vtx_z,"vtxZ/F");
    outtree->Branch("l1trg"  , &l1trg , "l1trg/i" );
    int p=0;
    while(reader.Next()) {
      //if(!isData[f]) {if(pdgid->size()<2 || motherid->size()<2) continue;}
      if (((*hlt) & 2) == 0) continue;   
      bool passIso=false;
      bool passIsoLoose=false;
        double ea = (isMC ? 0.3 : 0.45);
        std::vector<unsigned> goodmuons;

        int idx=0;

        float BSx = 0.092;
        float BSy = -0.06;
        float phi0 = -0.578;
        if(era==1) {BSx=0.0845; BSy=-0.0335; phi0=-0.3774;}
        if(era==2) {BSx=0.0835; BSy=-0.032; phi0=-0.366;}
        if(era==3) {BSx=0.0845; BSy=-0.032; phi0=-0.362;}
        if(era==4) {BSx=0.0835; BSy=-0.0295;phi0=-0.3396;}
        if(era==5) {BSx=0.0965; BSy=-0.066; phi0=-0.600;}
        if(era==6) {BSx=0.097; BSy=-0.064; phi0=-0.583;}
        if(era==7) {BSx=0.096; BSy=-0.0645; phi0=-0.5868;}
        if(era==8) {BSx=0.0965; BSy=-0.062; phi0=-0.571;}   //phi=arctan(BSy/BSx)

        if(isMC){ BSx=0.0107; BSy=0.0417; phi0=1.320;}
        float r0=sqrt(BSx*BSx+BSy*BSy);
       
        for(int i=0; i<4;i++){
          muid[i]=0;
        }
        for (std::size_t i = 0; i < mpt->size(); i++) {
           if(idx>=4) continue;
           double pt=(*mpt)[i];
           if((*mpt)[i]<3 || fabs((*meta)[i])>1.9) continue;
           mupt[idx]=(*mpt)[i];
           mueta[idx]=(*meta)[i];
           muphi[idx]=(*mphi)[i];
           mutrkiso[idx]=(*tkiso)[i];
           muphits[idx]=(*nphits)[i];
           mumhits[idx]=(*nmhits)[i];
           muchi[idx]=(*chi2)[i];
           munlayer[idx]=(*ntklayers)[i];
           Dxy[idx]=(*dxy)[i]+r0*sin((*mphi)[i]-phi0);
           Dz[idx]=(*dz)[i];
           muid[idx]=0;
           charge[idx]=(*mcharge)[i];
           if((*nphits)[i] > 0 && (*ntklayers)[i] > 5 && (*chi2)[i] < 10. && (*tkiso)[i]<0.15) muid[idx]=1;
           idx++;
        }
        
        nmuon=idx;
        if(nmuon<2) continue;

//        if(nmuon==2 && mupt[2]>0) cout  << "nm=" << idx  << " 4 muon id" <<muid[0]<<"," <<muid[1]<<"," << mupt[2]<<"," <<muid[3] << endl;
        count2++;
      //only check first muon for tag
        unsigned idx1 = 0;
        unsigned idx2 = 1;
        int ntag=0;
	//	mass4  = mmmm.M();
        nvtx = *nverts;
        
	
        bool passPVconstraint = false;
	bool passPVconstraintTight = false;
          
	if(nvtx>0. ){
	  if(  sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) )  < 0.2 ) passPVconstraint = true;
          vchi2 = (*vtxchi2)[0];
          PVd=sqrt( ((*vtxX)[0] - BSx)*((*vtxX)[0] - BSx) + ((*vtxY)[0] - BSy)*((*vtxY)[0] - BSy) );
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
        if(nvtx<=0.) continue;
         if(era<=4){if(!(l1Result->at(4) || l1Result->at(11) || l1Result->at(12) || l1Result->at(10))) continue;}
         else{ if(!(l1Result->at(2) ||l1Result->at(7) || l1Result->at(11) || l1Result->at(6))) continue;}

         vtx_z=(*vtxZ)[0];
         
         l1trg=0;
         if(era<=4){
           if(l1Result->at(4)) l1trg+=1;
           if(l1Result->at(10)) l1trg+=2;
           if(l1Result->at(11)) l1trg+=4;
           if(l1Result->at(12)) l1trg+=8;
         }
         else{ 
           if(l1Result->at(2)) l1trg+=1;
           if(l1Result->at(6)) l1trg+=2;
           if(l1Result->at(7)) l1trg+=4;
           if(l1Result->at(11)) l1trg+=8;
         }        

	 count++;     	
         outtree->Fill();	
      }
   }
   cout << "2m=" << count2 << "\n final=" << count << endl;
   tof->cd();
   outtree->Write("fitter_tree",TObject::kOverwrite);
   tof->Close();
    outfile->Close();
}
