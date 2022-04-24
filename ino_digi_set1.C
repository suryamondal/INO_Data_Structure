

// #define isDebug


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <bitset>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"

#include "CauTree.h"
#include "RPCEve.h"

#include "RunInfo.h"
#include "DigiStore.h"
#include "CauStore.h"
#include "HealthStore.h"


using namespace std;

/*
  

*/

const int nside = 2;
const int nlayer = 10;
const int nstrip = 64;
const int nTDC = 8;
const int CAUselectline = 2;
const int maxCAUhit = 1;


void getRPCId(UShort_t  rpcId,
	      Int_t  &module,
	      Int_t  &xrow,
	      Int_t  &yrow,
	      Int_t  &zlay) {
  // rpcId : 2 bit for INO module  
  //         3 bit for X-row
  //         3 bit for Y row
  //         8 bit for Z-layer
  zlay   = (rpcId    )&0xFF;//0b11111111;
  yrow   = (rpcId>> 8)&0b111;
  xrow   = (rpcId>>11)&0b111;
  module = (rpcId>>14)&0b11;
};

UShort_t constructRPCId(Int_t module,
			Int_t xrow,
			Int_t yrow,
			Int_t zlay) {
  UShort_t rpcId = module;
  rpcId        <<= 3;
  rpcId         +=  xrow;
  rpcId        <<= 3;
  rpcId         +=  yrow;
  rpcId        <<= 8;
  rpcId         +=  zlay;
  return rpcId;
};

UInt_t formTDC(Int_t  tdcval,
	       Bool_t isTrail,	// leading or trailing
	       Bool_t side,	// x or y
	       Int_t  tdcno) {	// tdc no 0-7
  UInt_t tdcout = tdcval;
  tdcout      <<= 1;
  tdcout       += isTrail;
  tdcout      <<= 4;
  tdcout       += (tdcno + (side?nTDC:0));
  return tdcout;
};

UInt_t formCAUmap(UShort_t rpcId,
		  UShort_t cauNo) {
  /* returns <CAU channel no 16bit> <rpc Id 16bit>  */
  UInt_t mapitem = cauNo;
  mapitem      <<= 16;
  mapitem       += rpcId;
  return mapitem;
};

ULong64_t formCAU(Int_t    cauTDC,
		  Bool_t   isTrail, // leading or trailing
		  UShort_t rpcId) {
  ULong64_t tdcout = cauTDC;
  tdcout         <<= 1;
  tdcout          += isTrail;
  tdcout         <<= 16;
  tdcout          += rpcId;
  return tdcout;
};




int main(int argc, char** argv) {

  RunInfo *runinfo = new RunInfo();
  runinfo->TriggerInfo = "5 out of 10";
  runinfo->Creator = "User";
  runinfo->CreatedOn.Set();
#ifdef isDebug
  cout << " time " << runinfo->CreatedOn << endl;
#endif	// #ifdef isDebug

  
  vector<UInt_t> CAUmap[CAUselectline];
  for(int ij=0;ij<CAUselectline;ij++) {CAUmap[ij].clear();}
  
  /* Entry <CAU channel no 16bit> <rpc Id 16bit>  */
  /* layer map for selectline=0 */
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 0),  0));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 1),  2));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 2),  4));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 3),  6));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 4),  8));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 5), 10));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 6), 12));
  CAUmap[0].push_back(formCAUmap(constructRPCId(0, 0, 0, 7), 14));
  /* layer map for selectline=1 */
  CAUmap[1].push_back(formCAUmap(constructRPCId(0, 0, 0, 8),  0));
  CAUmap[1].push_back(formCAUmap(constructRPCId(0, 0, 0, 9),  1));
  
#ifdef isDebug
    for(int ij=0;ij<CAUselectline;ij++) {
      for(int jk=0;jk<(CAUmap[ij].size());jk++) {
	cout << " " << ij << " " << jk
	     << " " << bitset<32>(CAUmap[ij][jk]) << endl;
      }
    }
#endif  // #ifdef isDebug
  
  
  Long64_t start_s = clock();
  
  char outfil[300] = {}; 
  char outfilx[300] = {};
  int len = strlen(argv[1]);
  strncpy(outfil, argv[1], len-4); // rre file
  // sprintf(outfilx, "/media/jim/INO2_mical_SSD/miniICAL_magnetic_data_GMA/GMA_%s.root", outfil);
  sprintf(outfilx, "/media/jim/INO2_mical_SSD/temp_surya/Anal_Surya_Format/data/GMA_%s.root", outfil);
  // sprintf(outfilx, "../rredata/GMA_%s.root", outfil);
  TFile *f1=new TFile(outfilx,"RECREATE");
  
  // TFile* f1=new TFile("abc.root", "recreate");

  TTree* TEve =new TTree("RPCtree", "INO_digi");
  DigiStore* digiall = new DigiStore();
  TEve->Branch("EventData","DigiStore",&digiall, 16000,2);
  
  TTree* TCau =new TTree("CAUtree", "INO_digi");
  CauStore* cauall = new CauStore();
  TCau->Branch("CauData","CauStore",&cauall, 16000,2);
  
  TTree* THealth =new TTree("HEALTHtree", "INO_digi");
  HealthStore* healthall = new HealthStore();
  TCau->Branch("HealthData","HealthStore",&healthall, 16000,2);
  
  DigiOpt1* digi1;
  // CauOpt1* cau1;
  vector<UInt_t> tdc;
  // vector<ULong64_t> cau_tdc, cau_tdc_ref1;
  
  char datafile[300] = {};
  strncpy(datafile,argv[1],300);
  Long64_t nentrymx = stoi(argv[3]);
  Long64_t nentrymn = stoi(argv[2]);
  char infile[300] = {};
  // sprintf(infile, "/media/jim/INO2_mical_SSD/miniICAL_magnetic_data/%s", datafile);
  sprintf(infile, "/media/jim/INO2_mical_SSD/temp_surya/Anal_Surya_Format/data/%s", datafile);
  // sprintf(infile, "../rredata/%s", datafile);
  TFile *fileIn = new TFile(infile, "read");
  
  if(!fileIn->IsZombie()) {
    
    TTree *event_tree=(TTree*)fileIn->Get("evetree");
    RPCEve *event = new RPCEve(event_tree);
    TTree *cau_tree=(TTree*)fileIn->Get("cautree");
    CauTree *cau = new CauTree(cau_tree);
    cau->Loop();
    event->Loop();
  
    Long64_t Cnentry = cau_tree->GetEntries();
    cout << " cau entry " << Cnentry << endl;
    for(Long64_t iev=0;iev<Cnentry;iev++) {
    
      if(iev%1000==0) {
	Long64_t stop_s = clock();
	cout << " cau " << iev
	     << " time " << (stop_s-start_s)/double(CLOCKS_PER_SEC)
	     << endl;
      }
      
      cau_tree->GetEntry(iev);
    
      // cauall->ClearCau();
      if((cau->select_line)>=CAUselectline) {continue;}
      
      cauall->CauNum = cau->CauNum;
      cauall->Cautime = *cau->Cautime;
      cauall->select_line = cau->select_line;
      cauall->raw_trig_cnt = cau->raw_trig_cnt;
      cauall->acpt_trig_cnt = cau->acpt_trig_cnt;
      cauall->cau_status = cau->cau_status;
      cauall->cau_ref_l[0] = cau->cau_ref1_l;
      cauall->cau_ref_l[1] = cau->cau_ref2_l;
      cauall->cau_ref_t[0] = cau->cau_ref1_t;
      cauall->cau_ref_t[1] = cau->cau_ref2_t;
      
      cauall->cau_tdc.clear();
      cauall->cau_tdc_ref1.clear();
      // ULong64_t tdcId1;
      
      for(int nt=0;nt<int(CAUmap[cau->select_line].size());nt++) {
	
	UShort_t rpcId1 = CAUmap[cau->select_line][nt]&0xFFFF;
	UShort_t cauNo1 = CAUmap[cau->select_line][nt]>>16;
	// cout << " " << rpcId1 << " " << cauNo1 << endl;
	
	for(int ix=0;
      	    ix<TMath::Min(maxCAUhit,int(cau->cau_tdc_l[cauNo1]->size()));
      	    ix++) {
	  // tdcId1 = cau->cau_tdc_l[cauNo1]->at(ix);
	  // tdcId1 <<= 1;
	  // tdcId1 += 0;		// leading
	  // tdcId1 <<= 16;
	  // tdcId1 += rpcId1;
	  // cauall->cau_tdc.push_back(tdcId1);
	  cauall->cau_tdc.push_back(formCAU(cau->cau_tdc_l[cauNo1]->at(ix),0,rpcId1));
	}
	for(int ix=0;
      	    ix<TMath::Min(maxCAUhit,int(cau->cau_tdc_t[cauNo1]->size()));
      	    ix++) {
	  // tdcId1 = cau->cau_tdc_t[cauNo1]->at(0);
	  // tdcId1 <<= 1;
	  // tdcId1 += 1;		// trailing
	  // tdcId1 <<= 16;
	  // tdcId1 += rpcId1;
	  // cauall->cau_tdc.push_back(tdcId1);
	  cauall->cau_tdc.push_back(formCAU(cau->cau_tdc_t[cauNo1]->at(ix),1,rpcId1));
	}
	
	for(int ix=0;
      	    ix<TMath::Min(maxCAUhit,int(cau->cau_tdc_ref1_l[cauNo1]->size()));
      	    ix++) {
	  // cout << " something " << rpcId1 << " " << ix << endl;
	  // tdcId1 = cau->cau_tdc_ref1_l[cauNo1]->at(ix);
	  // tdcId1 <<= 1;
	  // tdcId1 += 0;		// leading
	  // tdcId1 <<= 16;
	  // tdcId1 += rpcId1;
	  // cauall->cau_tdc_ref1.push_back(tdcId1);
	  cauall->cau_tdc_ref1.push_back(formCAU(cau->cau_tdc_ref1_l[cauNo1]->at(ix),0,rpcId1));
	}
	for(int ix=0;
      	    ix<TMath::Min(maxCAUhit,int(cau->cau_tdc_ref1_t[cauNo1]->size()));
      	    ix++) {
	  // tdcId1 = cau->cau_tdc_ref1_t[cauNo1]->at(ix);
	  // tdcId1 <<= 1;
	  // tdcId1 += 1;		// trailing
	  // tdcId1 <<= 16;
	  // tdcId1 += rpcId1;
	  // cauall->cau_tdc_ref1.push_back(tdcId1);
	  cauall->cau_tdc_ref1.push_back(formCAU(cau->cau_tdc_ref1_t[cauNo1]->at(ix),1,rpcId1));
	}
      }	// for(int nt=0;nt<int(CAUmap[cau->select_line].size());nt++) {
      
      TCau->Fill();
      
    } // for(Long64_t iev=0;iev<Cnentry;iev++) {
    
    
    Long64_t nentry = event_tree->GetEntries();
    nentry = min(nentry,nentrymx);
    cout << datafile <<" has "<<nentry<<" events "<<endl;
  
    for(Long64_t iev=nentrymn
	  ;iev<nentry;iev++) {

#ifdef isDebug
      cout << " " << iev << endl;
#endif  // #ifdef isDebug

      if(iev%100000==0) {
	Long64_t stop_s = clock();
	cout << " iev " << iev
	     << " time " << (stop_s-start_s)/double(CLOCKS_PER_SEC)
	     << endl;
      }
    
      fileIn->cd();
      event_tree->GetEntry(iev);
    
      digiall->ClearDigi();
    
      digiall->ENum  = event->ENum[0];
      digiall->REnum = event->REnum[0];
      digiall->CEnum = event->CEnum;
      digiall->EventTime = *event->EveTS[0];
      // cout << " time " << digiall->EventTime.AsDouble() << endl;
    
      if(iev == 0)        {runinfo->StartTime = *event->EveTS[0];}
      if(iev == nentry-1) {runinfo->StopTime  = *event->EveTS[0];}
    
    
      for(int ij=0;ij<nlayer;ij++) {
	
	UShort_t rpcId = constructRPCId(0,0,0,ij);
	ULong64_t xystrp[2] = {0};
	tdc.clear();
	
	// rpcId += 0;		// INO module
	// rpcId <<= 3;
	// rpcId += 0;		// X-row
	// rpcId <<= 3;
	// rpcId += 0;		// Y-row
	// rpcId <<= 8;
	// rpcId += ij;		// Z-layer
	
	for(int kl=nstrip-1; kl>=0; kl--) {
	  xystrp[0] <<= 1;
	  xystrp[0] += int(event->xLayer[ij]->TestBitNumber(kl));
	  xystrp[1] <<= 1;
	  xystrp[1] += int(event->yLayer[ij]->TestBitNumber(kl));
	} // for(int kl=nstrip-1; kl>=0; kl--) {
	
	tdc.push_back(UInt_t(event->tdc_ref_l[ij]));
	tdc.push_back(UInt_t(event->tdc_ref_t[ij]));
	
#ifdef isDebug
	cout << " tdc_ref " << event->tdc_ref_l[ij]
	     << " " << event->tdc_ref_t[ij] << endl;
	cout << " tdc_ref " << tdc[0] << " " << tdc[1] << endl;
#endif  // #ifdef isDebug
	
	/* tdc leading edge */
	/* X side */
	for(int nt=0;nt<nTDC;nt++) {
	  for(int ix=0;
	      ix<int(event->vxtdc_l[ij][nt]->size());
	      ix++) {
	    // UInt_t tdcId = event->vxtdc_l[ij][nt]->at(ix);
	    // tdcId <<= 1;
	    // tdcId += 0;		// leading
	    // tdcId <<= 4;
	    // tdcId += nt;		// x-side
	    // tdc.push_back(tdcId);
	    // cout << " " << tdcId << " " << formTDC(event->vxtdc_l[ij][nt]->at(ix),0,0,nt) << endl;
	    tdc.push_back(formTDC(event->vxtdc_l[ij][nt]->at(ix),0,0,nt));
	  }
	}	// for(int nt=0;nt<nTDC;nt++) {
	/* Y side */
	for(int nt=0;nt<nTDC;nt++) {
	  for(int ix=0;
	      ix<int(event->vytdc_l[ij][nt]->size());
	      ix++) {
	    // UInt_t tdcId = event->vytdc_l[ij][nt]->at(ix);
	    // tdcId <<= 1;
	    // tdcId += 0;		// leading
	    // tdcId <<= 4;
	    // tdcId += (nt + nTDC);		// x-side
	    // tdc.push_back(tdcId);
	    tdc.push_back(formTDC(event->vytdc_l[ij][nt]->at(ix),0,1,nt));
	  }
	}	// for(int nt=0;nt<nTDC;nt++) {

	/* tdc trailing edge */
	/* X side */
	for(int nt=0;nt<nTDC;nt++) {
	  for(int ix=0;
	      ix<int(event->vxtdc_t[ij][nt]->size());
	      ix++) {
	    // UInt_t tdcId = event->vxtdc_t[ij][nt]->at(ix);
	    // tdcId <<= 1;
	    // tdcId += 1;		// trailing
	    // tdcId <<= 4;
	    // tdcId += nt;		// x-side
	    // tdc.push_back(tdcId);
	    // cout << " " << tdcId << " " << formTDC(event->vxtdc_t[ij][nt]->at(ix),1,0,nt) << endl;
	    tdc.push_back(formTDC(event->vxtdc_t[ij][nt]->at(ix),1,0,nt));
	  }
	}	// for(int nt=0;nt<nTDC;nt++) {
	/* Y side */
	for(int nt=0;nt<nTDC;nt++) {
	  for(int ix=0;
	      ix<int(event->vytdc_t[ij][nt]->size());
	      ix++) {
	    // UInt_t tdcId = event->vytdc_t[ij][nt]->at(ix);
	    // tdcId <<= 1;
	    // tdcId += 1;		// trailing
	    // tdcId <<= 4;
	    // tdcId += (nt + nTDC);		// x-side
	    // tdc.push_back(tdcId);
	    // cout << " " << tdcId << " " << formTDC(event->vytdc_t[ij][nt]->at(ix),1,1,nt) << endl;
	    tdc.push_back(formTDC(event->vytdc_t[ij][nt]->at(ix),1,1,nt));
	  }
	}	// for(int nt=0;nt<nTDC;nt++) {
      
	if(xystrp[0]==0 && xystrp[1]==0
	   && int(tdc.size())<=2) {continue;}
	
#ifdef isDebug
	cout << " " << bitset<16>(rpcId) << endl
	     << " " << bitset<64>(xystrp[0])
	     << " " << bitset<64>(xystrp[1])
	     << " " << int(tdc.size())
	     << endl;
#endif  // #ifdef isDebug
	
	// cout << " test " << endl;
	digi1 = digiall->AddDigi1Store();
	// cout << " test1 " << endl;
	digi1->rpcId = rpcId;
	digi1->strips[0] = xystrp[0];
	digi1->strips[1] = xystrp[1];
	digi1->tdc.clear();
	// digi1->tdc = tdc;
	for(int pp=0;pp<int(tdc.size());pp++) {
	  digi1->tdc.push_back(tdc[pp]);
#ifdef isDebug
	  cout << " " << int(digi1->tdc[pp]);
#endif  // #ifdef isDebug
	}
#ifdef isDebug
	cout << endl;
#endif  // #ifdef isDebug
	
      }   // for(int ij=0;ij<nlayer;ij++) {
      
      // digiall->CompressDigi();
      
      TEve->Fill();
      
    } // for(Long64_t iev=0;iev<nentry;iev++) {
  
  } // if(!fileIn->IsZombie()) {
  
  f1->cd();
  runinfo->Write("run");
  TCau->Write();
  TEve->Write();
  // f1->Write();
  f1->Purge();
  f1->Close();
}

