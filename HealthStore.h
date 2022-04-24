#ifndef ROOT_HealthStore
#define ROOT_HealthStore

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include <vector>

#include "MonOpt1.h"

//#include <pair>
using namespace std;

/* 

 */


class HealthStore : public TObject {
 public:
  
  TTimeStamp               HealthTime;
  
  UInt_t NMonOpt1;		// Number of MonStore not sets!
  TClonesArray  *fMonOpt1;	//->array with all MonStore
  static TClonesArray *fgMonOpt1;
  
  // vector<UShort_t>         rpcId;
  // vector<UInt_t>           TPHdata;
  // vector<UShort_t>         xHVdata, yHVdata; // in Volt
  // vector<UShort_t>         xHCdata, yHCdata; // in nano-Amp
  // vector<vector<UInt_t>>   NOISEdata;	// one vector<UInt_t> per layer
    
public:
  HealthStore() {};
  virtual ~HealthStore() {};
  void ClearMon();
  
  MonOpt1  *AddMon1Store();
  TClonesArray *GetMonOpt1() const {return fMonOpt1;}
  
  ClassDef(HealthStore,1);  //Event structure 
};

#endif
 
