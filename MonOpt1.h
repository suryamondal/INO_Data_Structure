#ifndef ROOT_Mon_opt1
#define ROOT_Mon_opt1

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <vector>
using namespace std;

class MonOpt1 : public TObject {
 public:
  
  UShort_t         rpcId;
  UInt_t           TPHdata;
  Long64_t         HVdata;
  UShort_t         NOISEdata[128]; // noise/sec !!! careful about overflow

  // UShort_t         xHVdata, yHVdata; // in Volt
  // UShort_t         xHCdata, yHCdata; // in nano-Amp
  // vector<UInt_t>   NOISEdata;

  MonOpt1(){};
  // MonOpt1(Float_t random){random=0;};
  
  virtual ~MonOpt1(){};
  ClassDef(MonOpt1,1);
};

#endif
