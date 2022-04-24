#ifndef ROOT_CauStore
#define ROOT_CauStore

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include <vector>

#include "CauOpt1.h"

//#include <pair>
using namespace std;

/* 

 */


class CauStore : public TObject {
 public:
  
  TTimeStamp      Cautime;
  Int_t           CauENum;

  UInt_t NCauOpt1;
  TClonesArray  *fCauOpt1;
  static TClonesArray *fgCauOpt1;
  
  // Int_t           select_line;
  // ULong64_t       raw_trig_cnt;
  // ULong64_t       acpt_trig_cnt;
  // Int_t           cau_status;
  // ULong64_t       cau_ref_l[2];	// ref1 & ref2
  // ULong64_t       cau_ref_t[2];	// ref1 & ref2

  // vector<ULong64_t> cau_tdc_ref1;

  // vector<ULong64_t> cau_tdc;
  // element is <cau tdc time> <l or t: 1bit> <RPC Id: 16bit>
  
public:
  CauStore() {};
  virtual ~CauStore() {};
  void ClearCau();// clears previous track objects..
  
  CauOpt1  *AddCau1Store();
  TClonesArray *GetCauOpt1() const {return fCauOpt1;}

  ClassDef(CauStore,1);  //Event structure 
};

#endif
 
