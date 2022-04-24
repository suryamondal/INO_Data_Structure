#ifndef ROOT_Cau_opt1
#define ROOT_Cau_opt1

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <vector>
using namespace std;

class CauOpt1 : public TObject {
 public:
  
  UShort_t     rpcId;
  UInt_t       cau_tdc;
  
  CauOpt1(){};
  // CauOpt1(Float_t random){random=0;};
  
  virtual ~CauOpt1(){};
  ClassDef(CauOpt1,1);
};

#endif
