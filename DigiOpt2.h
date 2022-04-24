#ifndef ROOT_digi_opt2
#define ROOT_digi_opt2

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <vector>
using namespace std;
class DigiOpt2 : public TObject {
 public:
  
  UInt_t SiPM_data;
  //     2 bit for SiPM
  //     8 bit for bars (-Y to +Y for top, down to up for rest)
  //     2 bit for layers (inside to outside)
  //     4 bit for sides (Front:0, Top:1, back:2, left:3, right:4)
  //    20 bit for timing
  //    10 bit for Low gain
  //    10 bit for High gain

  DigiOpt2(){};
  // DigiOpt2(Float_t random){random=0;};
  
  virtual ~DigiOpt2(){};
  ClassDef(DigiOpt2,1);
};

#endif
