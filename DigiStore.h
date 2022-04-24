#ifndef ROOT_DigiStore
#define ROOT_DigiStore

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include "DigiOpt1.h"		// RPC
#include "DigiOpt2.h"		// SiPM

using namespace std;

/* 

 */


class DigiStore : public TObject {
  
public:

  TTimeStamp EventTime;
  UInt_t ENum;			// Event Number
  // UInt_t REnum;			// Event Number
  // UInt_t CEnum;
  
  UInt_t NDigiOpt1;		// Number of DigiStore not sets!
  TClonesArray  *fDigiOpt1;	//->array with all DigiStore
  static TClonesArray *fgDigiOpt1;

  UInt_t NDigiOpt2;		// Number of DigiStore not sets!
  TClonesArray  *fDigiOpt2;	//->array with all DigiStore
  static TClonesArray *fgDigiOpt2;
  
  
public:

  DigiStore();
  virtual ~DigiStore();
  void ClearDigi();// clears previous track objects..
  
  DigiOpt1  *AddDigi1Store();
  TClonesArray *GetDigiOpt1() const {return fDigiOpt1;}
  
  DigiOpt2  *AddDigi2Store();
  TClonesArray *GetDigiOpt2() const {return fDigiOpt2;}
  
  ClassDef(DigiStore,1);  //Event structure 
};

#endif
 
