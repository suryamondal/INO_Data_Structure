#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMath.h"   

#include "CauStore.h"

using namespace std;

ClassImp(CauOpt1);


TClonesArray *CauStore::fgCauOpt1 = 0;

CauStore::CauStore() {
  // Create an Hit object.
  if (!fgCauOpt1) fgCauOpt1 = new TClonesArray("CauOpt1", 30000);
  fCauOpt1 = fgCauOpt1;
  NCauOpt1=0;
}

//______________________________________________________________________________
CauStore::~CauStore() {
  
}

//______________________________________________________________________________
CauOpt1 * CauStore::AddCau1Store() {
  // Add a new track to the list of tracks for this event.
  TClonesArray &digi = *fCauOpt1;
  CauOpt1 *digipos = new(digi[NCauOpt1++]) CauOpt1();
  return digipos;
}

void CauStore::ClearCau() {
  NCauOpt1 = 0; // !!! reset counter !!!
  if(fCauOpt1) {
    // fCauOpt1->Clear("C"); // will also call Track::Clear
    fCauOpt1->Delete(); // will also call Track::Clear
  }
}

