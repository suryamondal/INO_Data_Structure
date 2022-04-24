#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMath.h"   

#include "HealthStore.h"

using namespace std;

ClassImp(MonOpt1);

TClonesArray *HealthStore::fgMonOpt1 = 0;

HealthStore::HealthStore() {
  // Create an Hit object.
  if (!fgMonOpt1) fgMonOpt1 = new TClonesArray("MonOpt1", 1000);
  fMonOpt1 = fgMonOpt1;
  NMonOpt1=0;
}

//______________________________________________________________________________
HealthStore::~HealthStore() {
  
}

//______________________________________________________________________________
MonOpt1 * HealthStore::AddMon1Store() {
  // Add a new track to the list of tracks for this event.
  TClonesArray &digi = *fMonOpt1;
  MonOpt1 *digipos = new(digi[NMonOpt1++]) MonOpt1();
  return digipos;
}

void HealthStore::ClearMon() {
  NMonOpt1 = 0; // !!! reset counter !!!
  if(fMonOpt1) {
    // fMonOpt1->Clear("C"); // will also call Track::Clear
    fMonOpt1->Delete(); // will also call Track::Clear
  }
}

