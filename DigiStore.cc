#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMath.h"   

#include "DigiStore.h"

using namespace std;

ClassImp(DigiOpt1);

/* Option 2 */
ClassImp(DigiOpt2);


TClonesArray *DigiStore::fgDigiOpt1 = 0;
TClonesArray *DigiStore::fgDigiOpt2 = 0;

DigiStore::DigiStore() {
  // Create an Hit object.
  if (!fgDigiOpt1) fgDigiOpt1 = new TClonesArray("DigiOpt1", 1000);
  fDigiOpt1 = fgDigiOpt1;
  NDigiOpt1=0;

  /* Option 2 */
  if (!fgDigiOpt2) fgDigiOpt2 = new TClonesArray("DigiOpt2", 1000);
  fDigiOpt2 = fgDigiOpt2;
  NDigiOpt2=0;
}

//______________________________________________________________________________
DigiStore::~DigiStore() {
  
}

//______________________________________________________________________________
DigiOpt1 * DigiStore::AddDigi1Store() {
  // Add a new track to the list of tracks for this event.
  TClonesArray &digi = *fDigiOpt1;
  DigiOpt1 *digipos = new(digi[NDigiOpt1++]) DigiOpt1();
  return digipos;
}

/* Option 2 */
DigiOpt2 * DigiStore::AddDigi2Store() {
  // Add a new track to the list of tracks for this event.
  TClonesArray &digi = *fDigiOpt2;
  DigiOpt2 *digipos = new(digi[NDigiOpt2++]) DigiOpt2();
  return digipos;
}

void DigiStore::ClearDigi() {
  NDigiOpt1 = 0; // !!! reset counter !!!
  if(fDigiOpt1) {
    // fDigiOpt1->Clear("C"); // will also call Track::Clear
    fDigiOpt1->Delete(); // will also call Track::Clear
  }
  ENum = 0;			// Event Number
  REnum = 0;			// Event Number
  CEnum = 0;

  /* Option 2 */
  NDigiOpt2 = 0; // !!! reset counter !!!
  if(fDigiOpt2) {
    // fDigiOpt2->Clear("C"); // will also call Track::Clear
    fDigiOpt2->Delete(); // will also call Track::Clear
  }
}

