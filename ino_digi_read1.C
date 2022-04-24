

/*
  20210320: Removing TimeGroup
          : Renaming TimeGroup1 to TimeGroup
	  : Removing hit[nside][nlayer] array


File: GMA_RPCv4t_evtraw-20181217_112935.root

Evt: 38, 75, 117, 201, 245

16, 17, 21, 38, 39

Cluster seperation test: 10, 21, 25, 28, 

Good Event: 6, 10, 115
Scattering: 15,


 */


#define isDebug
// #define isSpclDebug

#define isMiniICAL



#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <bitset>

#include "TTimeStamp.h"
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
#include "TVector2.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMinuit.h"

#include "RunInfo.h"
#include "DigiStore.h"
#include "CauStore.h"
#include "HealthStore.h"


using namespace std;

/*
  
  
*/

const int        nside         =   2;
const int        nlayer        =  10;
const int        nxrow         =   1;
const int        nyrow         =   1;
const int        nmodule       =   1;
const int        nstrip        =  64;
const int        nTDC          =   8;
const double     tdc_least     =   0.1; // in ns

const int        blockM        =   nstrip/4;
const double     strpwidth     = 0.03; // in m

const int        MaxTotHit     =   nstrip/2;
const int        MaxTotPos     =   3;
const int        MaxMulti      =   3;
const int        MinLayHit     =   7;
const int        MinClusterSep =   0;

const double     maxtime       =  22.e3;      // in ns
const double     cval          =   0.29979;   // light speed in m/ns


TVector3 ttxx, ttyy;
double calPointDist(TVector3 xx, TVector3 yy) {
  double dist = sqrt(pow(xx.X()-yy.X(),2.) +
		     pow(xx.Y()-yy.Y(),2.) +
		     pow(xx.Z()-yy.Z(),2.));
  return dist;
};


#ifdef isMiniICAL
const double   gapThickness      = 0.008; // 
const double   ironThickness     = 0.056; // 
const double   airGap            = 0.045; // 
const double   rpcXdistance      = 2.;	  // 
const double   rpcYdistance      = 2.1;	  // 
const double   rpcXOffset        = 0.;	  // 
const double   rpcYOffset        = 0.;	  // 
const double   moduleDistance    = 0.;	  // 

TVector3 getRPCpos(Int_t module,
		   Int_t xrow,
		   Int_t yrow,
		   Int_t zlay) {
  double zpos = (zlay + 1)*(ironThickness + airGap) - airGap/2.;
  double ypos = rpcYdistance*yrow + rpcYOffset;
  double xpos = rpcXdistance*xrow + rpcXOffset;
  TVector3 pos(xpos,ypos,zpos);
  return pos;			// in meter
};


#endif


void getRPCId(UShort_t  rpcId,
	      Int_t    &module,
	      Int_t    &xrow,
	      Int_t    &yrow,
	      Int_t    &zlay) {
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

void getTDC(UInt_t  tdc,
	    Bool_t &isTrail,	// leading or trailing
	    Bool_t &side,	// x or y
	    Int_t  &tdcno,	// tdc no 0-7
	    Int_t  &tdcval) {
  tdcno     = (tdc   )&0b1111;
  isTrail   = (tdc>>4)&0b1;
  tdcval    = (tdc>>5);
  side      = (tdcno>=nTDC)?1:0;
  if(tdcno>=nTDC) {tdcno-=nTDC;}
};

void getCAU(ULong64_t  tdc,
	    Bool_t    &isTrail,	// leading or trailing
	    UShort_t  &rpcId,
	    Int_t     &tdcval) {
  rpcId     = (tdc    )&0xFF;
  isTrail   = (tdc>>16)&0b1;
  tdcval    = (tdc>>17);
};

struct CauInfo {
  TTimeStamp       Cautime;
  vector<UShort_t> rpcId;
  vector<Double_t> val;
};

struct rawstrp {
  Int_t            strp;
  vector<Double_t> tdc[2];	// leading & trailing
};

struct rawlayer {
  int module, xrow, yrow, zlay;
  vector<rawstrp>         hit[nside];
  double tdc_ref[2];		// leading and trailing
  vector<vector<rawstrp>> cluster[nside];
};

struct PosInSpace {
  rawlayer     rawhitInfo;
  TVector3     RawHitPos;
  TVector3     PosCorrection;
};

struct GroupInfo {
  vector<rawlayer>   alllayer;
  vector<PosInSpace> PosForTrack;
  vector<vector<PosInSpace>> HoughCluster;
};


void calHitPos(PosInSpace &pos) {

  int nm = pos.rawhitInfo.module;
  int nx = pos.rawhitInfo.xrow;
  int ny = pos.rawhitInfo.yrow;
  int nl = pos.rawhitInfo.zlay;
  
  pos.RawHitPos = getRPCpos(nm, nx, ny, nl)*(1./strpwidth);
  
  int multi[nside];
  double calps[nside];
  double calps1[nside];
  for(int nj=0;nj<nside;nj++) {
    multi[nj] = pos.rawhitInfo.cluster[nj][0].size();
    calps[nj] = calps1[nj] = (pos.rawhitInfo.cluster[nj][0][0].strp
			      + multi[nj]*0.5);
    // cout << " lay " << pos.rawhitInfo.zlay << " " << nj
    // 	 << " " << calps[nj] << endl;
  } // for(int nj=0;nj<nside;nj++) {
  
  // TVector3 posxx(calps[0],calps[1],0);
  TVector3 posxx(calps1[0],calps1[1],0);
  pos.RawHitPos += posxx;

  // int blkx = int(posxx.X()*blockM/nstrip);
  // int blky = int(posxx.Y()*blockM/nstrip);
  // pos.PosCorrection.SetX(blockposoff[nm][nx][ny][nl][0][blkx][blky]);
  // pos.PosCorrection.SetY(blockposoff[nm][nx][ny][nl][1][blkx][blky]);
  // pos.PosCorrection.SetZ(0.);
  
};


void LinearVectorFit(bool              isTime,
		     vector<TVector3>  pos,
		     vector<TVector2>  poserr,
		     vector<bool>      occulay,
		     TVector2         &slope,
		     TVector2         &inter,
		     TVector2         &chi2,
		     vector<TVector3> &ext,
		     vector<TVector2> &exterr) {
  
  double szxy[nside] = {0};
  double   sz[nside] = {0};
  double  sxy[nside] = {0};
  double   sn[nside] = {0};
  double  sz2[nside] = {0};
  
  double     slp[nside] = {-10000,-10000};
  double  tmpslp[nside] = {-10000,-10000};
  double intersect[nside] = {-10000,-10000};
  double    errcst[nside] = {-10000,-10000};
  double    errcov[nside] = {-10000,-10000};
  double    errlin[nside] = {-10000,-10000};
  
  for(int ij=0;ij<int(pos.size());ij++) {
    if(int(occulay.size()) && !occulay[ij]) {continue;}
    // cout << " ij " << ij << endl;
    double xyzval[3] = {pos[ij].X(),
			pos[ij].Y(),
			pos[ij].Z()};
    double xyerr[2]  = {poserr[ij].X(),
			poserr[ij].Y()};
    for(int nj=0;nj<nside;nj++) {
      szxy[nj] += xyzval[2]*xyzval[nj]/xyerr[nj];
      sz[nj]   += xyzval[2]/xyerr[nj];
      sz2[nj]  += xyzval[2]*xyzval[2]/xyerr[nj];
      sxy[nj]  += xyzval[nj]/xyerr[nj];
      sn[nj]   += 1/xyerr[nj];
    }   // for(int nj=0;nj<nside;nj++) {
  } // for(int ij=0;ij<int(pos.size());ij++){
  
  for(int nj=0;nj<nside;nj++) {
    if(sn[nj]>0. && sz2[nj]*sn[nj] - sz[nj]*sz[nj] !=0.) { 
      slp[nj] = (szxy[nj]*sn[nj] -
		   sz[nj]*sxy[nj])/(sz2[nj]*sn[nj] - sz[nj]*sz[nj]);
      tmpslp[nj] = slp[nj]; 
      if(isTime) { //time offset correction
        // if(fabs((cval*1.e-9)*slope+1)<3.30) { 
	tmpslp[nj] = -1./cval;
	// }
      }
      intersect[nj] = sxy[nj]/sn[nj] - tmpslp[nj]*sz[nj]/sn[nj];

      double determ = (sn[nj]*sz2[nj] - sz[nj]*sz[nj]);
      errcst[nj] = sz2[nj]/determ;
      errcov[nj] = -sz[nj]/determ;
      errlin[nj] = sn[nj]/determ;
    }
  } // for(int nj=0;nj<nside;nj++) {
  slope.SetX(tmpslp[0]);
  slope.SetY(tmpslp[1]);
  inter.SetX(intersect[0]);
  inter.SetY(intersect[1]);
  
  // theta = atan(sqrt(pow(tmpslp[0],2.)+pow(tmpslp[1],2.)));
  // phi = atan2(tmpslp[1],tmpslp[0]);
  
  double sumx = 0, sumy = 0;
  ext.clear(); exterr.clear();
  TVector3 xxt;
  TVector2 xxtt;
  for(int ij=0;ij<int(pos.size());ij++){
    xxt.SetX(tmpslp[0]*pos[ij].Z()+intersect[0]);
    xxt.SetY(tmpslp[1]*pos[ij].Z()+intersect[1]);
    xxt.SetZ(pos[ij].Z());
    ext.push_back(xxt);
    xxtt.SetX(errcst[0] + 2*errcov[0]*pos[ij].Z()+
	     errlin[0]*pos[ij].Z()*pos[ij].Z());
    xxtt.SetY(errcst[1] + 2*errcov[1]*pos[ij].Z()+
	     errlin[1]*pos[ij].Z()*pos[ij].Z());
    exterr.push_back(xxtt);
    // cout << " " << int(exterr.size())
    // 	 << " " << 1./exterr.back().X()
    // 	 << " " << 1./exterr.back().Y() << endl;
    if(int(occulay.size())==0 || occulay[ij]) {
      sumx += pow(xxt.X()-pos[ij].X(), 2.)/poserr[ij].X(); 
      sumy += pow(xxt.Y()-pos[ij].Y(), 2.)/poserr[ij].Y(); 
    }
  } // for(int ij=0;ij<int(pos.size());ij++){
  chi2.SetX(sumx);
  chi2.SetY(sumy);
};


/*
  This retuns the vector for minimum distance between two lines
  represented by (pos0, dir0) and (pos1, dir1).
  Note: Here dir0 and dir1 are unit vectors.
  
  This function can also be used to get minimum distance between
  a line and a point by making dir0=dir1.
*/
TVector3 getMinDist(TVector3 pos0,
		    TVector3 dir0,
		    TVector3 pos1,
		    TVector3 dir1) {
  
  double l0 = dir0.X();
  double m0 = dir0.Y();
  double n0 = dir0.Y();
  double l1 = dir1.X();
  double m1 = dir1.Y();
  double n1 = dir1.Z();
  
  double alpha =  (l0*(pos0.X() - pos1.X()) +
		   m0*(pos0.Y() - pos1.Y()) +
		   n0*(pos0.Z() - pos1.Z()));
	      
  double beta  =  (l1*(pos0.X() - pos1.X()) +
		   m1*(pos0.Y() - pos1.Y()) +
		   n1*(pos0.Z() - pos1.Z()));
  
  double gamma =  l0*l1 + m0*m1 + n0*n1;
  
  double r0 = (-alpha +  beta*gamma)/(1. - gamma*gamma);
  double r1 = ( beta  - alpha*gamma)/(1. - gamma*gamma);
  
  double px0 = pos0.X() + l0*r0;
  double py0 = pos0.Y() + m0*r0;
  double pz0 = pos0.Z() + n0*r0;
	      
  double qx1 = pos1.X() + l1*r1;
  double qy1 = pos1.Y() + m1*r1;
  double qz1 = pos1.Z() + n1*r1;
  
  TVector3 minDistV(qx1-px0,qy1-py0,qz1-pz0);
  return minDistV;
  
};				// TVector3 getMinDist(







int main(int argc, char** argv) {

  const char *sideMark[nside] = {"x","y"};

  Long64_t start_s = clock();
  
  CauStore  *cau   = 0;
  DigiStore *event = 0;
  
  CauInfo tmpCau;
  vector<CauInfo> allCAUevents;
  double cauCorrections[nmodule][nxrow][nyrow][nlayer] = {0};
  
  rawstrp            tmphit;
  vector<rawstrp>    tmpCluster;
  rawlayer           tmplayer;
  vector<rawlayer>   allrawlay;

  GroupInfo          tmpevt;
  vector<GroupInfo>  TimeGroup;
  
  PosInSpace         tmpPos;
  
  Double_t           tdc_ref[nlayer][2]; // leading and trailing
  Double_t           ttExecTimeVal;
  
  
  char name[300];
  
  TH1D *strptimeraw = new TH1D("strptimeraw","strptimeraw",
			       int(maxtime/10.),0.,maxtime);
  TH1D *hbinIns = new TH1D("hbinIns","hbinIns",
			   50,0.5,50.5);
  
  char outfil[300] = {}; 
  char outfilx[300] = {};
  int len = strlen(argv[1]);
  strncpy(outfil, argv[1], len-5); // root file
  sprintf(outfilx, "../rredata/temp/%s_o_%i.root", outfil, stoi(argv[4]));
  // sprintf(outfilx, "..//miniICAL_data_GMA/temp/%s_o_%i.root", outfil, stoi(argv[4]));
  TFile *f1 = new TFile(outfilx,"RECREATE");



  
  /*
    C : a character string terminated by the 0 character
    B : an 8 bit signed integer (Char_t)
    b : an 8 bit unsigned integer (UChar_t)
    S : a 16 bit signed integer (Short_t)
    s : a 16 bit unsigned integer (UShort_t)
    I : a 32 bit signed integer (Int_t)
    i : a 32 bit unsigned integer (UInt_t)
    F : a 32 bit floating point (Float_t)
    f : a 24 bit floating point with truncated mantissa (Float16_t)
    D : a 64 bit floating point (Double_t)
    d : a 24 bit truncated floating point (Double32_t)
    L : a 64 bit signed integer (Long64_t)
    l : a 64 bit unsigned integer (ULong64_t)
    G : a long signed integer, stored as 64 bit (Long_t)
    g : a long unsigned integer, stored as 64 bit (ULong_t)
    O : [the letter o, not a zero] a boolean (Bool_t)
  */
  
  
  
  char datafile[300] = {};
  strncpy(datafile,argv[1],300);
  Long64_t nentrymx = stoi(argv[3]);
  Long64_t nentrymn = stoi(argv[2]);
  char infile[300]  = {};

  sprintf(infile, "../rredata/%s", datafile);
  TFile *fileIn = new TFile(infile, "read");
  
  if(!fileIn->IsZombie()) {
    


    
    /******  Run Info   ******/
    
    RunInfo *runinfo = (RunInfo*)fileIn->Get("run");
    cout << " created by: " << runinfo->Creator << endl;
    cout << " on: " << runinfo->CreatedOn << endl;
    cout << " trigger: " << runinfo->TriggerInfo << endl;








    /******  CAU Data    ******/
    
    TTree *cau_tree = (TTree*)fileIn->Get("CAUtree");
    cau_tree->SetBranchAddress("CauData", &cau);
    Long64_t Cnentry = cau_tree->GetEntries();
    cout << " cau entry " << Cnentry << endl;

    allCAUevents.clear();
    
    for(Long64_t iev=0;iev<Cnentry;iev++) {
      
      tmpCau.rpcId.clear();
      tmpCau.val.clear();
      
      fileIn->cd();
      cau_tree->GetEntry(iev);

      
      // CauInfo tmpCau;
      // vector<CauInfo> allCAUevents;
      
      // cout << " " << cau->Cautime
      // 	   << " " << cau->cau_tdc_ref1[0]//.size()
      // 	   << " " << cau->cau_tdc[0]//.size()
      // 	   << endl;
      // cout << " " << cau->Cautime << endl;
      
      tmpCau.Cautime = cau->Cautime;
      for(int ij=0;ij<int(cau->cau_tdc_ref1.size());ij++) {
	int tdcCnt; UShort_t rpcId; Bool_t lead;
	getCAU(cau->cau_tdc_ref1[ij],lead,rpcId,tdcCnt);
	// cout << " " << rpcId << " " << lead << " " << tdcCnt << endl;
	if(lead==0) {		// using only leading
	  if(int(tmpCau.rpcId.size())==0
	     || (int(tmpCau.rpcId.size())
		 && tmpCau.rpcId.back()!=rpcId)) {
	    tmpCau.rpcId.push_back(rpcId);
	    tmpCau.val.push_back(tdcCnt*tdc_least);
	  }
	}
      }	// for(int ij=0;ij<int(cau->cau_tdc_ref1.size());ij++) {
      allCAUevents.push_back(tmpCau);
    } // for(Long64_t iev=0;iev<Cnentry;iev++) {

// #ifdef isDebug
//     for(Long64_t iev=0;iev<int(allCAUevents.size());iev++) {
//       cout << " " << allCAUevents[iev].Cautime << endl;
//       for(int ij=0;ij<int(allCAUevents[iev].rpcId.size());ij++) {
// 	cout << "   " << allCAUevents[iev].rpcId[ij]
// 	     << " " << allCAUevents[iev].val[ij]
// 	     << endl;
//       }
//     } // for(Long64_t iev=0;iev<Cnentry;iev++) {
// #endif  // #ifdef isDebug
    


    




    /******  Event Data    ******/
    
    TTree *event_tree = (TTree*)fileIn->Get("RPCtree");
    event_tree->SetBranchAddress("EventData", &event);
    
    Long64_t nentry = event_tree->GetEntries();
    nentry = min(nentry,nentrymx);
    cout << datafile <<" has "<<nentry<<" events "<<endl;
    
    for(Long64_t iev=nentrymn;iev<nentry;iev++) {
      
#ifdef isDebug
      cout << endl << endl << " " << iev << endl;
#endif  // #ifdef isDebug
      
      if(iev%1000==0) {
    	Long64_t stop_s = clock();
    	cout << " iev " << iev
    	     << " time " << (stop_s-start_s)/Double_t(CLOCKS_PER_SEC)
    	     << endl;
      }
      
      ttExecTimeVal = 100000000;
      // for(int nj=0;nj<nside;nj++) {
      // 	for(int ij=0;ij<nlayer;ij++) {
      // 	  hit[nj][ij].clear();
      // 	}
      // }
      allrawlay.clear();
      
      fileIn->cd();
      event_tree->GetEntry(iev);
      
      TTimeStamp eventTime = event->EventTime;
#ifdef isDebug
      cout << " time " << eventTime << endl;
#endif  // #ifdef isDebug


      
      TClonesArray *tmpTObj = (TClonesArray *)event->GetDigiOpt1();
      // DigiOpt1 *tmpTObj = (DigiOpt1 *)event->GetDigiOpt1();
      // int objentry = event->GetDigiOpt1()->GetEntries();
      int objentry = tmpTObj->GetEntries();
#ifdef isDebug
      cout << " event entries : " << objentry << endl;
#endif  // #ifdef isDebug
      
      for(int iobj=0;iobj<objentry;iobj++) {
	
	DigiOpt1 *tmpDigiOpt1 = (DigiOpt1 *)tmpTObj->At(iobj);
#ifdef isDebug
	cout << " " << bitset<16>(tmpDigiOpt1->rpcId) << endl
	     << " " << bitset<64>(tmpDigiOpt1->strips[0])
	     << " " << bitset<64>(tmpDigiOpt1->strips[1])
	     << " " <<        int(tmpDigiOpt1->tdc.size())
	     << endl;
#endif  // #ifdef isDebug
	
	Int_t module, xrow, yrow, zlay;
	getRPCId(tmpDigiOpt1->rpcId,module,xrow,yrow,zlay);
#ifdef isDebug
	cout << " rpcId " << module
	     << " " << xrow
	     << " " << yrow
	     << " " << zlay
	     << endl;
#endif  // #ifdef isDebug
	tmplayer.module = module;
	tmplayer.xrow = xrow;
	tmplayer.yrow = yrow;
	tmplayer.zlay = zlay;
	
	tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
	for(int nj=0;nj<nside;nj++) {
	  tmplayer.hit[nj].clear();
	  for(int kl=nstrip-1; kl>=0; kl--) {
	    if(((tmpDigiOpt1->strips[nj])>>kl)&0b1) {
	      tmphit.strp = kl;
	      // hit[nj][zlay].push_back(tmphit);
	      tmplayer.hit[nj].push_back(tmphit);
	    }
	  } // for(int kl=nstrip-1; kl>=0; kl--) {
	} // for(int nj=0;nj<nside;nj++) {
	
	tmplayer.tdc_ref[0] = tdc_least*tmpDigiOpt1->tdc[0];
	tmplayer.tdc_ref[1] = tdc_least*tmpDigiOpt1->tdc[1];
	
#ifdef isDebug
	cout << " tdc_ref " << tmplayer.tdc_ref[0]
	     << " "         << tmplayer.tdc_ref[1] << endl;
#endif  // #ifdef isDebug
	for(int kl=2;kl<int(tmpDigiOpt1->tdc.size());kl++) {
	  Int_t tdcno, tdcval; Bool_t side; Bool_t lead;
	  getTDC(tmpDigiOpt1->tdc[kl],lead,side,tdcno,tdcval);
#ifdef isDebug
	  cout << " tdc " << lead
	       << " " << side
	       << " " << tdcno
	       << " " << tdcval*tdc_least
	       << endl;
#endif  // #ifdef isDebug
	  
	  for(int ki=0;ki<int(tmplayer.hit[side].size());ki++) {
	    if(tmplayer.hit[side][ki].strp%nTDC==tdcno) {
	      tmplayer.hit[side][ki].tdc[lead].push_back(tdcval*tdc_least);
	      if(lead==0) {
		Double_t tdctime = tdcval*tdc_least - tmplayer.tdc_ref[0];
		if(ttExecTimeVal>tdctime) {
		  ttExecTimeVal = tdctime;}
	      }
	    }
	  } // for(int ki=0;ki<int(tmplayer.hit[side][zlay].size());ki++) {
	} // for(int kl=2;kl<int(tmpDigiOpt1->tdc.size());kl++) {
	allrawlay.push_back(tmplayer);
      }	// for(int iobj=0;iobj<objentry;iobj++) {

      
      /* Making all tdc_l +ve w.r.t. tdc_ref_l */
      for(int nj=0;nj<nside;nj++) {
	for(int nl=0;nl<int(allrawlay.size());nl++) {
	  for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
	    for(int ntd=0;ntd<2;ntd++) {
	      for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[ntd].size());nht++) {
		allrawlay[nl].hit[nj][nh].tdc[ntd][nht] -= (ttExecTimeVal-100.);
	      }
	    } // for(int ntd=0;ntd<2;ntd++) {
	  }
	} // for(int nl=0;nl<nlayer;nl++) {
      }	// for(int nj=0;nj<nside;nj++) {
      
#ifdef isDebug
      cout << " min time " << ttExecTimeVal << endl;
      for(int nj=0;nj<nside;nj++) {
	cout << " " << sideMark[nj] << endl;
	for(int nl=0;nl<int(allrawlay.size());nl++) {
	  cout << "  l" << allrawlay[nl].zlay << " ";
	  for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
	    cout << " " << allrawlay[nl].hit[nj][nh].strp;
	    cout << " (";
	    for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[0].size());nht++) {
	      cout << " " << allrawlay[nl].hit[nj][nh].tdc[0][nht] - allrawlay[nl].tdc_ref[0];
	    }
	    cout << " ) ";
	  }
	  cout << endl;
	}
      }
#endif  // #ifdef isDebug
      
      
      
      /*****  Grouping Hits as per time info *****/
      /*****  Using only leading time stamp  *****/
      
      TimeGroup.clear();
      
      strptimeraw->Scale(0);
      for(int nj=0;nj<nside;nj++) {
	for(int nl=0;nl<int(allrawlay.size());nl++) {
	  for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
	    for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[0].size());nht++) {
	      double tdctime = allrawlay[nl].hit[nj][nh].tdc[0][nht]-allrawlay[nl].tdc_ref[0];
	      strptimeraw->Fill(tdctime);
	    }
	  }
	} // for(int nl=0;nl<int(allrawlay.size());nl++) {
      }	// for(int nj=0;nj<nside;nj++) {
      
      strptimeraw->SetBinContent(strptimeraw->GetNbinsX(),0);
      int bnx = 0;
      while(bnx<strptimeraw->GetNbinsX()) {
	if((strptimeraw->GetBinContent(bnx+1))>0) {
	  double bStart = strptimeraw->GetBinLowEdge(bnx+1);
	  while((strptimeraw->GetBinContent(bnx+1))>0) {
	    bnx++;
	  }
	  double bEnd = (strptimeraw->GetBinWidth(bnx)
			 + strptimeraw->GetBinLowEdge(bnx));
	  // cout << " " << bStart << " " << bEnd << endl;
	  tmpevt.alllayer.clear();
	  for(int nl=0;nl<int(allrawlay.size());nl++) {
	    tmplayer = allrawlay[nl];
	    for(int nj=0;nj<nside;nj++) {
	      tmplayer.hit[nj].clear(); tmplayer.hit[nj].clear();
	      for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
	  	tmphit.strp = allrawlay[nl].hit[nj][nh].strp;
	  	tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
	  	for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[0].size());nht++) {
	  	  double tdctime = allrawlay[nl].hit[nj][nh].tdc[0][nht]-allrawlay[nl].tdc_ref[0];
	  	  if(tdctime>=bStart && tdctime<=bEnd) {
	  	    tmphit.tdc[0].push_back(allrawlay[nl].hit[nj][nh].tdc[0][nht]);
	  	    if(nht<int(allrawlay[nl].hit[nj][nh].tdc[1].size())) {
	  	      tmphit.tdc[1].push_back(allrawlay[nl].hit[nj][nh].tdc[1][nht]);
	  	    }
	  	  }
	  	} // for(int nht=0;nht<int(hit[nj][nl][nh].tdc[0].size());nht++)
	  	if(int(tmphit.tdc[0].size())>0) {
	  	  tmplayer.hit[nj].push_back(tmphit);
	  	}
	      }	// for(int nh=0;nh<int(hit[nj][nl].size());nh++) {
	    }   // for(int nj=0;nj<nside;nj++) {
	    if(int(tmplayer.hit[0].size())>0
	       || int(tmplayer.hit[1].size())>0) {
	      tmpevt.alllayer.push_back(tmplayer);
	    }
	  } // for(int nl=0;nl<int(allrawlay.size());nl++) {
	  TimeGroup.push_back(tmpevt);
	} else {bnx++;}
      }	// while(1) {
      
#ifdef isDebug
      cout << " TimeGroup " << " " << TimeGroup.size() << endl;
      for(int isj=0;isj<int(TimeGroup.size());isj++) {
	cout << " TimeGroup " << isj << endl;
	for(int nj=0;nj<nside;nj++) {
	  cout << "   " << sideMark[nj] << endl;
	  for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	    tmplayer = TimeGroup[isj].alllayer[nl];
	    cout << "    l" << tmplayer.zlay << " ";
	    for(int nh=0;nh<int(tmplayer.hit[nj].size());nh++) {
	      cout << " " << tmplayer.hit[nj][nh].strp;
	      cout << " (";
	      for(int nht=0;nht<int(tmplayer.hit[nj][nh].tdc[0].size());nht++) {
		cout << " " << tmplayer.hit[nj][nh].tdc[0][nht] - tmplayer.tdc_ref[0];
	      }
	      cout << " ) ";
	    }
	    cout << endl;
	  }
	}
      }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {
#endif  // #ifdef isDebug
      
      // cout << " iev " << iev << endl;

      
      
      /*****  Forming Cluster of Strips in each Layer *****/
      
      for(int isj=0;isj<int(TimeGroup.size());isj++) {
	
	for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	  tmplayer = TimeGroup[isj].alllayer[nl];
	  
	  for(int nj=0;nj<nside;nj++) {
	    // TimeGroup[isj].alllayer[nl].cluster[nj].clear();
	    tmplayer.cluster[nj].clear();
	    int totHit = tmplayer.hit[nj].size();
	    // if(totHit>MaxTotHit) {continue;}
	    
	    tmpCluster.clear();
	    int tempX = 0;
	    for(int jk=0;jk<totHit;jk++) {
	      int lowX = tmplayer.hit[nj].back().strp;
	      tmpCluster.push_back(tmplayer.hit[nj].back());
	      tmplayer.hit[nj].insert(tmplayer.hit[nj].begin(),
				      tmplayer.hit[nj].back());
	      tmplayer.hit[nj].pop_back();
	      tempX++;
	      if((lowX+1!=tmplayer.hit[nj].back().strp)
		 ||(jk+1==totHit)) {
		// if(tempX<=MaxMulti) {
		//   tmplayer.cluster[nj].push_back(tmpCluster);
		// }
		tmplayer.cluster[nj].push_back(tmpCluster);
		tempX=0;
		tmpCluster.clear();
	      }
	    } // for(int jk=0;jk<totHit;jk++) {
	    // if(tmplayer.cluster[nj].size()>MaxTotPos) {tmplayer.cluster[nj].clear();}
	    
	    /*** checking gap between clusters ***/
	    for(int jk=0;jk<int(tmplayer.cluster[nj].size())-1;jk++) {
	      for(int ki=jk+1;ki<int(tmplayer.cluster[nj].size());) {
		// cout << " some1 " << jk << " " << ki << endl;
		int startPos = tmplayer.cluster[nj][jk].back().strp;
		int endPos   = tmplayer.cluster[nj][ki].front().strp;
		if(endPos-startPos<=MinClusterSep+1) {
		  for(int kii=startPos+1;kii<endPos;kii++) {
		    // filling gap strp without tdc entry
		    tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
		    tmphit.strp = kii;
		    tmplayer.cluster[nj][jk].push_back(tmphit);
		  }
		  for(int kii=0;kii<int(tmplayer.cluster[nj][ki].size());kii++) {
		    tmplayer.cluster[nj][jk].push_back(tmplayer.cluster[nj][ki][kii]);
		  }
		  tmplayer.cluster[nj].erase(tmplayer.cluster[nj].begin()+ki);
		} else {ki++;}
		// cout << " something " << iev
		//      << " " << tmplayer.cluster[nj][jk].back().strp
		//      << " " << tmplayer.cluster[nj][ki].front().strp
		//      << endl;
	      }	// for(int ki=jk+1;ki<int(tmplayer.cluster[nj].size());) {
	      // if(int(tmplayer.cluster[nj][jk].size())>MaxMulti) {
	      // 	tmplayer.cluster[nj].erase(tmplayer.cluster[nj].begin()+jk);
	      // } else {jk++;}
	    }	// for(int jk=0;jk<int(tmplayer.cluster[nj].size())-1;jk++) {
	  } // for(int nj=0;nj<nside;nj++) {
	  TimeGroup[isj].alllayer[nl] = tmplayer;
	} // for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	
#ifdef isDebug
	cout << " Cluster Pos " << isj << endl;
	for(int nj=0;nj<nside;nj++) {
	  cout << " " << sideMark[nj] << endl;
	  for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	    tmplayer = TimeGroup[isj].alllayer[nl];
	    cout << "  l" << tmplayer.zlay;
	    for(int jk=0;jk<int(tmplayer.cluster[nj].size());jk++) {
	      cout << " (";
	      for(int ki=0;ki<int(tmplayer.cluster[nj][jk].size());ki++) {
		cout << " " << tmplayer.cluster[nj][jk][ki].strp;
	      }
	      cout << " )";
	    }
	    cout << endl;
	  }
	} // for(int nj=0;nj<nside;nj++) {
#endif	// #ifdef isDebug
	
// #ifdef isSpclDebug
// 	// cout << " TimeGroup Pos " << isj << endl;
// 	for(int nj=0;nj<nside;nj++) {
// 	  // cout << " " << sideMark[nj] << endl;
// 	  for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
// 	    tmplayer = TimeGroup[isj].alllayer[nl];
// 	    // cout << "  l" << tmplayer.zlay;
// 	    for(int jk=0;jk<int(tmplayer.cluster[nj].size());jk++) {
// 	      if(int(tmplayer.cluster[nj][jk].size())<5) {continue;}
// 	      cout << " (";
// 	      for(int ki=0;ki<int(tmplayer.cluster[nj][jk].size());ki++) {
// 		cout << " iev " << iev
// 		     << " " << sideMark[nj]
// 		     << " l" << tmplayer.zlay
// 		     << " " << tmplayer.cluster[nj][jk][ki].strp
// 		     << endl;
// 	      }
// 	      cout << " )" << endl;
// 	    }
// 	    // cout << endl;
// 	  }
// 	} // for(int nj=0;nj<nside;nj++) {
// #endif	// #ifdef isDebug
	
      }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {


      
      

      /*****  Forming the vector of PosInSpace *****/
      
      for(int isj=0;isj<int(TimeGroup.size());isj++) {
	TimeGroup[isj].PosForTrack.clear();
	for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	  // tmplayer = TimeGroup[isj].alllayer[nl];
	  tmplayer = TimeGroup[isj].alllayer[nl];
	  tmpPos.rawhitInfo = TimeGroup[isj].alllayer[nl];
	  
	  if(int(tmplayer.hit[0].size())>MaxTotHit ||
	     int(tmplayer.hit[1].size())>MaxTotHit ||
	     int(tmplayer.cluster[0].size())>MaxTotPos ||
	     int(tmplayer.cluster[1].size())>MaxTotPos) {continue;}
	  
	  for(int xn=0;xn<int(tmplayer.cluster[0].size());xn++) {
	    for(int yn=0;yn<int(tmplayer.cluster[1].size());yn++) {
	      // cout << " in " << isj << " " << tmplayer.zlay
	      // 	   << " " << xn << " " << yn << endl;
	      if(int(tmplayer.cluster[0][xn].size())>MaxMulti ||
		 int(tmplayer.cluster[1][yn].size())>MaxMulti) {continue;}
	      tmpPos.rawhitInfo.cluster[0].clear();
	      tmpPos.rawhitInfo.cluster[0].push_back(tmplayer.cluster[0][xn]);
	      tmpPos.rawhitInfo.cluster[1].clear();
	      tmpPos.rawhitInfo.cluster[1].push_back(tmplayer.cluster[1][yn]);	      
	      // cout << " out " << isj << " " << tmplayer.zlay
	      // 	   << " " << xn << " " << yn << endl;
	      // tmpPos.rawhitInfo.hit[0] = tmplayer.cluster[0][xn];
	      // tmpPos.rawhitInfo.hit[1] = tmplayer.cluster[1][yn];
	      calHitPos(tmpPos);
	      TimeGroup[isj].PosForTrack.push_back(tmpPos);
	      // cout << " size " << TimeGroup[isj].PosForTrack.size() << endl;
	    } // for(int yn=0;yn<int(tmplayer.cluster[1].size());yn++) {
	  } // for(int xn=0;xn<int(tmplayer.cluster[0].size());xn++) {
	} // for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	
#ifdef isDebug
	cout << " TimeGroup Cluster List " << isj << endl;
	for(int nj=0;nj<nside;nj++) {
	  cout << " " << sideMark[nj] << endl;
	  for(int nl=0;nl<int(TimeGroup[isj].PosForTrack.size());nl++) {
	    tmpPos = TimeGroup[isj].PosForTrack[nl];
	    cout << "  l" << tmpPos.rawhitInfo.zlay;
	    for(int ki=0;ki<int(tmpPos.rawhitInfo.cluster[nj].back().size());ki++) {
	      cout << " " << tmpPos.rawhitInfo.cluster[nj][0][ki].strp;
	      cout << " (";
	      for(int nht=0;nht<int(tmpPos.rawhitInfo.cluster[nj][0][ki].tdc[0].size());nht++) {
		cout << " " << tmpPos.rawhitInfo.cluster[nj][0][ki].tdc[0][nht] - tmpPos.rawhitInfo.tdc_ref[0];
	      }
	      cout << " )";
	    }
	    cout << endl;
	  }
	} // for(int nj=0;nj<nside;nj++) {
#endif	// #ifdef isDebug
	
      }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {




      
      
      
      /*****  Analysis and Stuff *****/
      
      for(int isj=0;isj<int(TimeGroup.size());isj++) {
	// cout << " TimeGroup " << isj << endl;
	
      }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {
      
      
    } // for(Long64_t iev=nentrymn;iev<nentry;iev++) {
    
    
  } // if(!fileIn->IsZombie()) {


  
  fileIn->Close();
  
  
  return 0;
}

