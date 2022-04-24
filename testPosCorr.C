
#include "mIcal_correction_20210325.txt"


void testPosCorr() {

  const int        nside         =   2;
  const int        nlayer        =  10;
  const int        nxrow         =   1;
  const int        nyrow         =   1;
  const int        nmodule       =   1;
  const int        nstrip        =  64;
  const int        nmxhits       =  6;

  const int        blockM        = nstrip/4;

  const double     maxPosDev     = 2.;

  const char *sideMark[nside] = {"x","y"};
  
  char name[300];
  
  double poserrsq1[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  double blockposoff1[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  
  TFile *f2 = new TFile("../rredata/temp/GMA_RPCv4t_evtraw-20181224_192913_corr_iter5.root","read");
  
  TH1D *layPosR[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  TH1D *blockPosR[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      sprintf(name,"layPosR_m%i_xr%i_yr%i_%s%i_mul%i",
		      nm,nx,ny,sideMark[nj],nl,nh+1);
	      layPosR[nm][nx][ny][nl][nj][nh] = (TH1D *)f2->Get(name);
	      int status = layPosR[nm][nx][ny][nl][nj][nh]->Fit("gaus","SQ");
	      if(status==0) {
		TF1 *f1 = (TF1*)layPosR[nm][nx][ny][nl][nj][nh]->GetFunction("gaus");
		poserrsq1[nm][nx][ny][nl][nj][nh] = pow(f1->GetParameter(2),2.);
	      } else {poserrsq1[nm][nx][ny][nl][nj][nh] = -1000;}
	      // cout << poserrsq[nm][nx][ny][nl][nj][nh]
	      // 	   << ",\t// m" << nmod
	      // 	   << " xr" << nx << " yr" << ny << " z" << nl
	      // 	   << " " << sideMark[nj] << " mul" << nh+1
	      // 	   << endl;
	    }
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		sprintf(name,"blockPosR_m%i_xr%i_yr%i_%s%i_%i_%i",
			nm,nx,ny,sideMark[nj],nl,nsx,nsy);
		blockPosR[nm][nx][ny][nl][nj][nsx][nsy] = (TH1D *)f2->Get(name);
		int status = blockPosR[nm][nx][ny][nl][nj][nsx][nsy]->Fit("gaus","SQ");
	      if(status==0) {
		TF1 *f1 = (TF1*)blockPosR[nm][nx][ny][nl][nj][nsx][nsy]->GetFunction("gaus");
		blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] = f1->GetParameter(1);
	      } else {blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] = 0.;}
	      // if(fabs(blockposoff1[nm][nx][ny][nl][nj][nsx][nsy])>maxPosDev) {blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] = 0.;}
	      // cout << blockposoff[nm][nx][ny][nl][nj][nsx][nsy]
	      // 	   << ",\t// m" << nmod
	      // 	   << " xr" << nx << " yr" << ny << " z" << nl
	      // 	   << " " << sideMark[nj]
	      // 	   << " blk " << nsx << " " << nsy
	      // 	   << endl;
	    }}}}}}}

  ofstream file1("testTime.txt");
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      file1 << "   " << poserrsq1[nm][nx][ny][nl][nj][nh]
		    << ",\t// m" << nm
		    << " xr" << nx << " yr" << ny << " z" << nl
		    << " " << sideMark[nj] << " mul" << nh+1
		    << endl;
	      
	    }}}}}}
  file1 << endl;
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		double outval = blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] + blockposoff[nm][nx][ny][nl][nj][nsx][nsy];
		if(fabs(outval)>maxPosDev) {outval = 0.;}
		file1 << "   "
		      << outval
		      << ",\t// m" << nm
		      << " xr" << nx << " yr" << ny << " z" << nl
		      << " " << sideMark[nj]
		      << " blk " << nsx << " " << nsy
		      << endl;
	      }}}}}}}
  

}



void testTimeCorr() {

  const int        nside         =   2;
  const int        nlayer        =  10;
  const int        nxrow         =   1;
  const int        nyrow         =   1;
  const int        nmodule       =   1;
  const int        nstrip        =  64;
  const int        nmxhits       =  6;

  const int        blockM        = nstrip/4;
  
  const double     maxPosDev     = 20.;

  const char *sideMark[nside] = {"x","y"};
  
  char name[300];
  
  double mulposoff1[nside][nmxhits];
  double poserrsq1[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  double blockposoff1[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  
  TFile *f2 = new TFile("../rredata/temp/GMA_RPCv4t_evtraw-20181224_192913_iter12t.root","read");
  
  TH1D *timeR[nside][nmxhits];
  for(int nj=0;nj<nside;nj++) {
    for(int nh=0;nh<nmxhits;nh++) {
      sprintf(name,"layTimeR_m0_xr0_yr0_x0_mul%i",nh+1);
      timeR[nj][nh] = (TH1D *)((TH1D *)f2->Get(name))->Clone();;
      timeR[nj][nh]->Scale(0);
    }}
  
  TH1D *layTimeR[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  TH1D *blockTimeR[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      sprintf(name,"layTimeR_m%i_xr%i_yr%i_%s%i_mul%i",
		      nm,nx,ny,sideMark[nj],nl,nh+1);
	      layTimeR[nm][nx][ny][nl][nj][nh] = (TH1D *)f2->Get(name);
	      timeR[nj][nh]->Add(layTimeR[nm][nx][ny][nl][nj][nh]);
	      int status = layTimeR[nm][nx][ny][nl][nj][nh]->Fit("gaus","SQ");
	      if(status==0) {
		TF1 *f1 = (TF1*)layTimeR[nm][nx][ny][nl][nj][nh]->GetFunction("gaus");
		poserrsq1[nm][nx][ny][nl][nj][nh] = pow(f1->GetParameter(2),2.);
	      } else {poserrsq1[nm][nx][ny][nl][nj][nh] = -1000;}
	      // cout << poserrsq[nm][nx][ny][nl][nj][nh]
	      // 	   << ",\t// m" << nmod
	      // 	   << " xr" << nx << " yr" << ny << " z" << nl
	      // 	   << " " << sideMark[nj] << " mul" << nh+1
	      // 	   << endl;
	    }
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		sprintf(name,"blockTimeR_m%i_xr%i_yr%i_%s%i_%i_%i",
			nm,nx,ny,sideMark[nj],nl,nsx,nsy);
		blockTimeR[nm][nx][ny][nl][nj][nsx][nsy] = (TH1D *)f2->Get(name);
		int status = blockTimeR[nm][nx][ny][nl][nj][nsx][nsy]->Fit("gaus","SQ");
	      if(status==0) {
		TF1 *f1 = (TF1*)blockTimeR[nm][nx][ny][nl][nj][nsx][nsy]->GetFunction("gaus");
		blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] = f1->GetParameter(1);
	      } else {blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] = 0.;}
	      // cout << blockposoff[nm][nx][ny][nl][nj][nsx][nsy]
	      // 	   << ",\t// m" << nmod
	      // 	   << " xr" << nx << " yr" << ny << " z" << nl
	      // 	   << " " << sideMark[nj]
	      // 	   << " blk " << nsx << " " << nsy
	      // 	   << endl;
	    }}}}}}}

  for(int nj=0;nj<nside;nj++) {
    for(int nh=0;nh<nmxhits;nh++) {
      int status = timeR[nj][nh]->Fit("gaus","SQ");
      if(status==0) {
	TF1 *f1 = (TF1*)timeR[nj][nh]->GetFunction("gaus");
	mulposoff1[nj][nh] = f1->GetParameter(1);
      } else {mulposoff1[nj][nh] = 0.;}
    }}
  
  ofstream file1("testTime.txt");
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      file1 << "   " << poserrsq1[nm][nx][ny][nl][nj][nh]
		    << ",\t// m" << nm
		    << " xr" << nx << " yr" << ny << " z" << nl
		    << " " << sideMark[nj] << " mul" << nh+1
		    << endl;
	      
	    }}}}}}
  file1 << endl;
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		double outval = blockposoff1[nm][nx][ny][nl][nj][nsx][nsy] + blocktimeoff[nm][nx][ny][nl][nj][nsx][nsy];
		if(fabs(outval)>maxPosDev) {outval = 0.;}
		file1 << "   "
		      << outval
		      << ",\t// m" << nm
		      << " xr" << nx << " yr" << ny << " z" << nl
		      << " " << sideMark[nj]
		      << " blk " << nsx << " " << nsy
		      << endl;
	      }}}}}}}
  file1 << endl;
  for(int nj=0;nj<nside;nj++) {
    for(int nh=0;nh<nmxhits;nh++) {
      file1 << "   " << mulposoff1[nj][nh]
	    << ",\t// " << sideMark[nj]
	    << " mul" << nh+1
	    << endl;
    }}
  

}
