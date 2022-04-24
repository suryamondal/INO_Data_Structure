
void testTDC() {

  const int        nside         =   2;
  const int        nlayer        =  10;
  const int        nxrow         =   1;
  const int        nyrow         =   1;
  const int        nmodule       =   1;
  const int        nstrip        =  64;

  const char *sideMark[nside] = {"x","y"};

  char name[300];

  // TFile *f2 = new TFile("../rredata/GMA_RPCv4t_evtraw-20181217_112935_corr.root","read");
  TFile *f2 = new TFile("../rredata/temp/GMA_RPCv4t_evtraw-20181224_192913_corr_iter3.root","read");
  
  double TDC_Shift[nside][nlayer][nstrip];
  double lowValTDCg[nside][nlayer];// = 99999;
  TH1D *rawtdcStrp[nside][nlayer][nstrip];
  // for(int nj=0;nj<nside;nj++) {
  //   for(int nl=0;nl<nlayer;nl++) {
  //     for(int ns=0;ns<nstrip;ns++) {
  // 	sprintf(name,"rawtdcStrp_%s%i_s%i",sideMark[nj],nl,ns);
  // 	rawtdcStrp[nj][nl][ns] = new TH1D(name,name,240,80,140);
  //     }}}

  f2->cd();
  for(int nj=0;nj<nside;nj++) {
    for(int nl=0;nl<nlayer;nl++) {
      lowValTDCg[nj][nl] = 99999;
      for(int ns=0;ns<nstrip;ns++) {
  	sprintf(name,"rawtdcStrp_%s%i_s%i",sideMark[nj],nl,ns);
  	rawtdcStrp[nj][nl][ns] = (TH1D *)f2->Get(name);
	TDC_Shift[nj][nl][ns] = rawtdcStrp[nj][nl][ns]->GetMean();
	if(TDC_Shift[nj][nl][ns]!=0 &&
	   lowValTDCg[nj][nl]>TDC_Shift[nj][nl][ns]) {
	  lowValTDCg[nj][nl] = TDC_Shift[nj][nl][ns];}
	cout << " " << nj << " " << nl << " " << ns
	     << " " << TDC_Shift[nj][nl][ns] << endl;
      }
      cout << " lowValTDCg " << lowValTDCg[nj][nl] << endl;
    }}

  cout << endl;
  for(int nmod=0;nmod<nmodule;nmod++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  for(int nj=0;nj<nside;nj++) {
	    
	    for(int ns=0;ns<nstrip;ns++) {
	      if(TDC_Shift[nj][nl][ns]!=0) {
		TDC_Shift[nj][nl][ns] -= lowValTDCg[nj][nl];
	      } else {TDC_Shift[nj][nl][ns] = -100.;}
	      
	      // cout << " " << nj << " " << nl << " " << ns
	      //      << " " << TDC_Shift[nj][nl][ns] << endl;
	      cout << "   " << TDC_Shift[nj][nl][ns]
		   << ",\t// m" << nmod
		   << " xr" << nx << " yr" << ny << " z" << nl
		   << " " << sideMark[nj] << " str" << ns
		   << endl;
	    }}}}}}
  


}
