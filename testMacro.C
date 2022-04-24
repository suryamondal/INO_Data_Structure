void extTest() {

  int count = 100000;
  double sigma1 = 2.;
  double sigma2 = 1.;

  double sigma3 = sqrt(sigma1*sigma1-sigma2*sigma2);
  cout << " sigma3 " << sigma3 << endl; 
  
  TH1D *h1 = new TH1D("h1","h1",100,-10,10);
  
  for(int ij=0;ij<count;ij++) {
    h1->Fill(gRandom->Gaus(0,sigma1),1./pow(gRandom->Gaus(0,sigma2),2.));
  }

  // h1->Draw("hist");
  h1->Fit("gaus");

}
