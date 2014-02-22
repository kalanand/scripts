#include <vector>
#include "TMatrixD.h"
#include "CondFormats/JetMETObjects/src/CorrelationGenerator.cc"
void Test()
{
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(0000);
  gStyle->SetOptFit(000); 
  gStyle->SetPalette(1);
  
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();  


  double coeff[4] = {4,1.5,1.5,9};
  TArrayD data(4);
  TMatrixD A(2,2);
  for(int i=0;i<4;i++)
    data[i] = coeff[i];
  A.Use(2,2,data.GetArray()); 
  vector<double> mY;
  mY.push_back(4);
  mY.push_back(7);
  CorrelationGenerator *Cor = new CorrelationGenerator(mY,A);
  
  
  TH2F *hY = new TH2F("hY","hY",500,-50,50,500,-50,50);
  vector<double> v;
  v = Cor->getOrthSigmas();
  cout<<v[0]<<" "<<v[1]<<endl;
  v = Cor->getOrthMeans();
  cout<<v[0]<<" "<<v[1]<<endl;
  v = Cor->getEigenVector(0);
  cout<<v[0]<<" "<<v[1]<<endl;
  v = Cor->getEigenVector(1);
  cout<<v[0]<<" "<<v[1]<<endl;
  
  for(int i=0;i<20000;i++)
    {
      v = Cor->getRandomSet();
      hY->Fill(v[0],v[1]);
      
    }  
    
    
 
  TH1D *hY1 = (TH1D*)hY->ProjectionX("hY1"); 
  TH1D *hY2 = (TH1D*)hY->ProjectionY("hY2");  
    
  cout<<hY->GetCorrelationFactor()<<" "<<hY->GetCovariance()<<endl;
    
    
  TCanvas *c = new TCanvas("can","can",600,600);
  
  hY->Draw("colz");  
  
  TCanvas *c1 = new TCanvas("can1","can1",900,600);
  c1->Divide(2,1);
  c1->cd(1);
  hY1->Draw();
  c1->cd(2);
  hY2->Draw();
  
}

