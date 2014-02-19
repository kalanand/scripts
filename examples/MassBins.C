#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <cmath>
void MassBins()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(000); 
  gStyle->SetPalette(1);
  double m,x,dm,center,r,sigma;
  double min = 1;
  int Nbins,i;
  double Boundaries[200],v1[200],v2[200],pull[200];
  m = min;
  Nbins = 0;
  while (m<14000)
    {
      x = m;
      Boundaries[Nbins] = x;
      //m = 1.*TMath::Nint(2*Solution(x)-x);
      m = 2*Solution(x)-x;
      //if (m-x>100)
        //m = 1.*TMath::Nint(m);
      //else 	
        m = 0.1.*TMath::Nint(m*10);
      dm = m-x;	
      center = 0.5*(m+x);
      r = Resolution(center);
      sigma = r*center;
      v1[Nbins] = dm/center;
      v2[Nbins] = r;
      pull[Nbins] = (dm-sigma)/sigma;
      cout<<x<<" "<<m<<" "<<center<<" "<<dm<<" "<<sigma<<" "<<dm/center<<" "<<r<<endl;
      Nbins++;
    } 
   Boundaries[Nbins] = 14000;
   TH1F *h1 = new TH1F("h1","h1",Nbins,Boundaries);
   TH1F *h2 = new TH1F("h2","h2",Nbins,Boundaries);
   TH1F *h3 = new TH1F("h3","h3",Nbins,Boundaries);
   for(i=0;i<Nbins;i++)
     {
       h1->SetBinContent(i+1,v1[i]);
       h1->SetBinError(i+1,0); 
       h2->SetBinContent(i+1,v2[i]);
       h2->SetBinError(i+1,0);
       h3->SetBinContent(i+1,100*pull[i]);
       h3->SetBinError(i+1,0);
     }  
    ofstream MassFile;
    MassFile.open("DefaultMassBins.h");
    MassFile << "#ifndef DEFAULTMASSBINS_H\n";
    MassFile << "#define DEFAULTMASSBINS_H\n";
    MassFile << "namespace\n";
    MassFile << "{\n";
    MassFile << "  const int nMassBins = "<<Nbins-1<<";\n";
    MassFile << "  double massBoundaries[nMassBins+1] = {";
    for(i=0;i<Nbins;i++)
      MassFile << Boundaries[i] << ",";
    MassFile << Boundaries[i] << "};\n";
    MassFile << "}\n";
    MassFile << "#endif"; 
     
   TCanvas *c1 = new TCanvas("resolution","resolution",900,600);  
   gPad->SetLogx();
   h1->SetMarkerStyle(20);
   h2->SetMarkerColor(2);
   h2->SetMarkerStyle(21);
   h2->SetLineColor(2);
   h1->SetTitle("");
   h1->GetXaxis()->SetRangeUser(200,10000);
   h1->GetXaxis()->SetTitle("M_{jj} (GeV)");
   h1->GetYaxis()->SetTitle("Relative Resolution");
   TLegend *leg = new TLegend(0.7,0.7,0.85,0.85);
   leg->AddEntry(h1,"#Deltam/m","L");
   leg->AddEntry(h2,"#sigma(m)/m","L");
   leg->SetLineColor(0);
   leg->SetFillColor(0);
   h1->Draw();
   h2->Draw("same");   
   leg->Draw();
   
   TCanvas *c2 = new TCanvas("pull","pull",900,600); 
   gPad->SetLogx();
   h3->Draw();
   h3->SetTitle("");
   h3->GetXaxis()->SetRangeUser(200,10000);
   h3->GetXaxis()->SetTitle("M_{jj} (GeV)");
   h3->GetYaxis()->SetTitle("(#Deltam-#sigma)/#sigma (%)");
}

double func(double u, double x)
{
  double A = 0.01686;
  double B = 2.598;
  double p = 0.6304;
  double result = (A-2.)*u + B*pow(u,1.-p) + 2*x; 
  return result; 
}

double funcDer(double u)
{
  double A = 0.01686;
  double B = 2.598;
  double p = 0.6304;
  double result = (A-2.) + (1.-p)*B*pow(u,-p); 
  return result; 
}

double Solution(double x)
{
  int i;
  double tmp,result,e;
  tmp = x;
  e = 100;
  for(i=0;i<10;i++)
    {
      result = tmp-func(tmp,x)/funcDer(tmp);
      e = fabs(result-tmp);
      tmp = result;
      if (e<1e-4) continue;
    }
  return result;  
}

double Resolution(double x)
{
  double A = 0.01686;
  double B = 2.598;
  double p = 0.6304;
  return A+B/pow(x,p);
}
