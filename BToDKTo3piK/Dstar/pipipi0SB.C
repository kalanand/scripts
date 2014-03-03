#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>

// plots the signal and sideband regions along with backgroung 
// contributions from different types of backgroungs using MC truth 
//  information.


void pipipi0SB() {
  static Int_t nBins = 108;
//   static Int_t nBins = 216;

  TH1D hsig("hsig","",nBins,1.72,1.99);
  TH1D hbkg("hbkg","",nBins,1.72,1.99);
  TH1D href("href","",nBins,1.72,1.99);

  TH1D hb1("hb1","",nBins,1.72,1.99);
  TH1D hb2("hb2","",nBins,1.72,1.99);
  TH1D hb3("hb3","",nBins,1.72,1.99);


  Double_t weight1= 179.018*13.0;
  Double_t weight2= 332.850*5.5;
  Double_t weight3= 215.286*5.5;
  Double_t weight4= 162.902*21.0;
  Double_t sumweight = weight1+weight2+weight3+weight4;
  weight1 /= sumweight;
  weight2 /= sumweight;
  weight3 /= sumweight;
  weight4 /= sumweight;

  TCut chaincut("abs(delta_M-0.1455)<0.0006 && Dmass>1.72 &&  Dmass<1.99");
  TCut bgCut((chaincut)&&"!(Flag==1||Flag==2||Flag==11||Flag==101) && (hh_mass<0.489 || hh_mass>0.508)");
  TCut refCut((chaincut)&&"(Flag==10 || Flag==20) && abs(lundmcD0)>1 && Dmass<1.826");
  TCut norefCut((chaincut)&&"!((Flag==10 || Flag==20) && abs(lundmcD0)>1)");

  TChain chain1("ntp1");
  chain1.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/ccbar*.root");
  chain1.SetWeight(weight1,"global");
  TChain chain2("ntp1");
  chain2.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bbbar*.root");
  chain2.SetWeight(weight2,"global");
  TChain chain3("ntp1");
  chain3.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bplusbminus*.root");
  chain3.SetWeight(weight3,"global");
  TChain chain4("ntp1");
  chain4.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/uds*.root");
  chain4.SetWeight(weight4,"global");

  TChain chain("chain");
  chain.Add( &chain1 );
  chain.Add( &chain2 );
  chain.Add( &chain3 );
  chain.Add( &chain4 );
  TTree* myTree = (TTree*) &chain;



  myTree->Draw("Dmass>>hsig",(chaincut)&&"(hh_mass<0.489 || hh_mass>0.508)", "goff"); 
  myTree->Draw("Dmass>>hbkg",bgCut, "goff"); 
  myTree->Draw("Dmass>>href",refCut, "goff"); 
  chain2.Draw("Dmass>>hb2",bgCut, "goff"); 
  chain3.Draw("Dmass>>+hb2",bgCut, "goff"); 
  chain4.Draw("Dmass>>hb3",bgCut, "goff");

  hb2.SetLineColor(3);
  hb3.SetLineColor(0);
  hb1.SetMinimum(0);
  hb2.SetMinimum(0);
  hb3.SetMinimum(0);
  hb1.Scale(2);
  hb2.Scale(2);
  hb3.Scale(2);

  hsig.SetFillColor(2);
  hbkg.SetFillColor(4);
  href.SetFillColor(6);

  hsig.SetMinimum(0);
  hbkg.SetMinimum(0);
  href.SetMinimum(0);

  TCanvas canvas("canvas", "canvas", 880, 680); 
  gStyle->SetOptStat(0);
  hsig.Draw();  
  hbkg.Draw("same");
  href.Draw("same");
  //  TLine* sigline1 = new TLine(1.8459, 50, 1.8459, 3000);
 TLine* sigline1 = new TLine(1.82, 50, 1.82, 1200);
  sigline1->Draw();
  //  TLine* sigline2 = new TLine(1.8739, 50, 1.8739, 3000);
 TLine* sigline2 = new TLine(1.90, 50, 1.90, 1200);
  sigline2->Draw();
  TLine* sideline1 = new TLine(1.93, 50, 1.93, 1000);
  sideline1->Draw();
  TLine* sideline2 = new TLine(1.985, 50, 1.985, 1000);
  sideline2->Draw();
  leg_hist = new TLegend(0.12,0.7,0.4,0.89);
  leg_hist->SetHeader("#pi^{-}#pi^{+}#pi^{0} mass plot for generic MC");
  leg_hist->AddEntry(&hsig,"signal","f");
  leg_hist->AddEntry(&hbkg,"comb. background","f");
  leg_hist->AddEntry(&href,"K^{-}#pi^{+}#pi^{0} events misreconstructed","f");
  leg_hist->Draw();
  canvas.SaveAs("generic_pipipi0-1.eps");
  canvas.SaveAs("generic_pipipi0-1.gif");
  canvas.Close();


  TCanvas canvas("canvas", "canvas", 880, 680); 
  gStyle->SetOptStat(0);
  hbkg.Draw();
  href.Draw("same");
  hsig.Draw("same"); 
  hbkg.Draw("same");
  href.Draw("same");
  hb2.Draw("same");
  hb3.Draw("same");
  //  TLine* sigline1 = new TLine(1.8459, 0, 1.8459, 600);
 TLine* sigline1 = new TLine(1.82, 0, 1.82, 400);
  sigline1->Draw();
  //  TLine* sigline2 = new TLine(1.8739, 0, 1.8739, 600);
 TLine* sigline2 = new TLine(1.90, 0, 1.90, 400);
  sigline2->Draw();
  TLine* sideline1 = new TLine(1.93, 0, 1.93, 400);
  sideline1->Draw();
  TLine* sideline2 = new TLine(1.985, 0, 1.985, 400);
  sideline2->Draw();
  leg_hist = new TLegend(0.7,0.7,0.96,0.9);
  leg_hist->SetHeader("#pi^{-}#pi^{+}#pi^{0} mass plot for generic MC");
  leg_hist->AddEntry(&hsig,"signal","f");
  leg_hist->AddEntry(&href,"K^{-}#pi^{+}#pi^{0} events misreconstructed","f");
  leg_hist->AddEntry(&hbkg,"all comb. background","f");
  leg_hist->AddEntry(&hb2,"b #bar{b}","l");
  leg_hist->AddEntry(&hb3,"u#bar{u}, d#bar{d}, s#bar{s}","l");
  leg_hist->Draw();
  canvas.SaveAs("generic_pipipi0-2.eps");
  canvas.SaveAs("generic_pipipi0-2.gif");

}
