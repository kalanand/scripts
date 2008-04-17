// $Id: dalitzVarCorr.cc,v 1.2 2006/04/25 17:51:47 fwinkl Exp $
// Checks correlations of the Dalitz variables with other (fit) variables

TCanvas *can = 0;

// Do all the event types
void dalitzVarCorr()
{
  gROOT->cd();
  TCut cut = TCut(cutBasic+cutDeltaE);

  ofstream of;
  of.open("dalitzVarCorr.txt");

  TString s;

  s = checkCorrAll(qqTree,cutqqBadD+cut);
  cout << s; of << s;
  s = checkCorrAll(qqTree,cutqqGoodD+cut);
  cout << s; of << s;
  s = checkCorrAll(bbTree,cutBBBadD+cut);
  cout << s; of << s;
  s = checkCorrAll(bbTree,cutDPiX+cut);
  cout << s; of << s;
  s = checkCorrAll(bbTree,cutDKX+cut);
  cout << s; of << s;
  s = checkCorrAll(dpiTree,cutDPiGoodD+cut);
  cout << s; of << s;
  of.close();
}

// Do all the different binning variables
TString checkCorrAll(TTree *tree, TCut cut)
{
   TVector3 c12, c13, c23;
   TString s = "% ";
   s += cut.GetName();
   s += "\n";
   s += checkCorr(tree,cut,"mes",5.24,5.26);
   s += checkCorr(tree,cut,"Deltae",-0.025,0.025);
   s += checkCorr(tree,cut,"nnout",0.2,0.7);
   s += checkCorr(tree,cut,"d0mass",1.845,1.875);
   return s;
}

// Do test for one event type and one binning variable
TString checkCorr(TTree *tree, TCut addCut, const char* var, Double_t r1, Double_t r2)
{
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  gStyle->SetOptStat(1001110);


  TString s;
  s.Form("%s<%.4g",var,r1);
  TCut cut1 = TCut(s);

  s.Form("%s>%.4g && %s<%.4g",var,r1,var,r2);
  TCut cut2 = TCut(s);

  s.Form("%s>%.4g",var,r2);
  TCut cut3 = TCut(s);
  
  const int BINS = 30;
  TString type = TString(addCut.GetName())+": ";
  TH2D *h1 = new TH2D("h1",type+cut1.GetTitle(),BINS,0,3,BINS,0,3);
  TH2D *h2 = new TH2D("h2",type+cut2.GetTitle(),BINS,0,3,BINS,0,3);  
  TH2D *h3 = new TH2D("h3",type+cut3.GetTitle(),BINS,0,3,BINS,0,3);

  h1->GetXaxis()->SetTitle(m12->GetTitle());
  h1->GetYaxis()->SetTitle(m13->GetTitle());
  h2->GetXaxis()->SetTitle(m12->GetTitle());
  h2->GetYaxis()->SetTitle(m13->GetTitle());
  h3->GetXaxis()->SetTitle(m12->GetTitle());
  h3->GetYaxis()->SetTitle(m13->GetTitle());
  h1->Sumw2();
  h2->Sumw2(); 
  h3->Sumw2(); 

  tree->Draw("d0ppmupmass**2:d0pppupmass**2>>h1",addCut+cut1,"goff");
  tree->Draw("d0ppmupmass**2:d0pppupmass**2>>h2",addCut+cut2,"goff");
  tree->Draw("d0ppmupmass**2:d0pppupmass**2>>h3",addCut+cut3,"goff");

  // KS test
  Double_t ks12 = h1->KolmogorovTest(h2);
  Double_t ks13 = h1->KolmogorovTest(h3);
  Double_t ks23 = h2->KolmogorovTest(h3);

  /*
  // Chi2 test
  TH2D *h12chi2 = new TH2D("h12chi2","Bin1,Bin2",BINS,0,3,BINS,0,3);
  TH2D *h13chi2 = new TH2D("h13chi2","Bin1,Bin3",BINS,0,3,BINS,0,3);
  TH2D *h23chi2 = new TH2D("h23chi2","Bin2,Bin3",BINS,0,3,BINS,0,3);
  h1->Scale(1000/h1->Integral());
  h2->Scale(1000/h2->Integral());
  h3->Scale(1000/h3->Integral());
  TVector3 r12chi2 = chi2test2d(h1,h2,h12chi2);
  TVector3 r13chi2 = chi2test2d(h1,h3,h13chi2);
  TVector3 r23chi2 = chi2test2d(h2,h3,h23chi2);
  */

  if (!can) { 
    can = new TCanvas("can","dalitzVarCorr",800,300);
    can->Divide(3,1);
  }

  can->cd(1);
  h1->Draw("colz");
  can->cd(2);
  h2->Draw("colz");
  can->cd(3);
  h3->Draw("colz");

  /*
  can->cd(4);
  h12chi2->Draw("colz");
  can->cd(5);
  h13chi2->Draw("colz");
  can->cd(6);
  h23chi2->Draw("colz");
  */

  TString filename = TString("dalitzVarCorr_")+TString(addCut.GetName())+"_"+TString(var);
  can->SaveAs(filename+".eps");  
  can->SaveAs(filename+".root");

  TString str = "";
  //str.Form("%8s %8s %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",addCut.GetName(),var,r12chi2[2],r13chi2[2],r23chi2[2],ks12,ks13,ks23);

  TString varTex = "\\"+TString(var);
  str.Form("\\%-8s & %8.5f & %8.5f & %8.5f \\\\ \n",var,ks12,ks13,ks23);

  return str;
}
