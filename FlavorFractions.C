{

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadRightMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  gStyle->SetOptFit(1111);
  tdrStyle->SetOptStat(0); 


  int nBins = 50;
  double min = 20.0;
  double max = 520.0;

  TChain chain("ZJet");
  chain.Add("Summer09-7TeV-ZeeJets-Pt_*.root");



  TH1D* flavAll = new TH1D("flavAll","", nBins, 20,520);
  chain.Draw("JetCorPt[2][0]>>flavAll","mZee>20 && abs(JetCorDphi[2][0])>2.5","goff");

  TH1D* flavGluon = new TH1D("flavGluon","", nBins, 20,520);
  chain.Draw("JetCorPt[2][0]>>flavGluon",
	     "mZee>20 && abs(JetCorDphi[2][0])>2.5 && abs(JetCorFlavor[2][0])==21","goff");

  TH1D* flavLF = new TH1D("flavLF","", nBins, 20,520);
  chain.Draw("JetCorPt[2][0]>>flavLF",
	     "mZee>20 && abs(JetCorDphi[2][0])>2.5 && abs(JetCorFlavor[2][0])<4","goff");


  TH1D* flavHF = new TH1D("flavHF","", nBins, 20,520);
  chain.Draw("JetCorPt[2][0]>>flavHF","mZee>20 && abs(JetCorDphi[2][0])>2.5 && (abs(JetCorFlavor[2][0])==4 || abs(JetCorFlavor[2][0])==5)","goff");



  flavGluon->Divide(flavAll);
  flavLF->Divide(flavAll);
  flavHF->Divide(flavAll);


  flavGluon->SetLineColor(6);
  flavLF->SetLineColor(4);
  flavHF->SetLineColor(2);
 

  flavGluon->SetLineWidth(4);
  flavLF->SetLineWidth(3);
  flavHF->SetLineWidth(4);


  flavGluon->SetLineStyle(4);
  flavHF->SetLineStyle(2);


  flavGluon->Smooth(10);
  flavLF->Smooth(10);
  flavLF->GetXaxis()->SetTitle("Corrected jet p_{T} (GeV)");
  flavLF->GetXaxis()->SetNdivisions(505);
  flavLF->GetYaxis()->SetTitle("Fraction of jets");
  flavLF->GetYaxis()->SetRangeUser(0.01, 1.2);




  TCanvas* c2 = new TCanvas("c2","Jet flavor", 500, 500);
  flavLF->Draw();
  flavGluon->Draw("same");
  flavHF->Draw("same");
  gPad->SetLogy(1);
  TLatex* CMS = new  TLatex();
  CMS->SetTextAlign(12);
  CMS->SetTextSize(0.04);
  CMS->SetNDC();
  CMS->DrawLatex(0.2, 0.2, "CMS Preliminary");
  CMS->DrawLatex(0.8, 0.8, "#color[4]{uds}");
  CMS->DrawLatex(0.8, 0.65, "#color[6]{gluon}");
  CMS->DrawLatex(0.8, 0.45, "#color[2]{b,c}");
  c2->Update();
  c2->SaveAs("flavor_composition_Zjets_7TeV.eps");
  c2->SaveAs("flavor_composition_Zjets_7TeV.gif");
  c2->SaveAs("flavor_composition_Zjets_7TeV.png");
  c2->SaveAs("flavor_composition_Zjets_7TeV.root");
}
