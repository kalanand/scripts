
const int nHists = 30;
char* histNames[nHists] = {"ZPt_Gen","ZPt_Reco","GenJetPt","CaloJetPt",
			   "PFJetPt","emPt_Gen","emPt_Reco","epPt_Gen",
			   "epPt_Reco","Zmass","ZEta_Gen","ZEta_Reco",
			   "GenJetEta","CaloJetEta","PFJetEta","GenJetPt2",
			   "CaloJetPt2","PFJetPt2","GenJetPt3","CaloJetPt3",
			   "PFJetPt3","CaloJetPtRatio2over1",
			   "CaloJetPtRatio3over2","PFJetPtRatio2over1",
			   "PFJetPtRatio3over2","GenJetPtRatio2over1",
			   "GenJetPtRatio3over2","nJetsCalo","nJetsPF",
			   "nJetsGen"};


void plotRatio() {


  // set proper style for plots
  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.18);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.16);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);


  TFile* f1 = new TFile("histograms_7TeV.root","read");
  TFile* f2 = new TFile("histograms_10TeV.root","read");
  gROOT->cd();
  TH1D ratioHist[nHists];
  for (int i=0; i<nHists; i++) {
    TH1D* numHist = (TH1D* ) f1->Get( histNames[i] );
    TH1D* denHist = (TH1D* ) f2->Get( histNames[i] );
    ratioHist[i] = (TH1D*) numHist->Clone(histNames[i]);
    ratioHist[i]->Divide( denHist );
    ratioHist[i]->GetYaxis()->SetTitle("Events @ 7 TeV / @ 10 TeV");
  }

  TFile* file = new TFile( "histograms_ratio.root", "RECREATE" );
  for (int i=0; i<nHists; i++) {
    file->WriteObject(ratioHist[i], histNames[i]);
  }

  TH1D* ZPt_Gen = (TH1D* ) file->Get("ZPt_Gen");
  TH1D* ZPt_Reco = (TH1D* ) file->Get("ZPt_Reco");
  TH1D* GenJetPt = (TH1D* ) file->Get("GenJetPt");
  TH1D* CaloJetPt = (TH1D* ) file->Get("CaloJetPt");
  TH1D* PFJetPt = (TH1D* ) file->Get("PFJetPt");

  TH1D* GenJetEta = (TH1D* ) file->Get("GenJetEta");
  TH1D* CaloJetEta = (TH1D* ) file->Get("CaloJetEta");
  TH1D* PFJetEta = (TH1D* ) file->Get("PFJetEta");
  
  TH1D* GenJetPt2 = (TH1D* ) file->Get("GenJetPt2");
  TH1D* CaloJetPt2 = (TH1D* ) file->Get("CaloJetPt2");
  TH1D* PFJetPt2 = (TH1D* ) file->Get("PFJetPt2");

  TH1D* GenJetPt3 = (TH1D* ) file->Get("GenJetPt3");
  TH1D* CaloJetPt3 = (TH1D* ) file->Get("CaloJetPt3");
  TH1D* PFJetPt3 = (TH1D* ) file->Get("PFJetPt3");

  TH1D* CaloJetPtRatio2over1 = (TH1D* ) file->Get("CaloJetPtRatio2over1");
  TH1D* CaloJetPtRatio3over2 = (TH1D* ) file->Get("CaloJetPtRatio3over2");
  TH1D* PFJetPtRatio2over1 = (TH1D* ) file->Get("PFJetPtRatio2over1");
  TH1D* PFJetPtRatio3over2 = (TH1D* ) file->Get("PFJetPtRatio3over2");
  TH1D* GenJetPtRatio2over1 = (TH1D* ) file->Get("GenJetPtRatio2over1");
  TH1D* GenJetPtRatio3over2 = (TH1D* ) file->Get("GenJetPtRatio3over2");


  TH1D* nJetsCalo = (TH1D* ) file->Get("nJetsCalo");
  TH1D* nJetsPF = (TH1D* ) file->Get("nJetsPF");
  TH1D* nJetsGen = (TH1D* ) file->Get("nJetsGen");


  makeplotTwo(*ZPt_Reco, *ZPt_Gen, "ZPt_spectrum", 1);
  makeplotThree(*CaloJetPt, *GenJetPt, *PFJetPt, "Jet_spectrum", 1);
  makeplotThree(*CaloJetPt2, *GenJetPt2, *PFJetPt2, "Jet_spectrum2", 1);
  makeplotThree(*CaloJetPt3, *GenJetPt3, *PFJetPt3, "Jet_spectrum3", 1);
  makeplotThree(*CaloJetEta, *GenJetEta, *PFJetEta, "Jet_eta_spectrum", 2);

  makeplotThree(*CaloJetPtRatio2over1, *GenJetPtRatio2over1, *PFJetPtRatio2over1, "JetPtRatio2over1", 0);
  makeplotThree(*CaloJetPtRatio3over2, *GenJetPtRatio3over2, *PFJetPtRatio3over2, "JetPtRatio3over2", 0);
  makeplotThree( *nJetsCalo, *nJetsGen, *nJetsPF, "JetMultiplicity", 1);

}














void makeplotThree(TH1& hist1, TH1& hist2, TH1& hist3, const char* plotname, int log) {

  hist1.SetLineColor(4);
  hist1.SetMarkerColor(4);
  hist3.SetLineColor(2);
  hist3.SetMarkerColor(2);

  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500);
  hist1.Draw( );
  hist2.Draw( "same" );
  hist3.Draw( "same" );
  TLegend *leg = new TLegend(0.48,0.7,0.89,0.92);
  leg->AddEntry( &hist1,"CaloJet","LP");
  leg->AddEntry( &hist3,"PF Jet","LP");
  leg->AddEntry( &hist2,"GenJet","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  can->SetLogy( log );
  if(log==2) {
    can->SetLogy( 1 );
    hist1.GetYaxis()->SetMoreLogLabels();
  }

  std::string plot("ratio-");
  plot.append(plotname);
  can->SaveAs( (plot+".eps").c_str() );
  can->SaveAs( (plot+".gif").c_str() );
  can->SaveAs( (plot+".root").c_str() );
  // delete can;
  cout << hist1.Integral() << endl;
}




void makeplotTwo(TH1& hist1, TH1& hist2, const char* plotname, int logy) {

  hist1.SetLineColor(4);
  hist1.SetMarkerColor(4);
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500);
  hist1.Draw( );
  hist2.Draw( "same" );
  TLegend *leg = new TLegend(0.55,0.8,0.89,0.92);
  leg->AddEntry( &hist1,"Reconstructed","LP");
  leg->AddEntry( &hist2,"Generated","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  can->SetLogy( logy );
  if(logy==2) hist1.GetYaxis()->SetMoreLogLabels();

  std::string plot("ratio-");
  plot.append(plotname);
  can->SaveAs( (plot+".eps").c_str() );
  can->SaveAs( (plot+".gif").c_str() );
  can->SaveAs( (plot+".root").c_str() );
  // delete can;
  cout << hist1.Integral() << endl;
}

