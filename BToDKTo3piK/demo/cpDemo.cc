void cpDemo() 
{
  const int nEvents = 20000;
  const int BINS = 200;

  BdkPdfDKDalitz& pdf = dalitzHolderN.sigGoodD0Type();

  dalitzCfg->setM23VetoMass(0,0);

  BdkDDalitzAmp * amp = pdf.pdfType()->dalitzAmp();

  TGraph* gDalitz = pdf.pdfType()->drawBoundary(200);
  gDalitz->SetLineWidth(2);

  const int nComps = amp->nComps();
  
  // Set only the amp of this resonance to 1, all others to 0:
  for (int j = 0; j < nComps; ++j){
    amp->ampRes(j)->setVal(0);
    amp->phaseRes(j)->setVal(0);
    pdf.fixAll();
  }
  
  amp->ampRes(0)->setVal(1);    // rho+
  //  amp->ampRes(3)->setVal(1);


  pdfOnResDK.xMinus()->setVal(0);
  pdfOnResDK.yMinus()->setVal(0);
  RooDataSet* data1 = pdf.generate(nEvents);

  TH2* h1 = data1->createHistogram(*m12,*m13,BINS,BINS);
  h1->SetTitle("No CP violation");
  h1->GetXaxis()->SetTitle(m12->GetTitle());
  h1->GetYaxis()->SetTitle(m13->GetTitle());


  pdfOnResDK.xMinus()->setVal(0.2*cos(0*TMath::DegToRad()));
  pdfOnResDK.yMinus()->setVal(0.2*sin(0*TMath::DegToRad()));
  RooDataSet* data2 = pdf.generate(nEvents);
  
  TH2* h2 = data2->createHistogram(*m12,*m13,BINS,BINS);
  h2->SetTitle("r_{B}=0.2, #gamma=0^{o}, #delta=0^{o}");
  h2->GetXaxis()->SetTitle(m12->GetTitle());
  h2->GetYaxis()->SetTitle(m13->GetTitle());


  pdfOnResDK.xMinus()->setVal(0.2*cos(180*TMath::DegToRad()));
  pdfOnResDK.yMinus()->setVal(0.2*sin(180*TMath::DegToRad()));
  RooDataSet* data3 = pdf.generate(nEvents);

  TH2* h3 = data3->createHistogram(*m12,*m13,BINS,BINS);
  h3->SetTitle("r_{B}=0.2, #gamma=180^{o}, #delta=0^{o}");
  h3->GetXaxis()->SetTitle(m12->GetTitle());
  h3->GetYaxis()->SetTitle(m13->GetTitle());
  
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.3,"Y");

  TCanvas* can = new TCanvas("can","can",1200,400);
  can->Divide(3,1);
  
  can->cd(1);
  h1->Draw("col");
  gDalitz->Draw("c same");
  
  can->cd(2);
  h2->Draw("col");
  gDalitz->Draw("c same");

  can->cd(3);
  h3->Draw("col");
  gDalitz->Draw("c same");

  can->SaveAs("cpDemo.eps");
  can->SaveAs("cpDemo.png");
}

