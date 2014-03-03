// $Id: testDalitzHist.cc,v 1.1 2006/04/04 05:39:32 fwinkl Exp $
// Test script for BdkDalitzHist class

void testDalitzHist() {

  const Int_t N = 10000;

  RooArgSet vars(*m12,*m13);
  m12->setBins(10);
  m13->setBins(10);

  // Create a flat Dalitz histogram
  TH2D* h = new TH2D("h","",m12->getBins(),0,3,m13->getBins(),0,3);
  h->Reset();

  TGraph* p = new TGraph(N);
  for (int i=0; i<N; i++) {    
    Double_t s12, s13;
    do {
      s12 = gRandom->Rndm()*3;
      s13 = gRandom->Rndm()*3;
    } while (!((BdkDalitzBase*)eff.getPdf())->inDalitz(s12,s13));
    h->Fill(s12,s13);
    p->SetPoint(i,s12,s13);
  }

  // Weigth bins according to their area inside Dalitz plot
  ((BdkDalitzBase*)eff.getPdf())->weightBins(h);

  // Create a histogram PDF and generate data
  BdkDalitzHist pdf("pdf","",BdkDalitz::D0,BdkDalitz::PPP0,*m12,*m13,*h);
  /*
  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  cfg->method1D().setLabel("RooSegmentedIntegrator1D");
  pdf.setIntegratorConfig(*cfg);
  pdf.forceNumInt();
  */
  RooDataSet *d = pdf.generate(vars,N);

  TCanvas* c = new TCanvas("c","c",800,800);
  c->Divide(2,2);
  c->cd(1);
  h->Draw("colz");
  drawBinning(h);
  TGraph* gr = ((BdkDalitzBase*)eff.getPdf())->drawBoundary();
  gr->Draw("c same");
  p->Draw("p same");

  c->cd(2);
  d->createHistogram(*m12,*m13)->Draw("colz");
  d->tree().Draw("m12:m13","","same");
  gr->Draw("c same");
  drawBinning(h);

  c->cd(3);
  RooPlot *p12 = m12->frame();
  d->plotOn(p12);
  pdf.plotOn(p12);
  p12->Draw();

  c->cd(4);
  RooPlot *p13 = m13->frame();
  d->plotOn(p13);
  pdf.plotOn(p13);
  p13->Draw();
}
