// $Id: testDalitzSmear.cc,v 1.1 2006/03/20 19:05:12 fwinkl Exp $
// Test script for Bdk(Pdf)DalitzSmear class

void testDalitzSmear() {

  BdkPdfDKDalitz& pdf = dalitzHolderP.sigGoodD0Type();

  // Resolution models
  BdkPdfGauss res12("res12","",*m12);
  BdkPdfGauss res13("res13","",*m13);
  res12.b()->setVal(0);
  res12.s()->setVal(0.1);
  res13.b()->setVal(0);
  res13.s()->setVal(0.1);

  
  BdkPdfDalitzSmear smear("smear","",pdf,res12,res13);
  smear.setEventBuffer(6000);

  // Generate original and smeared datasets
  RooDataSet *d1 = pdf.generate(5000);
  RooDataSet *d2 = smear.generate(5000);

  // Plot Dalitz plot and projections of original and smeared datasets
  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide(2,2);

  can->cd(1);
  TH2* h1 = d1->createHistogram(*m12,*m13);
  h1->Draw();
  ((BdkDalitzBase*)pdf.getPdf())->drawBoundary()->Draw("c same");
  
  can->cd(2);
  TH2* h2 = d2->createHistogram(*m12,*m13);
  h2->SetMarkerColor(kBlue);
  h2->Draw();
  ((BdkDalitzBase*)pdf.getPdf())->drawBoundary()->Draw("c same");
  
  can->cd(3);
  RooPlot *p12 = m12->frame();
  d1->plotOn(p12);
  d2->plotOn(p12,MarkerColor(kBlue));
  p12->Draw();

  can->cd(4);
  RooPlot *p13 = m13->frame();
  d1->plotOn(p13);
  d2->plotOn(p13,MarkerColor(kBlue));
  p13->Draw();

}
