// $Id: testBlattWeisskopf.cc,v 1.1 2006/05/03 21:41:20 fwinkl Exp $
// Test script for the Blatt-Weisskopf penetration factors

void testBlattWeisskopf(Int_t nEvents = 1000)
{
  BdkPdfDKDalitz& pdf = dalitzHolderP.sigGoodD0Type();

  setVar(pdf.parameters(),"dalitzHolderP.sigGoodD0.pdf.dalitzAmp.resRadius",0);
  pdf.dalitzAmp()->calNorm();
  setRandomGenSeed(712);
  RooDataSet* data1 = pdf.generate(nEvents);

  setVar(pdf.parameters(),"dalitzHolderP.sigGoodD0.pdf.dalitzAmp.resRadius",1.5);
  pdf.dalitzAmp()->calNorm();
  setRandomGenSeed(712);
  RooDataSet* data2 = pdf.generate(nEvents);

  TCanvas* can = new TCanvas("can","",1200,400);
  can->Divide(3,1);

  can->cd(1);
  TH2F* h1 = data1->createHistogram(*m12,*m13);
  h1->Draw();

  can->cd(2);
  TH2F* h2 = data2->createHistogram(*m12,*m13);
  h2->Draw();

  can->cd(3);
  TH2F* hchi2 = h2->Clone("hchi2");
  chi2test2d(h1,h2,hchi2);
  hchi2->Draw("colz");
    
  delete data1;
  delete data2;
}
