void sigDemo()
{
  const int BINS = 100;

  readOnResDKPar();

  gStyle->SetOptStat(0);

  TGraph* gDalitz = dalitzHolderN.sigGoodD0Type().pdfType()->drawBoundary(200);
  gDalitz->SetLineWidth(2);

  gROOT->cd();

  RooDataSet* data1 = pdfOnResDK.generate();

  TH1F* h =  m12->createHistogram("h",YVar(*m13));
  h->SetTitle("");

  TCanvas* can1 = new TCanvas("can1","can1",500,500);
  h->Draw();
  data1->tree().Draw("m13:m12","","same");
  TGraph *graph1 = (TGraph*)gPad->GetPrimitive("Graph");
  graph1->SetMarkerStyle(21);
  graph1->SetMarkerSize(0.2);
  graph1->SetMarkerColor(kBlue);
  gDalitz->Draw("c same");

  setVar(pdfOnResDK.parameters(), "pdfOnResDK.sigBadD0Frac",0);
  setVar(pdfOnResDK.parameters(), "pdfOnResDK.DpiGoodD0NumEvts",0);
  setVar(pdfOnResDK.parameters(), "pdfOnResDK.qqBadD0NumEvts",0);     
  setVar(pdfOnResDK.parameters(), "pdfOnResDK.totBBNumEvts",0);

  RooDataSet* data2 = pdfOnResDK.generate();

  TCanvas* can2 = new TCanvas("can2","can2",500,500);
  h->Draw();
  data2->tree().Draw("m13:m12","","same");
  TGraph *graph2 = (TGraph*)gPad->GetPrimitive("Graph");
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(0.2);
  graph2->SetMarkerColor(kBlue);
  gDalitz->Draw("c same");

  TCanvas* can3 = new TCanvas("can3","can3",500,500);
  RooDataSet* data3 = read(dataTree);
  h->Draw();
  data3->tree().Draw("m13:m12","","same");
  TGraph *graph3 = (TGraph*)gPad->GetPrimitive("Graph");
  graph3->SetMarkerStyle(21);
  graph3->SetMarkerSize(0.2);
  graph3->SetMarkerColor(kBlue);
  gDalitz->Draw("c same");


  TCanvas* can4 = new TCanvas("can4","can4",500,500);
  RooDataSet* data4 = dalitzHolderN.sigGoodD0Type().generate(10000);
  data4.append(*dalitzHolderP.sigGoodD0Type().generate(10000));
  h->Draw();
  data4->tree().Draw("m13:m12","","same");
  TGraph *graph4 = (TGraph*)gPad->GetPrimitive("Graph");
  graph4->SetMarkerStyle(21);
  graph4->SetMarkerSize(0.2);
  graph4->SetMarkerColor(kBlue);
  gDalitz->Draw("c same");


  can1->SaveAs("sigDemo-all.png");
  can2->SaveAs("sigDemo-sig.png");
  can3->SaveAs("sigDemo-data.png");
  can4->SaveAs("sigDemo-sigHigh.png");
}
