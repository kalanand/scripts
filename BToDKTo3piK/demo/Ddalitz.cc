// D Dalitz plot

void Ddalitz()
{
  data = dalitzHolderN.DpiGoodD0Type().generate(20000);

  TH2 * hist = data->createHistogram(*m12,*m13,200,200);

  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(m12->GetTitle());
  hist->GetYaxis()->SetTitle(m13->GetTitle());

  TGraph* gDalitz = dalitzHolderN.DpiGoodD0Type().pdfType()->drawBoundary(200);
  gDalitz->SetLineWidth(2);

  TCanvas* can = new TCanvas("Ddalitz", "D Dalitz", 500, 500);  
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.2,"Y");
  hist->Draw("col");
  gDalitz->Draw("c same");   // boundary
  
  can->SaveAs("Ddalitz.eps");
}
