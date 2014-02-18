{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(000); 
  gStyle->SetPalette(1);
 
  double rMin = 2.2;
  double rMax = 35.0;
  int rBin = (int)(rMax - rMin);


 TF1* func = new TF1("func","x*([0]+[1]/(pow(log10(x),[2])+[3]))", 
			rMin, rMax); 
 func->SetParameter(0, -0.215828);
 func->SetParameter(1, 0.336558);
 func->SetParameter(2, 0.0767208);
 func->SetParameter(3, -0.839992);
 func->SetLineColor(4);
 func->SetLineStyle(1);
 func->SetTitle("");
 func->GetXaxis()->SetTitle("Uncorrected p_{T} (GeV/c)");
 func->GetYaxis()->SetTitle("Corrected p_{T} (GeV/c)");


  TCanvas PtVsCorrectedPt("PtVsCorrectedPt","PtVsCorrectedPt",900,600);
  func->Draw("L"); 
  PtVsCorrectedPt.SetGridx();
  PtVsCorrectedPt.SetGridy();
  PtVsCorrectedPt.SaveAs("Zjet_PtVsCorrectedPt.gif");
  PtVsCorrectedPt.SaveAs("Zjet_PtVsCorrectedPt.eps");
  PtVsCorrectedPt.SaveAs("Zjet_PtVsCorrectedPt.root");
}
