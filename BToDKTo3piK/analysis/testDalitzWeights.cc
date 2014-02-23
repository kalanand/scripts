// $Id: testDalitzWeights.cc,v 1.2 2006/02/19 01:01:49 fwinkl Exp $
// Standalone test script for Dalitz plot bin reweighting


void testDalitzWeights() {

  RooRealVar* s12 = new RooRealVar("s12","",0,3);
  RooRealVar* s13 = new RooRealVar("s13","",0,3);

  // Define some fancy binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);
  bins12.addUniform(20,0,1);  
  bins13.addUniform(20,0,1);
  bins12.addUniform(10,1,2);
  bins13.addUniform(10,1,2);

  // Set the binning
  s12->setBinning(bins12);
  s13->setBinning(bins13);

  // Empty dataset
  RooDataSet* data = new RooDataSet("data","",RooArgSet(*s12,*s13));

  // Plotting
  TCanvas *can = new TCanvas("can","can",500,500);
  gStyle->SetOptStat(0);
  // Create histogram with binning
  TH2* h = data->createHistogram("h",*s12,Binning(bins12),
				  YVar(*s13,Binning(bins13)));

  // Dummy PDF so that we have access to the BdkDalitzBase functionality
  BdkPdf2DpolyDalitz pdf("pdf","",*s12,*s13,BdkDalitzBase::D0); 
    
  ((BdkDalitzBase*)pdf.getPdf())->weightBins(h,false);
  h->Draw("colz");
  drawBinning(h);
  
  TGraph *g = ((BdkDalitzBase*)pdf.getPdf())->drawBoundary();
  g->SetLineWidth(2);
  g->Draw("c same");
}


void drawBinning(TH2 *h) {

  if (!h) return;

  TLine line;
  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();

  for (Int_t i=ax->GetFirst()+1; i<=ax->GetNbins(); i++) {
    line.DrawLine(ax->GetBinLowEdge(i),ay->GetBinLowEdge(ay->GetFirst()),
                  ax->GetBinLowEdge(i),ay->GetBinUpEdge(ay->GetLast()));
  }

  for (Int_t i=ay->GetFirst()+1; i<=ay->GetNbins(); i++) {
    line.DrawLine(ax->GetBinLowEdge(ax->GetFirst()),ay->GetBinLowEdge(i),
                  ax->GetBinUpEdge(ax->GetLast()),ay->GetBinLowEdge(i));
  }
}
