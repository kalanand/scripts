// makes a Dalitz plot of the complex amplitude of the signal PDF

void plotAmp(int nBins) {
  double step = 3./nBins;
  BdkDKDalitz * pdf = dalitzHolderN.sigGoodD0Type().getPdf();

  TH2F * hReal = new TH2F("hReal", "Real", nBins, 0, 3, nBins, 0, 3);
  TH2F * hImag = new TH2F("hImag", "Imag", nBins, 0, 3, nBins, 0, 3);
  TH2F * hMag = new TH2F("hMag", "Mag", nBins, 0, 3, nBins, 0, 3);
  TH2F * hPhase = new TH2F("hPhase", "Phase", nBins, 0, 3, nBins, 0, 3);

  for (int ix = 0; ix < nBins; ++ix) {
    double x = (ix + 0.5) * step;
    m12->setVal(x);
    for (int iy = 0; iy < nBins; ++iy) {
      double y = (iy + 0.5) * step;
      m13->setVal(y);
      
      pdf->getVal();
      RooComplex amp = pdf->lastAmp();
      
      hReal->SetBinContent(ix, iy, amp.re());
      hImag->SetBinContent(ix, iy, amp.im());
      hMag->SetBinContent(ix, iy, amp.abs());
      hPhase->SetBinContent(ix, iy, atan2(amp.im(), amp.re()));
      
    }
  }

  gStyle->SetPalette(1);

  TCanvas * can = new TCanvas("can", "can", 1000, 1000);
  can->Divide(2,2);
  
  can->cd(1);
  hReal->Draw("colz");
  can->cd(2);
  hImag->Draw("colz");
  can->cd(3);
  hMag->Draw("colz");
  can->cd(4);
  hPhase->Draw("colz");
}
      
      
  
  


