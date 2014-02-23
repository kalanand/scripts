// $Id: plotDalitz3D.cc,v 1.1 2006/04/20 01:07:43 fwinkl Exp $

// Make a 3D Dalitz plot of data and pdf
void plotDalitz3D(RooAbsData* data, RooAbsPdf* pdf,
                  RooRealVar* x = m12, RooRealVar* y = m13,
                  Int_t binsx = -1, Int_t binsy = -1)
{
  if (binsx<=0) binsx = x->getBins();
  if (binsy<=0) binsy = y->getBins();

  TH1* hdata = data->createHistogram("hdata",*x,Binning(binsx),YVar(*y,Binning(binsy)));
  TH1* hpdf = pdf->createHistogram("hpdf",*x,Binning(binsx),YVar(*y,Binning(binsy)));
  
  hdata->SetLineColor(kBlue);
  hpdf->SetTitle("");
  hpdf->GetXaxis()->SetTitleOffset(1.6);
  hpdf->GetYaxis()->SetTitleOffset(1.6);
  hpdf->GetZaxis()->SetTitle("");
  
  hpdf->GetZaxis()->SetAxisColor(0);
  hpdf->GetZaxis()->SetLabelOffset(999);
  hpdf->GetZaxis()->SetTickLength(0);
  
  hpdf->Draw("surf4 fb bb");
  hdata->Draw("lego fb bb same");

}
