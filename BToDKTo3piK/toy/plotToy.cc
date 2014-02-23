// $Id: plotToy.cc,v 1.8 2006/06/23 20:45:50 fwinkl Exp $
// Plot toy results created by BdkBatchMCStudy
//
// If max<min for error or var then ranges are determined automatically


// These Gaussian are used in the fits and can be used
// to access the fit result
BdkPdfGauss* toyGaussPull;
BdkPdfGauss* toyGaussErr;
BdkPdfGauss* toyGaussVar;

TCanvas* plotToy(RooDataSet* data, const RooRealVar& var,
                 Bool_t fitPull = kTRUE, Bool_t fitErr = kTRUE,
                 Double_t errMin = 0, Double_t errMax = 0.3, Int_t errBins = 30,
                 Double_t pullMin = -4, Double_t pullMax = 4, Int_t pullBins = 40,
                 TCanvas* canvas = 0, Int_t pad = 1)
{
  return plotToy(data, var.GetName(), var.GetTitle(),
                 fitPull, fitErr,
                 errMin, errMax, errBins,
                 pullMin, pullMax, pullBins,
                 canvas, pad);
}


// Plot the pull and error distributions in 'data' for 'var'
TCanvas* plotToy(RooDataSet* data, const char* var, const char* title,
                 Bool_t fitPull = kTRUE, Bool_t fitErr = kTRUE,
                 Double_t errMin = 0, Double_t errMax = 0.3, Int_t errBins = 30,
                 Double_t pullMin = -4, Double_t pullMax = 4, Int_t pullBins = 40,
                 TCanvas* canvas = 0, Int_t pad = 1)

{
  if (!canvas) {
    canvas = new TCanvas(var,title,800,400);
    canvas->Divide(2,1);
  }
  gStyle->SetTitleYOffset(1.3);

  canvas->cd(pad);
  RooPlot* pullPlot = plotToyPull(data,var,title,fitPull,pullMin,pullMax,pullBins);  
  pullPlot->Draw();

  canvas->cd(++pad);  
  RooPlot* errPlot = plotToyErr(data,var,title,fitErr,errMin,errMax,errBins);
  errPlot->Draw();

  canvas->Update();
  return canvas;
}



TCanvas* plotToy2(RooDataSet* data, const RooRealVar& var,
                  Bool_t fitPull = kTRUE, Bool_t fitErr = kTRUE, Bool_t fitVar = kTRUE,
                  Double_t errMin = 0, Double_t errMax = 0.3, Int_t errBins = 30,
                  Double_t pullMin = -4, Double_t pullMax = 4, Int_t pullBins = 40,
                  Double_t varMin = 0, Double_t varMax = 100, Int_t varBins = 40,
                  TCanvas* canvas = 0, Int_t pad = 1)
{
  return plotToy2(data, var.GetName(), var.GetTitle(),
                  fitPull, fitErr, fitVar,
                  errMin, errMax, errBins,
                  pullMin, pullMax, pullBins,
                  varMin, varMax, varBins,
                  canvas, pad);
}


// Plot the pull, error and variable distributions in 'data' for 'var'
TCanvas* plotToy2(RooDataSet* data, const char* var, const char* title,
                  Bool_t fitPull = kTRUE, Bool_t fitErr = kTRUE, Bool_t fitVar = kTRUE,
                  Double_t errMin = 0, Double_t errMax = 0.3, Int_t errBins = 30,
                  Double_t pullMin = -4, Double_t pullMax = 4, Int_t pullBins = 40,
                  Double_t varMin = 0, Double_t varMax = 100, Int_t varBins = 40,
                  TCanvas* canvas = 0, Int_t pad = 1)

{
  if (!canvas) {
    gStyle->SetTitleYOffset(1.4);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.04);
    //    TGaxis::SetMaxDigits(3);
    canvas = new TCanvas(var,title,1200,300);
    canvas->Divide(3,1,0.005,0.005);
  }

  canvas->cd(pad);
  RooPlot* pullPlot = plotToyPull(data,var,title,fitPull,pullMin,pullMax,pullBins);  
  if (pullPlot) pullPlot->Draw();

  canvas->cd(++pad);  
  RooPlot* errPlot = plotToyErr(data,var,title,fitErr,errMin,errMax,errBins);
  if (errPlot) errPlot->Draw();

  canvas->cd(++pad);
  RooPlot* varPlot = plotToyVar(data,var,title,fitVar,varMin,varMax,varBins);
  if (varPlot) varPlot->Draw();

  canvas->Update();
  return canvas;
}


// Make the pull plot
RooPlot* plotToyPull(RooDataSet* data, const char* var, const char* title,
                     Bool_t fitPull = kTRUE, 
                     Double_t pullMin = -4, Double_t pullMax = 4, Int_t pullBins = 40)
{
  RooArgSet *args = data->get();
  RooRealVar *rPull = args->find(TString(var)+"pull");
  rPull->setRange(pullMin,pullMax);
  rPull->setBins(pullBins);

  RooPlot* pullPlot = rPull->frame();
  data->plotOn(pullPlot,MarkerSize(0.6),XErrorSize(0));
  if (fitPull) toyGaussPull = plotGauss(data,fitGauss(data, rPull),pullPlot);
  else toyGaussPull = 0;

  pullPlot->SetTitle("");
  pullPlot->GetXaxis()->SetTitle(title+TString(" pull"));
  pullPlot->GetYaxis()->SetTitle("");
  return pullPlot;
}

// Make the error plot
RooPlot* plotToyErr(RooDataSet* data, const char* var, const char* title,
                    Bool_t fitErr = kTRUE,
                    Double_t errMin = 0, Double_t errMax = 0.3, Int_t errBins = 30)
{
  RooArgSet *args = data->get();
  RooRealVar* r = (RooRealVar*)args->find(TString(var));
  if (!r->hasError()) return 0;    // no asymmetric errors at this point

  RooRealVar *rErr = args->find(TString(var)+"err");

  // Set reasonable range
  if (errMax<errMin) rErr->setRange(getMean(*data,*rErr) - 5*getRMS(*data,*rErr),
                                    getMean(*data,*rErr) + 5*getRMS(*data,*rErr));
  else rErr->setRange(errMin,errMax);

  if (fitErr) toyGaussErr = fitGauss(data, rErr);

  // Get better range in case we fit the error distribution
  if (errMax<errMin && fitErr) {
    errMin = toyGaussErr->b()->getVal() - 7*fabs(toyGaussErr->s()->getVal());
    errMax = toyGaussErr->b()->getVal() + 7*fabs(toyGaussErr->s()->getVal());
  }
    
  rErr->setRange(errMin,errMax);
  rErr->setBins(errBins);

  RooPlot* errPlot = rErr->frame(); 
  data->plotOn(errPlot,MarkerSize(0.6),XErrorSize(0));
  if (fitErr) plotGauss(data, toyGaussErr, errPlot);
  else toyGaussErr = 0;

  errPlot->SetTitle("");
  errPlot->GetXaxis()->SetTitle(title+TString(" error"));
  errPlot->GetYaxis()->SetTitle("");
  return errPlot;
}

// Plot of the fit variable itself
RooPlot* plotToyVar(RooDataSet* data, const char* var, const char* title,
                    Bool_t fitVar = kTRUE,
                    Double_t varMin = 0, Double_t varMax = 0.3, Int_t varBins = 30)
{
  RooArgSet *args = data->get();
  RooRealVar *rVar = args->find(var);

  if (varMax<varMin) rVar->setRange(getMean(*data,*rVar) - 5*getRMS(*data,*rVar),
                                    getMean(*data,*rVar) + 5*getRMS(*data,*rVar));
  else rVar->setRange(varMin, varMax);

  if (fitVar) toyGaussVar = fitGauss(data, rVar);

  if (varMax<varMin && fitVar) {
    varMin = toyGaussVar->b()->getVal() - 7*fabs(toyGaussVar->s()->getVal());
    varMax = toyGaussVar->b()->getVal() + 7*fabs(toyGaussVar->s()->getVal());
  }

  rVar->setRange(varMin,varMax);
  rVar->setBins(varBins);

  RooPlot* varPlot = rVar->frame(); 
  data->plotOn(varPlot,MarkerSize(0.6),XErrorSize(0));
  if (fitVar) plotGauss(data, toyGaussVar, varPlot);
  else toyGaussVar = 0;
  
  varPlot->SetTitle("");
  varPlot->GetXaxis()->SetTitle(title);
  varPlot->GetYaxis()->SetTitle("");
  return varPlot;
}


// Fit a Gaussian
BdkPdfGauss* fitGauss(RooDataSet* d, RooRealVar* r)
{
  BdkPdfGauss* gauss = new BdkPdfGauss(TString("gauss.")+r->GetName(),
                                       TString("Gauss ")+r->GetTitle(),*r);
  // Some reasonable start values
  gauss->b()->setVal(0.5*(r->getMax()+r->getMin()));
  gauss->s()->setVal(0.1*(fabs(r->getMax())+fabs(r->getMin())));
  
  gauss->b()->SetName("#mu");
  gauss->s()->SetName("#sigma");
  cout <<"-----------------------------------------------------------------------"<<endl;
  cout <<"     Fitting "<<r->GetName()<<endl;
  cout <<"-----------------------------------------------------------------------"<<endl;
  gauss->getPdf()->fitTo(*d,"m");

  return gauss;
}

BdkPdfGauss* plotGauss(RooDataSet* d, BdkPdfGauss* gauss, RooPlot* p)
{
  gauss->getPdf()->plotOn(p);
  gauss->getPdf()->paramOn(p,d,"",2,"NELU",0.62,0.95,0.85);
  ((TPaveText*)p->findObject("TPave"))->SetBorderSize(0);
  ((TPaveText*)p->findObject("TPave"))->SetFillStyle(kNone);
  return gauss;
}
