
Int_t contourPoints = 30;

// Plot x/y contours for pdfOnResDK
// use generateAndFit=false if you just did a fit
void fullContours(Bool_t generateAndFit = kTRUE)
{
  readOnResDKPar();

  RooRealVar xm(*dalitzHolderN.sigGoodD0Type().x());
  RooRealVar ym(*dalitzHolderN.sigGoodD0Type().y());
  RooRealVar xp(*dalitzHolderP.sigGoodD0Type().x());
  RooRealVar yp(*dalitzHolderP.sigGoodD0Type().y());

  xm.SetTitle("x-");
  ym.SetTitle("y-");  
  xp.SetTitle("x+");
  yp.SetTitle("y+");  

  xm.setRange(-1,1);
  ym.setRange(-1,1);
  xp.setRange(-1,1);
  yp.setRange(-1,1);

  RooDataSet* data = 0;
  if (generateAndFit) {
    // Generate data
    gROOT->cd();
    data = pdfOnResDK.generate();
    
    // Fit data
    fitOption = "mer";
    optOption = "c";
    fit(pdfOnResDK, *data);
  }


  // Mark final point
  TMarker *fitPointN = new TMarker(dalitzHolderN.sigGoodD0Type().x()->getVal(),
                                   dalitzHolderN.sigGoodD0Type().y()->getVal(),
                                   8);
  TMarker *fitPointP = new TMarker(dalitzHolderP.sigGoodD0Type().x()->getVal(),
                                   dalitzHolderP.sigGoodD0Type().y()->getVal(),
                                   8);

  
  
  // 1- and 2-sigma contours
  RooArgList floatParamList(pdfOnResDK.parametersFree());
  Int_t ixm = floatParamList.index(floatParamList.find(xm.GetName()));
  Int_t iym = floatParamList.index(floatParamList.find(ym.GetName()));
  Int_t ixp = floatParamList.index(floatParamList.find(xp.GetName()));
  Int_t iyp = floatParamList.index(floatParamList.find(yp.GetName()));

  /*
  pdfOnResDK.fixAll();
  dalitzHolderN.sigGoodD0Type().x()->setConstant(false);
  dalitzHolderN.sigGoodD0Type().y()->setConstant(false);
  minuit->migrad();
  */
  TGraph* gm1 = contourGraph(ixm, iym, 1, contourPoints);
  TGraph* gm2 = contourGraph(ixm, iym, 2, contourPoints);
  /*
  pdfOnResDK.fixAll();
  dalitzHolderP.sigGoodD0Type().x()->setConstant(false);
  dalitzHolderP.sigGoodD0Type().y()->setConstant(false);
  minuit->migrad();
  */
  TGraph* gp1 = contourGraph(ixp, iyp, 1, contourPoints);
  TGraph* gp2 = contourGraph(ixp, iyp, 2, contourPoints);

  
  // NLL plots
  //  RooPlot* px = pdfOnResDK.getPdf()->plotNLLOn(xm.frame(), data);
  //  RooPlot* py = pdfOnResDK.getPdf()->plotNLLOn(ym.frame(), data);

  TCanvas* can = new TCanvas("can","Signal toy MC",1200,400);
  can->Divide(3,1,0.002,0.002);

  can->cd(1);
  TH2F *frame = xm.createHistogram("contourPlot", ym);

  frame->SetStats(kFALSE);
  frame->SetTitle("a) Confidence regions");
  frame->Draw();

  gm1->SetFillColor(kBlue);
  gm2->SetFillColor(38);  
  gm2->Draw("lf");
  gm1->Draw("lf");

  gp1->SetFillColor(kBlue);
  gp2->SetFillColor(38);  
  gp2->Draw("f");
  gp1->Draw("f");


  fitPointN->Draw();
  fitPointP->Draw();
  /*
  can->cd(2);
  px->SetTitle("b) NLL projection for x-");
  px->Draw();

  can->cd(3);
  py->SetTitle("c) NLL projection for y-");
  py->Draw();

  */
  can->SaveAs("fullContours.eps");
  can->SaveAs("fullContours.root");
  delete data;
}


// Plot confidence regions and NLL for B- signal only
void sigContours(Int_t nEvents = 100)
{
  const BdkPdfDKDalitz& pdf = dalitzHolderN.sigGoodD0Type();

  readOnResDKPar();

  RooRealVar xm(*pdf.x());  
  RooRealVar ym(*pdf.y());
  xm.SetTitle("x-");
  ym.SetTitle("y-");  
  xm.setRange(-1,1);
  ym.setRange(-1,1);

  // Generate data
  gROOT->cd();
  RooDataSet* data = pdf.generate(nEvents);

  // Mark initial point
  TMarker *genPoint = new TMarker(pdf.x()->getVal(), 
                                  pdf.y()->getVal(), 
                                  24);

  // Fit data
  fitOption = "r";
  optOption = "c";
  fit(pdf, *data);
  
  // Mark final point
  TMarker *fitPoint = new TMarker(pdf.x()->getVal(), 
                                  pdf.y()->getVal(), 
                                  8);

  // 1- and 2-sigma contours
  TGraph* g1 = contourGraph(pdf, xm, ym, 1, contourPoints);
  TGraph* g2 = contourGraph(pdf, xm, ym, 2, contourPoints);
  
  // NLL plots
  RooPlot* px = pdf.getPdf()->plotNLLOn(xm.frame(), data);
  RooPlot* py = pdf.getPdf()->plotNLLOn(ym.frame(), data);

  TCanvas* can = new TCanvas("can","Signal toy MC",1200,400);
  can->Divide(3,1,0.002,0.002);

  can->cd(1);
  xm.setRange(-0.5,0.5);
  ym.setRange(-0.5,0.5);
  TH2F *frame = xm.createHistogram("contourPlot", ym);

  frame->SetStats(kFALSE);
  frame->SetTitle("a) Confidence regions");
  frame->Draw();

  g1->SetFillColor(kBlue);
  g2->SetFillColor(38);
 
  g2->Draw("f");
  g1->Draw("f");
  fitPoint->Draw();
  genPoint->Draw();

  can->cd(2);
  px->SetTitle("b) NLL projection for x-");
  px->Draw();

  can->cd(3);
  py->SetTitle("c) NLL projection for y-");
  py->Draw();


  can->SaveAs("sigContours.eps");
  can->SaveAs("sigContours.root");
  delete data;
}
