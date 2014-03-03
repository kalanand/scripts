// $Id: plotContour.cc,v 1.1 2006/05/03 17:37:43 fwinkl Exp $
// plot n1- and n2-sigma contours of pdf for var1 and var2


// use data to fit
void plotContour(BdkPdfAbsBase& pdf,
		 RooRealVar& var1, RooRealVar& var2,
		 RooAbsData& data, Double_t n1 = 1, Double_t n2 = 2,
		 Bool_t plotLegend = kTRUE, Bool_t drawInitialPoint = kFALSE)
{
  // generated point
  TMarker *genPoint = new TMarker(var1.getVal(), var2.getVal(), 24);

  // fit
  fit(pdf, data);

  // plot contour
  plotContour(pdf, var1, var2, n1, n2, plotLegend);

  if (drawInitialPoint) genPoint->Draw();
  else delete genPoint;
}
 

// use a previously done fit
// adapted from RooMinuit::contour()
void plotContour(BdkPdfAbsBase& pdf,
		 RooRealVar& var1, RooRealVar& var2,
		 Double_t n1 = 1, Double_t n2 = 2,
		 Bool_t plotLegend = kTRUE)
{
  const Int_t contourPoints = 30;

  // create and draw a frame
  TH2F *frame = var1.createHistogram("contourPlot", var2, "-log(likelihood)") ;
  frame->SetStats(kFALSE);
  frame->SetTitle("");
  
  TGraph* graph1 = 0;
  if(n1 > 0) {
    graph1 = contourGraph(pdf, var1, var2, n1, contourPoints);
    graph1->SetName("contourN1");
    TString s;
    s.Form("%1.1f #sigma",n1);
    graph1->SetTitle(s);
  }
  
  TGraph* graph2 = 0;
  if(n2 > 0) {
    graph2 = contourGraph(pdf, var1, var2, n2, contourPoints);
    graph2->SetName("contourN2");
    TString s;
    s.Form("%1.1f #sigma",n2);
    graph2->SetTitle(s);
  }
  

  // Draw all objects
  frame->Draw();

  // Legend
  TLegend* leg = new TLegend(0.71,0.85,0.86,0.6);
  leg->SetFillStyle(kNone);

  if (graph2) {
    leg->AddEntry(graph2,graph2->GetTitle(),"f");
    graph2->SetFillColor(38);
    graph2->Draw("f") ;
  }
  if (graph1) {
    leg->AddEntry(graph1,graph1->GetTitle(),"f");
    graph1->SetFillColor(kBlue);
    graph1->Draw("f") ;
  }

  if (plotLegend) leg->Draw();
  else delete leg;  
}



// Draw nSigma contour graph for var1 and var2
TGraph* contourGraph(BdkPdfAbsBase& pdf,
                     RooRealVar& var1, RooRealVar& var2,
                     Double_t nSigma, Int_t contourPoints = 30)
{
  // Verify that both variables are floating parameters of PDF
  RooArgList floatParamList(pdf.parametersFree());
  Int_t index1 = floatParamList.index(floatParamList.find(var1.GetName()));
  if(index1 < 0) {
    cout << "plotContour ERROR: " 
	 << var1.GetName() << " is not a floating parameter of PDF " 
         << pdf.GetName() << endl ;
    return 0;
  }

  Int_t index2 = floatParamList.index(floatParamList.find(var2.GetName()));
  if(index2 < 0) {
    cout << "plotContour ERROR: " 
	 << var2.GetName() << " is not a floating parameter of PDF " 
         << pdf.GetName() << endl ;
    return 0;
  }

  return contourGraph(index1, index2, nSigma, contourPoints);
}



// Draw nSigma contour graph for variables with index index1 and index2
TGraph* contourGraph(Int_t index1, Int_t index2, 
                     Double_t nSigma, Int_t contourPoints = 30)
{
  if (!gMinuit) return 0;
  if (nSigma<=0) return 0;
  if (contourPoints<=0) return 0;

  // remember our original value of ERRDEF
  Double_t errdef = gMinuit->fUp;

  gMinuit->SetErrorDef(nSigma*nSigma*errdef);
  TGraph* g = (TGraph*)gMinuit->Contour(contourPoints, index1, index2);

  // restore the original ERRDEF
  gMinuit->SetErrorDef(errdef);

  return g;
}
