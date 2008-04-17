// $Id: prlPlotNLL.cc,v 1.1 2006/09/13 22:39:34 fwinkl Exp $
// rho/theta NLLs in one plot


// ugly global variable use by plot()

int nCont = -1;


void prlPlotNLL(int xnBins = 20, double xmin = 0.5, double xmax = 1.2,
                int ynBins = 20, double ymin = 90, double ymax = 270,
                int nContours = 3,
                const char * filePrefix = 0) {

  const Bool_t drawCont = kTRUE;
  nCont = nContours;
    
  RooRealVar* xN;
  RooRealVar* yN;
  RooRealVar* xP;
  RooRealVar* yP;

  if (pdfOnResDK.xMinus()!=0) {  // cartesian coordinates
    xN = pdfOnResDK.xMinus(); xN->SetTitle("x^{-}");
    yN = pdfOnResDK.yMinus(); yN->SetTitle("y^{-}");
    xP = pdfOnResDK.xPlus();  xP->SetTitle("x^{+}");
    yP = pdfOnResDK.yPlus();  yP->SetTitle("y^{+}");
  }
  else {   // polar coordinates
    xN = pdfOnResDK.rhoMinus();   xN->SetTitle("#rho^{-}");
    yN = pdfOnResDK.thetaMinus(); yN->SetTitle("#theta^{-}");
    xP = pdfOnResDK.rhoPlus();    xP->SetTitle("#rho^{+}");
    yP = pdfOnResDK.thetaPlus();  yP->SetTitle("#theta^{+}");
  }

  // Now set the ranges and binnings (use a named range)
  xN->setRange("plotNLL",xmin,xmax); xN->setBins(xnBins);
  yN->setRange("plotNLL",ymin,ymax); yN->setBins(ynBins);
  xP->setRange("plotNLL",xmin,xmax); xP->setBins(xnBins);
  yP->setRange("plotNLL",ymin,ymax); yP->setBins(ynBins);


  readCut = cutSigReg;
  data = read(dataTree);

  BdkOnResNLLYields nllYields(pdfOnResDK, -0.038);
  nllYields.Print();  
  ///  nll.setVerbose("V");

  // Get the RRV's and their initial values:
  
  
  RooAbsCollection* orig = pdfOnResDK.cpParams().snapshot(false);

  gROOT->SetStyle("BABAR");
  gStyle->SetPalette(1);
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetNdivisions(506,"XY");

  // Plot the NLL's:

  TString prefix;
  TString prePrefix = filePrefix;
  
  // Prepare the data and the NLLVar:
  RooAbsReal * nllVar = 0;

  nllVar = new RooNLLVar("nllVar","-log(likelihood)",
                         *(pdfOnResDK.getPdf()), 
                         *data,
                         RooArgSet());

  // Plot with both:
  
  prefix = prePrefix + "Yields-Shape-";
  
  RooFormulaVar sumNLL("nll+penalty", "nll+penalty",
                       TString(nllVar->GetName()) + "+" +
                       TString(nllYields.GetName()), 
                       RooArgSet(*nllVar, nllYields));
  
  plotAll(prefix, sumNLL, *xN, *yN, *xP, *yP, orig, drawCont);
}



void plotAll(const TString & prefix,
	     RooAbsReal & nll,
	     RooRealVar& xN, RooRealVar& yN, 
	     RooRealVar& xP, RooRealVar& yP,
	     const RooAbsCollection * orig,
	     Bool_t drawCont) {

  TH2F * LXnYn = plot(prefix, nll, xN, yN, orig, drawCont);
  
  TH2F * LXpYp = plot(prefix, nll, xP, yP, orig, drawCont);
  
  TCanvas * can2 = new TCanvas(prefix + "NLL2D", prefix + "NLL2D", 
                               400, 400);

  LXnYn->SetLineStyle(7);
  LXnYn->SetLineColor(kBlue);  
  LXnYn->SetLineWidth(2);
  LXnYn->GetXaxis()->SetTitle("#rho_{#pm}");
  LXnYn->GetYaxis()->SetTitle("#theta_{#pm} (degree)");

  LXpYp->SetLineColor(kRed);
  LXpYp->SetLineWidth(2);

  LXnYn->Draw("cont3");
  LXpYp->Draw("cont3 same");

  // No-CP marker
  //  TMarker* nocp = new TMarker(0.85,180,5);
  //  nocp->SetMarkerSize(2.0);
  //  nocp->Draw();

  const char* range = "plotNLL";
  Double_t xlen = 0.03*(xN.getMax(range)-xN.getMin(range));
  Double_t ylen = 0.03*(yN.getMax(range)-yN.getMin(range));

  TLine l;
  l.SetLineWidth(2);
  l.DrawLine(0.85-xlen,180-ylen,0.85+xlen,180+ylen);
  l.DrawLine(0.85-xlen,180+ylen,0.85+xlen,180-ylen);

  // Labels
  TLatex* texBp = new TLatex(1.0,210,"B^{+}");
  texBp->SetTextSize(0.07);
  texBp->SetTextColor(kRed);
  texBp->Draw();

  TLatex* texBm = new TLatex(1.1,232,"B^{-}");
  texBm->SetTextSize(0.07);
  texBm->SetTextColor(kBlue);
  texBm->Draw();


  can2->Update();
  can2->SaveAs(prefix + "NLL2D-vs-XY.eps");
  can2->SaveAs(prefix + "NLL2D-vs-XY.root");
}




void resetVals(const RooAbsCollection* orig) {
  if (orig) pdfOnResDK.cpParams() = *(RooArgSet*)orig;
}


// make 2D plot:
TH2F * plot(TString prefix,
	    RooAbsReal & nll, RooRealVar & x, RooRealVar & y,
	    const RooAbsCollection * orig, Bool_t drawCont) {

  // Use the const status of x and y to tell us if to plot them or not:
  if (x.isConstant() || y.isConstant()) {
    return 0;
  }

  const char* range = "plotNLL";
  Double_t xmin = x.getMin(range);
  Double_t xmax = x.getMax(range);
  Double_t ymin = y.getMin(range);
  Double_t ymax = y.getMax(range);
  Int_t xnBins = x.getBins();
  Int_t ynBins = y.getBins();
  
  TString name = prefix + "NLL, " + y.GetTitle() + TString(" vs. ") + x.GetTitle();
  TString title = name;

  TH2F * result = new TH2F(name, title, xnBins, xmin, xmax, ynBins, ymin, ymax);
  result->GetXaxis()->SetTitle(x.GetTitle());
  result->GetYaxis()->SetTitle(y.GetTitle());

  double baseline = DBL_MAX;

  double xstep = (xmax - xmin) / xnBins;
  double ystep = (ymax - ymin) / ynBins;

  for (int ix = 1; ix <= xnBins; ++ix){
    double xx = xmin + (ix - 0.5) * xstep;
    x.setVal(xx);
    for (int iy = 1; iy <= ynBins; ++iy){
      double yy = ymin + (iy - 0.5) * ystep;
      y.setVal(yy);

      double val = nll.getVal();
      if (baseline > val) {
	baseline = val;
      }
      
      result->SetBinContent(ix, iy, val);
    }
  }

  for (int ix = 1; ix <= xnBins; ++ix){
    for (int iy = 1; iy <= ynBins; ++iy){
      double newVal = result->GetBinContent(ix, iy) - baseline;
      result->SetBinContent(ix, iy, newVal);
    }
  }
  cout << "Subtracted " << baseline << " from " << result->GetName() << endl;

  if (drawCont) {    
    int nLevels = (int)sqrt(2*result->GetMaximum());

    // Overwrite if requested
    if (nCont>=0) nLevels = nCont;

    cout << "Maximum "<< result->GetName() <<" NLL is "<<result->GetMaximum()
         << ". Using "<<nLevels<<" contours."<<endl;

    // 0 are the default contours of root
    if (nLevels==0) nLevels++;   

    TArrayD levels(nLevels);
    for (int l=1; l<=nLevels; l++) levels[l-1] = l*l/2.0;
    result->SetContour(nLevels, levels.GetArray());
  }

  result->SetStats(kFALSE);

  resetVals(orig);
  return result;
}


