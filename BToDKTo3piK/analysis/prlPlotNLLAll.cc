// $Id: prlPlotNLLAll.cc,v 1.1 2006/09/13 22:39:34 fwinkl Exp $
// All six NLL plots in PRL quality

// ugly global variable use by plot()
int nCont = -1;


void prlPlotNLLAll(int xnBins = 20, double xmin = 0.5, double xmax = 1.5,
                   int ynBins = 20, double ymin = 90, double ymax = 270,
                   int nContours = -1,
                   const char * filePrefix = 0) {

  const Bool_t nll1Donly = kFALSE;
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
  gStyle->SetTitleOffset(0.8,"Y");
  gStyle->SetTitleOffset(0.8,"X");
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
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
  
  plotAll(prefix, sumNLL, *xN, *yN, *xP, *yP, orig, drawCont, nll1Donly);
}



void plotAll(const TString & prefix,
	     RooAbsReal & nll,
	     RooRealVar& xN, RooRealVar& yN, 
	     RooRealVar& xP, RooRealVar& yP,
	     const RooAbsCollection * orig,
	     Bool_t drawCont,
             Bool_t nll1Donly) {

  if (nll1Donly==false) {
    TH2F * LXnYn = plot(prefix, nll, xN, yN, orig, drawCont);

    TH2F * LXnXp = plot(prefix, nll, xN, xP, orig, drawCont);

    TH2F * LXnYp = plot(prefix, nll, xN, yP, orig, drawCont);
    
    TH2F * LYnXp = plot(prefix, nll, yN, xP, orig, drawCont);
    
    TH2F * LYnYp = plot(prefix, nll, yN, yP, orig, drawCont);
    
    TH2F * LXpYp = plot(prefix, nll, xP, yP, orig, drawCont);

    TCanvas * can2 = new TCanvas(prefix + "NLL2D", prefix + "NLL2D", 
                                 1200, 800);
    can2->Divide(3,2,0.002,0.002);
    
    can2->cd(1);
    plotConts(LXnYn, drawCont);
    
    can2->cd(2);
    plotConts(LXnXp, drawCont);
    
    can2->cd(3);
    plotConts(LXnYp, drawCont);
    
    can2->cd(4);
    plotConts(LYnXp, drawCont);
    
    can2->cd(5);
    plotConts(LYnYp, drawCont);
    
    can2->cd(6);
    plotConts(LXpYp, drawCont);
    
    can2->Update();
    can2->SaveAs(prefix + "NLL2D-vs-XY.eps");
    can2->SaveAs(prefix + "NLL2D-vs-XY.root");
  }


  TH1F * LXn = plot(prefix, nll, xN, orig); 			    
  
  TH1F * LYn = plot(prefix, nll, yN, orig); 			    
  
  TH1F * LXp = plot(prefix, nll, xP, orig); 			    
  
  TH1F * LYp = plot(prefix, nll, yP, orig); 			    
  
  TCanvas * can1 = new TCanvas(prefix + "NLL1D", prefix + "NLL1D", 
			       800, 800);
  can1->Divide(2,2,0.002,0.002);
  
  can1->cd(1);
  if (LXn) LXn->Draw();
  
  can1->cd(2);
  if (LYn) LYn->Draw();
  
  can1->cd(3);
  if (LXp) LXp->Draw();
  
  can1->cd(4);
  if (LYp) LYp->Draw();
  
  can1->Update();
  can1->SaveAs(prefix + "NLL1D-vs-XY.eps");
  can1->SaveAs(prefix + "NLL1D-vs-XY.root");
}


void plotConts(TH2F* hist, Bool_t drawCont) {

  if (hist==0) return;

  const char * option = "colz";
  if (drawCont) {
    option = "cont4";
  }
    
  TH2F* h = hist;
  TH2F* conts = 0;
  if (drawCont) {
    h = (TH2F*)hist->Clone(hist->GetName()+TString("_clone"));
    conts = (TH2F*)hist->Clone(hist->GetName()+TString("_conts"));
  }

  if (drawCont) {
    // Remove Axis from original histogram
    h->GetXaxis()->SetAxisColor(0);
    h->GetXaxis()->SetLabelOffset(999);
    h->GetXaxis()->SetTickLength(0);
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetAxisColor(0);
    h->GetYaxis()->SetLabelOffset(999);
    h->GetYaxis()->SetTickLength(0); 
    h->GetYaxis()->SetTitle("");
  }
  h->Draw(option);

  if (drawCont) {
    // Put contour lines in a transparent pad on top
    TPad* pad = new TPad(hist->GetName()+TString("_pad"),"",0,0,1,1);
    pad->SetFillStyle(0);
    pad->SetFrameFillStyle(0);
    pad->Draw();
    pad->cd(); 
    conts->Draw("cont3");
  }
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

    TArrayD levels(nLevels+1);
    for (int l=0; l<=nLevels; l++) levels[l] = l*l/2.0;
    result->SetContour(nLevels+1, levels.GetArray());
  }

  result->SetStats(kFALSE);

  resetVals(orig);
  return result;
}


// make 1D plot:
TH1F * plot(TString prefix,
	    RooAbsReal & nll, RooRealVar & x,
	    const RooAbsCollection * orig) {

  TString name = prefix + TString("NLL, ") + x.GetTitle();
  TString title = name;

  double max = x.getMax("plotNLL");
  double min = x.getMin("plotNLL");
  int nBins = x.getBins();
  double step = (max - min) / nBins;

  TH1F * result = new TH1F(name, title, nBins, min, max);
  result->GetXaxis()->SetTitle(x.GetTitle());

  double baseline = DBL_MAX;


  for (int ix = 1; ix <= nBins; ++ix){
    double xx = min + (ix - 0.5) * step;
    x.setVal(xx);
    double val = nll.getVal();
    if (baseline > val) {
      baseline = val;
    }
    
    result->SetBinContent(ix, val);
  }

  for (int ix = 1; ix <= nBins; ++ix){
    double newVal = result->GetBinContent(ix) - baseline;
    result->SetBinContent(ix, newVal);
  }    
  cout << "Subtracted " << baseline << " from " << result->GetName() << endl;

  result->SetStats(kFALSE);
  resetVals(orig);
  return result;
}


