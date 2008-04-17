// $Id: plotNLL_RhoAmpsPhases2.cc,v 1.2 2006/06/28 23:03:17 abi Exp $


#include "../BToDKTo3piK/analysis/plotNLL.cc"

TH2F * histRhoAmpsPhases2 = 0;

void plotNLL_RhoAmpsPhases2(Bool_t varyPhase = kTRUE,
			    Bool_t signalOnly = kTRUE,
			    int nBins = 20, double min = -2, double max = 2,
			    double minRatio = 0, double maxRatio = 1,
			    double ratioStep = 0.1,
			    int nEvts = 1700,
			    int nEvtsNorm = 170,
			    const char * parFile = "../BToDKTo3piK/params/cross/plotNLL_RhoPhases2.par",
			    const char * option = "cont1") {

  TString bgd = "";

  // Get the PDF:
  BdkPdfAbsBase * pdf = &(dalitzHolderN.sigGoodD0Type());
  if (kFALSE == signalOnly) {
    pdf = &pdfOnResDK;
    bgd = "+Bgd-";
  }

  // Read the parameters we want to start with:
  if (parFile) {
    pdf->parameters().readFromFile(parFile);
    cout << "Read parameters from " << parFile << endl;
  }
  else {
    cout << "Using current PDF  parameters" << endl;
  }
    

  // Get the charged rho amps:
  RooRealVar* ampMinus = (RooRealVar*)pdf->parameters().
    find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_amp");
    
  RooRealVar* ampPlus = (RooRealVar*)pdf->parameters().
    find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho+_amp");
    
  // Loop over ratio of rho- to rho+ amps and plot NLL for the
  // different rho- phases:
  for (double ratio = minRatio; ratio <= maxRatio; ratio += ratioStep) {
    ampMinus->setVal(ampPlus->getVal() * ratio);

    TString theName = bgd;
    theName += "plotNLL_Rho-_amp=";
    char ratioString[100];
    sprintf(ratioString, "%0.3f", ratio);
    theName += ratioString;

    plotNLL_RhoPhases2(*pdf, varyPhase, signalOnly, nBins, min, max, 
		       nEvts, nEvtsNorm, option, theName); 
  }
}


// 2D-Plot the NLL, changing rho- phase. Otherwise takes whatever other
// parameter values are in place:
void plotNLL_RhoPhases2(BdkPdfAbsBase & pdf, 
			Bool_t varyPhase = kTRUE,
			Bool_t signalOnly = kTRUE,
			int nBins = 20, double min = -2, double max = 2,
			int nEvts = 1700,
			int nEvtsNorm = 170,
			const char * option = "cont1", 
			const TString & theName = "") {
  
  gStyle->SetPalette(1);
  TCanvas * c;

  if (varyPhase) {
    c = new TCanvas(theName, theName, 1450, 580);
    c->Divide(5,2);
  }
  else {
    c = new TCanvas(theName, theName, 700, 700);    
  } 

  // The rho- phase:
  RooRealVar* phase = (RooRealVar*)pdf.parameters().
    find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_phase");

  // Get x and y:
  RooRealVar * xN = pdfOnResDK.xMinus();
  RooRealVar * yN = pdfOnResDK.yMinus();
  RooRealVar * xP = pdfOnResDK.xPlus();
  RooRealVar * yP = pdfOnResDK.yPlus();

  const double orig[4] = {xN->getVal(), 
			  yN->getVal(), 
			  xP->getVal(), 
			  yP->getVal()};

  // Loop over rho- phase values:
  int i = 0;
  double minPhase = 0;
  double maxPhase = 180;
  double phaseStep = 20;
  if (kFALSE == varyPhase) {
    minPhase = phase->getVal();
    maxPhase = phase->getVal();
  }

  for (double x = minPhase; x<=maxPhase; x+=phaseStep) {
    phase->setVal(x);

    double lumiRatio = (double)nEvts / nEvtsNorm;

    RooAbsData * theData = 0;
    double chi2spacing = 1;

    if(signalOnly) {
      theData = pdf.generate(nEvts);
    }
    else {
      RooDataSet * tempData = ((BdkPdfOnRes&)pdf).generate();
      theData = tempData->reduce(cutMinus);
      cout << "Generated " << theData->numEntries() << " events" << endl;
      lumiRatio = 1;
      chi2spacing = 1./9.; // 1/3 sigma
    }
  
    // The nllvar:
    RooNLLVar nll("nll","-log(likelihood)", *(pdf.getPdf()), *theData, RooArgSet());
    
    if (varyPhase) c->cd(++i);

    // The histogram's name:
    TString prefix = theName;
    if (varyPhase) {
      prefix += "_phase=";
      char phaseString[100];
      sprintf(phaseString, "%0.0f", x);
      prefix += phaseString;
    }

    // Make the histogram (using plot function from plotNLL.cc):
    histRhoPhases2 = plot(prefix, prefix, nBins, min, max, nll, *xN, *yN, 
			  orig, kFALSE);
    
    // Set the contour levels, scaled to nEevtsNorm, regardless of how
    // many nEvts we actually generate:
    const int LENGTH = 20;
    double levels[LENGTH];
    for (int l = 0; l < LENGTH; ++l) {
      levels[l] = lumiRatio * l * l * chi2spacing;
    }

    histRhoPhases2->SetContour(LENGTH, levels);
    histRhoPhases2->Draw(option);

    c->Update();
  }


  c->SaveAs(prefix + ".eps");
  c->SaveAs(prefix + ".root");
}



