// $Id: plotNLL_RhoPhases2.cc,v 1.1 2006/06/26 23:09:16 abi Exp $


#include "../BToDKTo3piK/analysis/plotNLL.cc"

TH2F * histRhoAmpsPhases2 = 0;

void plotNLL_RhoAmpsPhases2(int nBins = 20, double min = -2, double max = 2,
			    const char * option = "cont3",
			    double minRatio = 0, double maxRatio = 1,
			    double ratioStep = 0.2,
			    Bool_t save = kTRUE) {

  // Read the parameters we want to start with:
  pdfOnResDK.parameters().
    readFromFile("../BToDKTo3piK/params/cross/plotNLL_RhoPhases2.par");

  // Get the charged rho amps:
  RooRealVar* ampMinus = (RooRealVar*)pdf.parameters().
    find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_amp");
    
  RooRealVar* ampPlus = (RooRealVar*)pdf.parameters().
    find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho+_amp");
    
  // Loop over ratio of rho- to rho+ amps:
  for (double ratio = minRatio; ratio += ratioStep; ratio < maxRatio) {
    ampMinus = ampPlus->getVal() * ratio;

    



// 2D-Plot the NLL, changing rho- phase. Otherwise takes whatever other
// parameter values are in place:
void plotNLL_RhoPhases2(int nBins = 20, double min = -2, double max = 2,
			const char * option = "cont3",
			const char * name) {
  
  gStyle->SetPalette(1);
  TCanvas* c = new TCanvas(name, name,1500,600);
  c->Divide(5,2);
 
  BdkPdfDKDalitz& pdf = dalitzHolderN.sigGoodD0Type();

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
  for (double x = 0; x<=180; x+=20) {

    phase->setVal(x);
    data = pdf.generate(200);
  
    // The nllvar:
    RooNLLVar nll("nll","-log(likelihood)", *(pdf.getPdf()), *data, RooArgSet());
    
    c->cd(++i);
    TString prefix = name + "_phase=";
    prefix += x;

    histRhoPhases2 = plot(prefix, prefix, nBins, min, max, nll, *xN, *yN, 
			  orig, kFALSE);
    
    double levels[] = {0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 
		       121, 144, 169, 196, 225};
    
    int length = sizeof(levels)/sizeof(double);
    histRhoPhases2->SetContour(length, levels);


    histRhoPhases2->Draw(option);

    delete data;
    c->Update();
  }


  c->SaveAs(prefix + ".eps");
  c->SaveAs(prefix + ".root");
}



