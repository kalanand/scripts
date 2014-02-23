// plots the Dalitz plot distribution of signal events for different 
// values of x, after reading pars from a file:
void plotDPXY(int nEvt = 3000,
	      double deltaX = 1,
	      const char * parFile = 
	      "../BToDKTo3piK/params/cross/plotNLL_RhoPhases2.par") {
  TCanvas * can = new TCanvas("plotDPXY", "plotDPXY", 1500, 1000);
  can->Divide(3,2);
  int cd = 0;

  // Plot 3 plots with parameters from plotNLL_RhoPhases2.par:
  plotDPXYFile(cd, nEvt, deltaX, can, "plotDPXY.eps");

  // Plot 3 plots with current parameters:
  plotDPXYFile(cd, nEvt, deltaX, can, "plotDPXY.eps", "../BToDKTo3piK/params/all.par");
}
  

void plotDPXYFile(int & cd,
		  int nEvt = 3000, 
		  double deltaX = 1,
		  TCanvas * can = 0,
		  const char * epsFile = "plotDPXYFile.eps",
		  const char * parFile = 
		  "../BToDKTo3piK/params/cross/plotNLL_RhoPhases2.par") {

  BdkPdfAbsBase & pdf = dalitzHolderN.sigGoodD0Type();

  if (parFile) {
    cout << "Reading parameters from " << parFile << endl;
    pdf.parameters().readFromFile(parFile);
  }
  else {
    cout << "Using current parameters" << endl;
  }

  char deltaXStringArray[10];
  sprintf(deltaXStringArray, "%0.2f", deltaX);
  TString deltaXString(deltaXStringArray);

  RooRealVar * xN = pdfOnResDK.xMinus();
  RooRealVar * yN = pdfOnResDK.yMinus();
  yN->setVal(0);
  
  if (0 == can) {
    can = new TCanvas("plotDPXYFile", "plotDPXYFile", 1500, 500);
    can->Divide(3,1);
  }

  can->cd(++cd);
  xN->setVal(0);
  TH2F * hist1 = pdf.generate(nEvt)->createHistogram(*m12, *m13);
  hist1->SetTitle("x=0");
  hist1->Draw();

  can->cd(++cd);
  xN->setVal(-deltaX);
  TH2F * hist2 = pdf.generate(nEvt)->createHistogram(*m12, *m13);
  hist2->SetTitle(TString("x=-") + deltaXString);
  hist2->Draw();

  can->cd(++cd);
  xN->setVal(deltaX);
  TH2F * hist3 = pdf.generate(nEvt)->createHistogram(*m12, *m13);
  hist3->SetTitle(TString("x=") + deltaXString);
  hist3->Draw();

  can->SaveAs(epsFile);
}
