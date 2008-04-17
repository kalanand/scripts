// The function sensitivity() produces plots of the sensitivities to gamma, 
// x and y for B- only. 


TH2F * gammaSens = 0;
TH2F * rhoSens = 0;
TH2F * thetaSens = 0;
TH2  * dataHist = 0;
RooDataSet * wgtData = 0;

// Calculates rho and theta given gamma, for a fixed value of delta and r_B.
// Note that we are running this file for B-.
void calcRhoTheta(double gamma, RooRealVar & rhoM, RooRealVar & thetaM) {
  // for r_B = 0.1 and the world average gamma=62 deg, the BaBar
  // Kspipi analysis yields  delta=102:
  static const double delta = 102;
  static const double r_B = 0.1;

  double phase = delta - gamma;
  rhoM.setVal(rhoFromCP(r_B, phase));
  thetaM.setVal(thetaFromCP(r_B, phase));
}

// Returns the 2nd derivative of the log of the PDF. The formula is 
// [pdf(x+zStep) - 2 pdf(x) + pdf(x-zStep)]/zStep^2.
//
// If positiveTerm=true, then 
// returns the positive term of the 2nd derivative of the log of the PDF.
// The second derivative is a sume over all events of 
// d/dz(d/dz (log P)), where P is the PDF of the event, and z is the 
// variable WRT which we are maximizing the log likelihood. This equals
// - d/dz(P'/P), where P' = dP/dz. This further equals
// (P'/P)^2 - P''/p.
// The _average_ over all events of the 2nd term is 0, since 
// integral{P (P''/P)} = integral{P''} = d^2/dz^2 (integral{P}), and 
// integral{P} = 1 by definition of normalization. So we are left with
// the positive piece of the 2nd derivative being (P'/P)^2:
//
Double_t secondDer(RooAbsPdf & pdf, 
		   Double_t m12Val, Double_t m13Val, 
		   RooRealVar & z, Double_t zStep,
		   Bool_t positiveTerm = kFALSE) {
  
  m12->setVal(m12Val);
  m13->setVal(m13Val);
  const Double_t zOrig = z.getVal();

  if (0 >= pdf.getVal()) {
    return 0;
  }
  
  const Double_t f = log(pdf.getVal());

  z.setVal(zOrig + zStep);
  const Double_t fUp = log(pdf.getVal());
      
  z.setVal(zOrig - zStep);
  const Double_t fDown = log(pdf.getVal());

  double result = 0;
  if (kFALSE == positiveTerm) {// full 2nd derivative of the NLL:
    result = (-1)*(fUp - 2 * f + fDown) / zStep / zStep;
  }
  else { // then return the 1st derivative squared. See comment above:
    const Double_t firstDer = (fUp - fDown) / 2 / zStep;
    result = firstDer * firstDer / f /  f;
  }

  z.setVal(zOrig);   // reset
  return result;
}
      
// Does the same as secondDer, but with respect to gamma:
Double_t secondDerGamma(RooAbsPdf & pdf, 
			Double_t m12Val, Double_t m13Val, 
			Double_t gamma, Double_t zStep,
			RooRealVar & rhoM, RooRealVar & thetaM,
			Bool_t positiveTerm = kFALSE) {
  
  m12->setVal(m12Val);
  m13->setVal(m13Val);

  if (0 >= pdf.getVal()) {
    return 0;
  }
  
  double rhoOrig = rhoM.getVal();
  double thetaOrig = thetaM.getVal();

  // get nominal point:
  calcRhoTheta(gamma, rhoM, thetaM);
  const Double_t f = log(pdf.getVal());

  calcRhoTheta(gamma + zStep, rhoM, thetaM); 
  const Double_t fUp = log(pdf.getVal());
      
  calcRhoTheta(gamma - zStep, rhoM, thetaM); 
  const Double_t fDown = log(pdf.getVal());

  double result = 0;
  if (kFALSE == positiveTerm) {// then calculate the full 2nd derivative:
    result = (-1)*(fUp - 2 * f + fDown) / zStep / zStep;
  }
  else { // then return the 1st derivative squared. See comment above:
    const Double_t firstDer = (fUp - fDown) / 2 / zStep;
    result = firstDer * firstDer / f /  f;
  }

  rhoM.setVal(rhoOrig);   // reset
  thetaM.setVal(thetaOrig);   // reset
  return result;
}
      


// Generates a data set created according to the PDF value, reweighted
// according to the value (or positive value, if positiveTerm=true) of
// the 2nd derivative:

void sensitivity(Bool_t positiveTerm = kFALSE,
		 int nBins = 20,
		 int nEvents = 100000, 
		 double zStep = 1.e-08) { 

  BdkPdfAbsBase & wrapper = dalitzHolderN.sigGoodD0Type();
  RooAbsPdf * pdf = dalitzHolderN.sigGoodD0Type().getPdf();
  RooDataSet * unweighted = wrapper.generate(nEvents);
  RooRealVar * rho = dalitzHolderN.sigGoodD0Type().rho();
  RooRealVar * theta = dalitzHolderN.sigGoodD0Type().theta();

  gammaSens = new TH2F("gammaSens", "gammaSens", nBins, 
		   m12->getMin(), m12->getMax(), 
		   nBins, m13->getMin(), m13->getMax());
  
  rhoSens = new TH2F("rhoSens", "rhoSens", nBins, 
		   m12->getMin(), m12->getMax(), 
		   nBins, m13->getMin(), m13->getMax());
  
  thetaSens = new TH2F("thetaSens", "thetaSens", 
		   nBins, m12->getMin(), m12->getMax(), 
		   nBins, m13->getMin(), m13->getMax());

  dataHist = new TH2F("dataHist", "dataHist", 
		   nBins, m12->getMin(), m12->getMax(), 
		   nBins, m13->getMin(), m13->getMax());

  for (int e = 0; e < unweighted->numEntries(); ++e){
    const RooArgSet * event = unweighted->get(e);
    
    Double_t tempM12 = 
      ((RooAbsReal*)(event->find(m12->GetName())))->getVal();

    Double_t tempM13 = 
      ((RooAbsReal*)(event->find(m13->GetName())))->getVal();

    double gamma = 62;

    double gammaWeight = 
      secondDerGamma(*pdf, tempM12, tempM13, gamma, zStep, *rho, *theta, 
		     positiveTerm);

    double rhoWeight = 
      secondDer(*pdf, tempM12, tempM13, *rho, zStep, positiveTerm);

    double thetaWeight = 
      secondDer(*pdf, tempM12, tempM13, *theta, zStep, positiveTerm);    

    //    cout << " " << tempM12 << " " << tempM13 << " " 
    //    	 << rhoWeight << " " << thetaWeight << endl;
    
    dataHist->Fill(tempM12, tempM13);
    rhoSens->Fill(tempM12, tempM13, rhoWeight);
    thetaSens->Fill(tempM12, tempM13, thetaWeight);
    gammaSens->Fill(tempM12, tempM13, gammaWeight);

  }

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TCanvas * can = new TCanvas("can", "can", 1200, 400);
  can->Divide(3,1);

  can->cd(1);
  rhoSens->Draw("contz");
  can->cd(2);
  thetaSens->Draw("contz");
  can->cd(3);
  gammaSens->Draw("contz");

  can->SaveAs("sensitivity.eps");
}

