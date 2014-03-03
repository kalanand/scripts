void IsospinErrorMatrix(){
  const double DEGTORAD = TMath::Pi()/180.0;
  gROOT->Reset();
  gStyle->Reset();
  gROOT->SetStyle("BABAR");
  int numBINS = 30;
  bool doNorm = true;

  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0}) [GeV^{2}/c^{4}]",0,3);
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]",0,3);
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+}) [GeV^{2}/c^{4}]",0,3); 
  fit_S23.setBins(numBINS);
  fit_S31.setBins(numBINS);
  fit_S12.setBins(numBINS);

  BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
  dalitzCfg->getParameters(RooArgSet())->readFromFile("../BToDKTo3piK/Dstar/dalitzCfg.par");
 

  // define signal pdf
  BdkPdfDDalitz dstar("dstar","dstar",fit_S23,fit_S31,BdkDalitzBase::D0);
  RooArgSet allPars = dstar.parameters();
  dstar.parameters().readFromFile("PPP_NominalFit.par");

  if(doNorm==true) {               
    dstar.parameters(); 
    BdkDDalitzAmp::normalizeAll(); 
  }


  RooArgSet fractions = (RooArgSet&) dstar.IntegralOverIsospin();
  RooRealVar* rrv_B_NR = fractions.find("isospin-norm-NR");
  RooRealVar* rrv_B_rhoP = fractions.find("isospin-norm-rhoP");
  RooRealVar* rrv_B_rho0 = fractions.find("isospin-norm-rho0");
  RooRealVar* rrv_B_rhoM = fractions.find("isospin-norm-rhoM");
  RooRealVar* rrv_B_F = fractions.find("isospin-norm-F");
  RooRealVar* rrv_phi_NR = fractions.find("isospin-phase-NR");
  RooRealVar* rrv_phi_rhoP = fractions.find("isospin-phase-rhoP");
  RooRealVar* rrv_phi_rho0 = fractions.find("isospin-phase-rho0");
  RooRealVar* rrv_phi_rhoM = fractions.find("isospin-phase-rhoM");
  RooRealVar* rrv_phi_F = fractions.find("isospin-phase-F");




  // Now Calculate B_{I, 0} coefficients

  double B_NR   = rrv_B_NR->getVal();
  double B_rhoP = rrv_B_rhoP->getVal();
  double B_rho0 = rrv_B_rho0->getVal();
  double B_rhoM = rrv_B_rhoM->getVal();
  double B_F    = rrv_B_F->getVal();
  double phi_NR   = rrv_phi_NR->getVal();
  double phi_rhoP = rrv_phi_rhoP->getVal();
  double phi_rho0 = rrv_phi_rho0->getVal();
  double phi_rhoM = rrv_phi_rhoM->getVal();
  double phi_F    = rrv_phi_F->getVal();

  RooComplex B_NRComplex   = RooComplex(B_NR * cos(phi_NR*DEGTORAD), 
					  B_NR * sin(phi_NR*DEGTORAD));
  RooComplex B_rhoPComplex = RooComplex(B_rhoP * cos(phi_rhoP*DEGTORAD), 
					  B_rhoP * sin(phi_rhoP*DEGTORAD));
  RooComplex B_rho0Complex = RooComplex(B_rho0 * cos(phi_rho0*DEGTORAD), 
					  B_rho0 * sin(phi_rho0*DEGTORAD));
  RooComplex B_rhoMComplex = RooComplex(B_rhoM * cos(phi_rhoM*DEGTORAD), 
					  B_rhoM * sin(phi_rhoM*DEGTORAD));
  RooComplex B_FComplex    = RooComplex(B_F * cos(phi_F*DEGTORAD), 
					  B_F * sin(phi_F*DEGTORAD));
  RooComplex Complex2 = RooComplex(2.0,0.0);


  // Now Calculate C_{I(I12),0} coefficients

  RooComplex C_32Complex = RooComplex(sqrt(10.0)/3.0, 0.0) * B_NRComplex;
  RooComplex C_21Complex = RooComplex(1.0/sqrt(6.0), 0.0)*
    (B_rhoPComplex - Complex2*B_rho0Complex + B_rhoMComplex);
  RooComplex C_11Complex = RooComplex(1.0/sqrt(2.0), 0.0)*(B_rhoPComplex - B_rhoMComplex);
  RooComplex C_01Complex = RooComplex(1.0/sqrt(3.0), 0.0)*
    (B_rhoPComplex + B_rho0Complex + B_rhoMComplex);
  RooComplex C_10Complex = RooComplex(sqrt(3.0)/sqrt(2.0), 0.0)*B_FComplex + 
    RooComplex(5.0/(3.0*sqrt(3.0)), 0.0) * B_NRComplex;

  
  double C_32Mag = C_32Complex.abs();
  double C_21Mag = C_21Complex.abs();
  double C_11Mag = C_11Complex.abs();
  double C_01Mag = C_01Complex.abs();
  double C_10Mag = C_10Complex.abs();
  
  double C_32Phas = TMath::ATan(C_32Complex.im()/C_32Complex.re())/DEGTORAD;
  double C_21Phas = TMath::ATan(C_21Complex.im()/C_21Complex.re())/DEGTORAD;
  double C_11Phas = TMath::ATan(C_11Complex.im()/C_11Complex.re())/DEGTORAD;
  double C_01Phas = TMath::ATan(C_01Complex.im()/C_01Complex.re())/DEGTORAD;
  double C_10Phas = TMath::ATan(C_10Complex.im()/C_10Complex.re())/DEGTORAD;
  


  TString parFile = "cov.par";
  TList list;
  list.AddAll(doPdf(dstar,parFile));


  cout << "***************************************************" << endl;
  cout << "***************************************************" << endl;
  cout << "Original Parameters Coefficients (Eq.16 of isospin note ):"<< endl;
  cout << "|B_NR|     =  " << B_NR      <<endl;
  cout << "|B_rhoP|   =  " << B_rhoP    <<endl;
  cout << "|B_rho0|   =  " << B_rho0    <<endl;
  cout << "|B_rhoM|   =  " << B_rhoM    <<endl;
  cout << "|B_F|      =  " << B_F       <<endl;
  cout << "phi_NR     =  " << phi_NR    <<endl;
  cout << "phi_rhoP   =  " << phi_rhoP  <<endl;
  cout << "phi_rho0   =  " << phi_rho0  <<endl;
  cout << "phi_rhoM   =  " << phi_rhoM  <<endl;
  cout << "phi_F      =  " << phi_F     <<endl;
  cout << "----------------------------------------------------" << endl;
  cout << "Original Isospin Coefficients (Eq.22 of isospin note):"<< endl;
  cout << "|C_32|     =  " <<  C_32Mag  <<endl;
  cout << "|C_21|     =  " <<  C_21Mag  <<endl;
  cout << "|C_11|     =  " <<  C_11Mag  <<endl;
  cout << "|C_01|     =  " <<  C_01Mag  <<endl;
  cout << "|C_10|     =  " <<  C_10Mag  <<endl;
  cout << "phi_C_32   =  " <<  C_32Phas <<endl;
  cout << "phi_C_21   =  " <<  C_21Phas <<endl;
  cout << "phi_C_11   =  " <<  C_11Phas <<endl;
  cout << "phi_C_01   =  " <<  C_01Phas <<endl;
  cout << "phi_C_10   =  " <<  C_10Phas <<endl;
  cout << "***************************************************" << endl;
  cout << "***************************************************" << endl;



  cout << list.GetEntries() << " different parameter configurations."<<endl;


  double BNR[60];
  double BrhoP[60];
  double Brho0[60];
  double BrhoM[60];
  double BF[60];
  double phiNR[60];
  double phirhoP[60];
  double phirho0[60];
  double phirhoM[60];
  double phiF[60];

  RooComplex BNRComplex[60];
  RooComplex BrhoPComplex[60];
  RooComplex Brho0Complex[60];
  RooComplex BrhoMComplex[60];
  RooComplex BFComplex[60];
  RooComplex C32Complex[60];
  RooComplex C21Complex[60];
  RooComplex C11Complex[60];
  RooComplex C01Complex[60];
  RooComplex C10Complex[60];

  double C32Mag[60];
  double C21Mag[60];
  double C11Mag[60];
  double C01Mag[60];
  double C10Mag[60];
  double C32Phas[60];
  double C21Phas[60];
  double C11Phas[60];
  double C01Phas[60];
  double C10Phas[60];

  double z[10][60];

  for (int i=0; i<list.GetEntries(); i++) {
  
    RooArgList* newpars = (RooArgList*)list.At(i);
    cout << "----------------------------------------------------------------"<<endl;
    cout << "    New Pars "<< i << ": " << newpars->GetName() << endl;
    cout << "----------------------------------------------------------------"<<endl;
    newpars->Print("v");


    setDalitzParam(dstar, *newpars);
    RooArgSet fract = (RooArgSet&) dstar.IntegralOverIsospin();
    RooRealVar* rrv_B_NR     = (RooRealVar*) fract.find("isospin-norm-NR");
    RooRealVar* rrv_B_rhoP   = (RooRealVar*) fract.find("isospin-norm-rhoP");
    RooRealVar* rrv_B_rho0   = (RooRealVar*) fract.find("isospin-norm-rho0");
    RooRealVar* rrv_B_rhoM   = (RooRealVar*) fract.find("isospin-norm-rhoM");
    RooRealVar* rrv_B_F      = (RooRealVar*) fract.find("isospin-norm-F");
    RooRealVar* rrv_phi_NR   = (RooRealVar*) fract.find("isospin-phase-NR");
    RooRealVar* rrv_phi_rhoP = (RooRealVar*) fract.find("isospin-phase-rhoP");
    RooRealVar* rrv_phi_rho0 = (RooRealVar*) fract.find("isospin-phase-rho0");
    RooRealVar* rrv_phi_rhoM = (RooRealVar*) fract.find("isospin-phase-rhoM");
    RooRealVar* rrv_phi_F    = (RooRealVar*) fract.find("isospin-phase-F");


    // Calculate B_{I, 0}, C_{I(I12),0} coefficients
    BNR[i]   = rrv_B_NR->getVal();
    BrhoP[i] = rrv_B_rhoP->getVal();
    Brho0[i] = rrv_B_rho0->getVal();
    BrhoM[i] = rrv_B_rhoM->getVal();
    BF[i]    = rrv_B_F->getVal();
    phiNR[i]   = rrv_phi_NR->getVal();
    phirhoP[i] = rrv_phi_rhoP->getVal();
    phirho0[i] = rrv_phi_rho0->getVal();
    phirhoM[i] = rrv_phi_rhoM->getVal();
    phiF[i]    = rrv_phi_F->getVal();

    BNRComplex[i]   = RooComplex(BNR[i] * cos(phiNR[i]*DEGTORAD), 
				 BNR[i] * sin(phiNR[i]*DEGTORAD));
    BrhoPComplex[i] = RooComplex(BrhoP[i] * cos(phirhoP[i]*DEGTORAD), 
				 BrhoP[i] * sin(phirhoP[i]*DEGTORAD));
    Brho0Complex[i] = RooComplex(Brho0[i] * cos(phirho0[i]*DEGTORAD), 
				 Brho0[i] * sin(phirho0[i]*DEGTORAD));
    BrhoMComplex[i] = RooComplex(BrhoM[i] * cos(phirhoM[i]*DEGTORAD), 
				 BrhoM[i] * sin(phirhoM[i]*DEGTORAD));
    BFComplex[i]    = RooComplex(BF[i] * cos(phiF[i]*DEGTORAD), 
				 BF[i] * sin(phiF[i]*DEGTORAD));

    C32Complex[i] = RooComplex(sqrt(10.0)/3.0, 0.0) * BNRComplex[i];
    C21Complex[i] = RooComplex(1.0/sqrt(6.0), 0.0)*
      (BrhoPComplex[i] - Complex2*Brho0Complex[i] + BrhoMComplex[i]);
    C11Complex[i] = RooComplex(1.0/sqrt(2.0), 0.0)*(BrhoPComplex[i] - BrhoMComplex[i]);
    C01Complex[i] = RooComplex(1.0/sqrt(3.0), 0.0)*
      (BrhoPComplex[i] + Brho0Complex[i] + BrhoMComplex[i]);
    C10Complex[i] = RooComplex(sqrt(3.0)/sqrt(2.0), 0.0)*BFComplex[i] + 
      RooComplex(5.0/(3.0*sqrt(3.0)), 0.0) * BNRComplex[i];


    C32Mag[i] = C32Complex[i].abs();
    C21Mag[i] = C21Complex[i].abs();
    C11Mag[i] = C11Complex[i].abs();
    C01Mag[i] = C01Complex[i].abs();
    C10Mag[i] = C10Complex[i].abs();

    C32Phas[i] = TMath::ATan(C32Complex[i].im()/C32Complex[i].re())/DEGTORAD;
    C21Phas[i] = TMath::ATan(C21Complex[i].im()/C21Complex[i].re())/DEGTORAD;
    C11Phas[i] = TMath::ATan(C11Complex[i].im()/C11Complex[i].re())/DEGTORAD;
    C01Phas[i] = TMath::ATan(C01Complex[i].im()/C01Complex[i].re())/DEGTORAD;
    C10Phas[i] = TMath::ATan(C10Complex[i].im()/C10Complex[i].re())/DEGTORAD;



    cout << "***************************************************" << endl;
    cout << "Isospin Coefficients : " << endl;
    cout << "|B_NR["   << i << "]|   =  " << BNR[i]     <<endl;
    cout << "|B_rhoP[" << i << "]|   =  " << BrhoP[i]   <<endl;
    cout << "|B_rho0[" << i << "]|   =  " << Brho0[i]   <<endl;
    cout << "|B_rhoM[" << i << "]|   =  " << BrhoM[i]   <<endl;
    cout << "|B_F["    << i << "]|   =  " << BF[i]      <<endl;
    cout << "phi_NR["   << i << "]   =  " << phiNR[i]   <<endl;
    cout << "phi_rhoP[" << i << "]   =  " << phirhoP[i] <<endl;
    cout << "phi_rho0[" << i << "]   =  " << phirho0[i] <<endl;
    cout << "phi_rhoM[" << i << "]   =  " << phirhoM[i] <<endl;
    cout << "phi_F["    << i << "]   =  " << phiF[i]    <<endl;
    cout << "|C_32["   << i << "]|   =  " << C32Mag[i]  <<endl;
    cout << "|C_21["   << i << "]|   =  " << C21Mag[i]  <<endl;
    cout << "|C_11["   << i << "]|   =  " << C11Mag[i]  <<endl;
    cout << "|C_01["   << i << "]|   =  " << C01Mag[i]  <<endl;
    cout << "|C_10["   << i << "]|   =  " << C10Mag[i]  <<endl;
    cout << "phi_C_32[" << i << "]   =  " << C32Phas[i] <<endl;
    cout << "phi_C_21[" << i << "]   =  " << C21Phas[i] <<endl;
    cout << "phi_C_11[" << i << "]   =  " << C11Phas[i] <<endl;
    cout << "phi_C_01[" << i << "]   =  " << C01Phas[i] <<endl;
    cout << "phi_C_10[" << i << "]   =  " << C10Phas[i] <<endl;
    cout << "***************************************************" << endl; 
 
    z[0][i] = C32Mag[i] - C_32Mag;
    z[1][i] = C32Phas[i] - C_32Phas;
    z[2][i] = C21Mag[i] - C_21Mag;
    z[3][i] = C21Phas[i] - C_21Phas;
    z[4][i] = C11Mag[i] - C_11Mag;
    z[5][i] = C11Phas[i] - C_11Phas;
    z[6][i] = C01Mag[i] - C_01Mag;
    z[7][i] = C01Phas[i] - C_01Phas;
    z[8][i] = C10Mag[i] - C_10Mag;
    z[9][i] = C10Phas[i] - C_10Phas;
  }

  // Now calculate correlation coefficients V_mn between the isospin vectors
  // C_m and C_n
  // index: 0 == |C32|,   2 == |C21|,    4 == |C11|,    6 == |C01|,    8 == |C10|
  //        1 == phi_C32, 3 == phi_C21,  5 == phi_C11,  7 == phi_C01,  9 == phi_C10

  TMatrixDSym mtx(10);
  TMatrixDSym mtxScaled(10);
  for(int i=0; i<10; i++) {
    for(int j=0; j<=i; j++) {
	mtx(i,j) = getErrCoeff(z[i], z[j]);
	mtx(j,i) = mtx(i,j);
    }
  }


  double Icoeff[10];
  Icoeff[0] = C_32Mag;
  Icoeff[1] = C_32Phas;
  Icoeff[2] = C_21Mag;
  Icoeff[3] = C_21Phas;
  Icoeff[4] = C_11Mag;
  Icoeff[5] = C_11Phas;
  Icoeff[6] = C_01Mag;
  Icoeff[7] = C_01Phas;
  Icoeff[8] = C_10Mag;
  Icoeff[9] = C_10Phas;

  double xi0, xj0;
  double result = 0.0;

  for(int i=0; i<10; i++) {

    if(i%2==0) xi0 = Icoeff[6];
    else xi0 = Icoeff[7];
    for(int j=0; j<=i; j++) {
      if(j%2==0) xj0 = Icoeff[6];
      else xj0 = Icoeff[7];

      result = 0.0;

      for(int m=0; m<10; m++) {
	for(int n=0; n<10; n++) {
	  
	  result += PartialDerivative(i, m, Icoeff[i], xi0)* 
	    PartialDerivative(j, n, Icoeff[j], xj0) * mtx[m][n];    
	}
      }
      mtxScaled(i,j) = result;
      mtxScaled(j,i) = mtxScaled(i,j);
    }
  }

  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "Isospin Error matrix : " << endl;
  cout << "index: 0 == |C32|,   2 == |C21|,    4 == |C11|,    6 == |C01|,    8 == |C10|" << endl;
  cout << "index: 1 == phi_C32, 3 == phi_C21,  5 == phi_C11,  7 == phi_C01,  9 == phi_C10" 
       << endl;
  mtx.Print();
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

  cout << "Normalized Isospin Error matrix : " << endl;
  mtxScaled.Print();
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;



  double tt1=0.0, tt2=0.0;
  TMatrixDSym corrmtx(10);

  for(int i=0; i<10; i++) {
    for(int j=0; j<=i; j++) {
      tt1 = mtxScaled(i,i);
      tt2 = mtxScaled(j,j);
      if( tt1==0.0 || tt2==0.0) corrmtx(i,j) = 0.0;
      else corrmtx(i,j) = mtxScaled(i,j)/sqrt(tt1*tt2);
      corrmtx(j,i) = corrmtx(i,j);
    }
  }

  cout << "Normalized Isospin Correlation matrix : " << endl;
  corrmtx.Print();
}




double PartialDerivative(int indexY, int indexX, double xi, double x0) {

  double dydx = 0.0;

  if(indexY%2==0 && !(indexY==6)) {
    if(indexX==6) dydx = -xi/(x0*x0);
    if(!(indexX==6) && indexX==indexY) dydx = 1.0/x0;
  }

  if(indexY%2==1 && !(indexY==7)) {
    if(indexX==7) dydx = -1.0;
    if(!(indexX==7) && indexX==indexY) dydx = 1.0;
  }

  return dydx;
}





// isospin error coefficients
//  V_mn = (1/N) *  (1/2) * Sum{ (x_+) * (y_+)  +  (x_-) * (y_-) }
double getErrCoeff(double x[], double y[]) {

  double coeff = 0.0;

  for(int k=0; k<60; k++) {
    
    if(k%2==1) continue;    
    coeff +=  x[k]   * y[k];  
    coeff +=  x[k+1] * y[k+1];
  }	

  return(coeff/60.0);  

}



// Set the pdf parameters to the new parameters 'newpars'
void setDalitzParam(const BdkPdfDDalitz& pdf, const RooArgList& pars)
{
  TIterator* iter = pars.createIterator();
  RooRealVar* r;
  while (r = (RooRealVar*)iter->Next()) {
    TString name(r->GetName());
    RooRealVar* rP = (RooRealVar*)pdf.parameters().find(name);
    if (rP) rP->setVal(r->getVal());    
  }
}





// Do the systematics for one pdf
TList* doPdf(const BdkPdfDDalitz& pdf, const char* parFile)
{
 
  RooArgList pars;
  RooArgSet allPars = pdf.parameters();  
  RooRealVar* var[30];
  TString strg = pdf.GetName() + TString(".pdf.dalitzAmp.");
  var[0]  = (RooRealVar*) allPars.find( strg + TString("F0_1370_amp") );
  var[1]  = (RooRealVar*) allPars.find( strg + TString("F0_1370_phase") );
  var[2]  = (RooRealVar*) allPars.find( strg + TString("F0_1500_amp") );
  var[3]  = (RooRealVar*) allPars.find( strg + TString("F0_1500_phase") );
  var[4]  = (RooRealVar*) allPars.find( strg + TString("F0_1710_amp") );
  var[5]  = (RooRealVar*) allPars.find( strg + TString("F0_1710_phase") );
  var[6]  = (RooRealVar*) allPars.find( strg + TString("F0_amp") );
  var[7]  = (RooRealVar*) allPars.find( strg + TString("F0_phase") );
  var[8]  = (RooRealVar*) allPars.find( strg + TString("F2_amp") );
  var[9]  = (RooRealVar*) allPars.find( strg + TString("F2_phase") );
  var[10] = (RooRealVar*) allPars.find( strg + TString("Nonres_amp") );
  var[11] = (RooRealVar*) allPars.find( strg + TString("Nonres_phase") );
  var[12] = (RooRealVar*) allPars.find( strg + TString("Rho-_amp") );
  var[13] = (RooRealVar*) allPars.find( strg + TString("Rho-_phase") );
  var[14] = (RooRealVar*) allPars.find( strg + TString("Rho0_amp") );
  var[15] = (RooRealVar*) allPars.find( strg + TString("Rho0_phase") );
  var[16] = (RooRealVar*) allPars.find( strg + TString("Rho1700+_amp") );
  var[17] = (RooRealVar*) allPars.find( strg + TString("Rho1700+_phase") );
  var[18] = (RooRealVar*) allPars.find( strg + TString("Rho1700-_amp") );
  var[19] = (RooRealVar*) allPars.find( strg + TString("Rho1700-_phase") );
  var[20] = (RooRealVar*) allPars.find( strg + TString("Rho17000_amp") );
  var[21] = (RooRealVar*) allPars.find( strg + TString("Rho17000_phase") );
  var[22] = (RooRealVar*) allPars.find( strg + TString("Rho2s+_amp") );
  var[23] = (RooRealVar*) allPars.find( strg + TString("Rho2s+_phase") );
  var[24] = (RooRealVar*) allPars.find( strg + TString("Rho2s-_amp") );
  var[25] = (RooRealVar*) allPars.find( strg + TString("Rho2s-_phase") );
  var[26] = (RooRealVar*) allPars.find( strg + TString("Rho2s0_amp") );
  var[27] = (RooRealVar*) allPars.find( strg + TString("Rho2s0_phase") );
  var[28] = (RooRealVar*) allPars.find( strg + TString("Sigma_amp") );
  var[29] = (RooRealVar*) allPars.find( strg + TString("Sigma_phase") );

  for(int i=0; i<30; i++) { pars.add(*var[i]); }

  cout << "-----------------------------------------"<<endl;
  cout << " Original isobar-model fit parameters : "<< endl;
  cout << "-----------------------------------------"<<endl;
  pars.Print("v");


  TMatrixDSym m = getErrorMatrix(pdf, parFile);
     
  TList* list =  pdfChangePars(pars, m);
  if (list==0) {
    cout << "Nothing to do for "<<pdf.GetName()<<endl;
    return;
  }
  
  //Correlation matrix
  TMatrixD mcorr = m;
  for (int i=0; i<mcorr.GetNcols(); i++) {
    for (int j=0; j<mcorr.GetNrows(); j++) {
      mcorr(i,j) /= ((RooRealVar&)pars[i]).getError()*
                    ((RooRealVar&)pars[j]).getError();
    }
  }

//   cout << endl << "Error matrix for "<<pdf.GetName()<<endl;
//   m.Print();
//   cout << endl << "Correlation matrix for "<<pdf.GetName()<<endl;
//   mcorr.Print();

  return list;
}






// Returns a TList of RooArgLists with parameters varied by +- sigmas
// The first RooArgList has the first parameter varied by +sigmas
// The second     "     "    "     "     "         "   by -sigmas  etc.
TList* pdfChangePars(const RooArgList& pars, 
                     const TMatrixDSym& m, 
                     Double_t sigmas = 1) 
{
  if (m.GetNoElements()==0) {
    cout << "Empty error matrix."<<endl;
    return 0;
  }


  if (pars.getSize()!=m.GetNrows()) {
    cout << "Number of parameters does not match dimensions of matrix."<<endl;
    return 0;
  }

  // This corrects a mistake in our par files (always has 0 if one floating parameter)
  if (m.GetNoElements()==1) {
    m(0,0) = ((RooRealVar&)pars[0]).getError();
    m.Sqr();
  }


  TVectorD p(pars.getSize());

  for (int i=0; i<pars.getSize(); i++) {
    p(i) = ((RooRealVar*)pars[i])->getVal();
  }


  // Get eigenvalues e and transformation matrix V
  TMatrixDSymEigen eigen(m);
  TMatrixD V = eigen.GetEigenVectors();

  // inverse transformation matrix
  TMatrixD Vinv = V;
  Vinv.InvertFast();

  // Get eigenvalue matrix
  TMatrixD diag = Vinv*m*V;

  TVectorD pe = Vinv*p;   // Vector p in eigen-basis

  //  p.Print();
  //  pe.Print();
  //  V.Print();
  //  diag.Print();

  TList* list = new TList();
  
  for (int i=0; i<pe.GetNrows(); i++) {
    TVectorD pPlus = pe;
    TVectorD pMinus = pe;
    pPlus(i) += sigmas*sqrt(diag(i,i));
    pMinus(i) -= sigmas*sqrt(diag(i,i));
  
    // back into original basis
    pPlus = V*pPlus;
    pMinus = V*pMinus;

    // create new RooArgLists with changes parameters
    RooArgList* newpars = (RooArgList*)pars.snapshot(false);
    newpars->setName(pars[i].GetName()+TString("_Plus"));
    list->Add(newpars);
    copyPars(*newpars, pPlus);
    
    newpars = (RooArgList*)pars.snapshot(false);
    newpars->setName(pars[i].GetName()+TString("_minus"));
    list->Add(newpars);
    copyPars(*newpars, pMinus);
  }
  return list;
}







// copy p to list
void copyPars(RooArgList& list, const TVectorD& p)
{
  for (int j=0; j<p.GetNrows(); j++) ((RooRealVar&)list[j]).setVal(p(j));
}









// Read error matrix for pdf from parFile.
// Returns 0x0 matrix if none is found.
TMatrixDSym getErrorMatrix(const BdkPdfAbsBase& pdf, const char* parFile)
{
  const Bool_t DEBUG = false;

  // State machine
  enum state {SM_READ, SM_SEEK, SM_MATRIX, SM_UNCORR, SM_STOP};
  state sm;

  ifstream if;
  if.open(parFile);

  // Run the state machine to find the error matrix
  ostringstream matrix;
  sm = SM_SEEK;
  while (!if.eof() && sm!=SM_STOP && sm!=SM_UNCORR) {
    string line;
    getline(if,line);
    TString s(line.c_str());

    if (sm==SM_READ) {
      if (s.Contains(pdf.GetName())) sm = SM_SEEK;
    }
    else if (sm==SM_SEEK) {
      if (s.Contains("//$")) sm = SM_MATRIX;
      else if (s.Contains("//") || s.Contains(pdf.GetName()) || s=="") sm = SM_SEEK;
      else sm = SM_READ;
    }
    else if (sm==SM_MATRIX) {
      if (s.Contains("//$")) sm = SM_STOP;
      else if (s.Contains("UNCORRELATED")) sm = SM_UNCORR;
      else matrix << s << endl;
    }
  }
  if.close();

  if (DEBUG) {
    cout << "State machine state = " << sm << endl;
    cout << "Error matrix for "<<pdf.GetName()<<endl;
    if (matrix.str()=="") cout << "Not found."<<endl;
    else cout << matrix.str() << endl;
  }

  if (sm==SM_UNCORR) {
    RooArgList pars(pdf.parametersFree());
    TMatrixDSym m(pars.getSize());
    for (int i=0; i<pars.getSize(); i++) {
      m(i,i) = sqr(((RooRealVar*)pars[i])->getError());
    }
    return m;
  }
  else return parseErrorMatrix(matrix.str());
}








// Parse error matrix from string
TMatrixDSym parseErrorMatrix(string matrix)
{
  TMatrixDSym m(0);
  if (matrix=="") return m;

  istringstream is(matrix);
  int j = 0;
  while (!is.eof()) {
    string line;
    getline(is,line);
    TString s(line);
    TObjArray* a = s.Tokenize(" ");
    int  N = a->GetEntries()-1;

    if (m.GetNoElements()==0) m.ResizeTo(N,N);
    // Format is: "// 1 2 3 4"
    for (int i=1; i<a->GetEntries(); i++) {  
      m(i-1,j) = atof(a->At(i)->GetName());
    }
    j++;
  }

  return m;
}
