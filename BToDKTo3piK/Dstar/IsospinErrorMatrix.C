void IsospinErrorMatrix(){

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


  // Now Calculate B_{I, 0} coefficients

  double B_NR   = rrv_B_NR->getVal();
  double B_rhoP = rrv_B_rhoP->getVal();
  double B_rho0 = rrv_B_rho0->getVal();
  double B_rhoM = rrv_B_rhoM->getVal();
  double B_F    = rrv_B_F->getVal();


  // Now Calculate C_{I(I12),0} coefficients

  double C_32 = B_NR;
  double C_21 = (1.0/sqrt(6.0))*(B_rhoP - 2.0*B_rho0 + B_rhoM);
  double C_11 = (1.0/sqrt(2.0))*(B_rhoP - B_rhoM);
  double C_01 = (1.0/sqrt(3.0))*(B_rhoP + B_rho0 + B_rhoM);
  double C_10 = B_F;


  TString parFile = "cov.par";
  TList list;
  list.AddAll(doPdf(dstar,parFile));


  cout << "***************************************************" << endl;
  cout << "***************************************************" << endl;
  cout << "Original Parameters Coefficients (Eq.16 of isospin note ):"<< endl;
  cout << "B_NR     =  " << B_NR    <<endl;
  cout << "B_rhoP   =  " << B_rhoP  <<endl;
  cout << "B_rho0   =  " << B_rho0  <<endl;
  cout << "B_rhoM   =  " << B_rhoM  <<endl;
  cout << "B_F      =  " << B_F     <<endl;
  cout << "----------------------------------------------------" << endl;
  cout << "Original Isospin Coefficients (Eq.22 of isospin note):"<< endl;
  cout << "C_32     =  " <<  C_32  <<endl;
  cout << "C_21     =  " <<  C_21  <<endl;
  cout << "C_11     =  " <<  C_11  <<endl;
  cout << "C_01     =  " <<  C_01  <<endl;
  cout << "C_10     =  " <<  C_10  <<endl;
  cout << "***************************************************" << endl;
  cout << "***************************************************" << endl;



  cout << list.GetEntries() << " different parameter configurations."<<endl;


  double BNR[60];
  double BrhoP[60];
  double Brho0[60];
  double BrhoM[60];
  double BF[60];
  double C32[60];
  double C21[60];
  double C11[60];
  double C01[60];
  double C10[60];
  double z[5][60];

  for (int i=0; i<list.GetEntries(); i++) {
  
    RooArgList* newpars = (RooArgList*)list.At(i);
    cout << "----------------------------------------------------------------"<<endl;
    cout << "    New Pars "<< i << ": " << newpars->GetName() << endl;
    cout << "----------------------------------------------------------------"<<endl;
    newpars->Print("v");


    setDalitzParam(dstar, *newpars);
    RooArgSet fract = (RooArgSet&) dstar.IntegralOverIsospin();
    RooRealVar* rrv_B_NR = (RooRealVar*) fract.find("isospin-norm-NR");
    RooRealVar* rrv_B_rhoP = (RooRealVar*) fract.find("isospin-norm-rhoP");
    RooRealVar* rrv_B_rho0 = (RooRealVar*) fract.find("isospin-norm-rho0");
    RooRealVar* rrv_B_rhoM = (RooRealVar*) fract.find("isospin-norm-rhoM");
    RooRealVar* rrv_B_F = (RooRealVar*) fract.find("isospin-norm-F");

    // Calculate B_{I, 0}, C_{I(I12),0} coefficients
    BNR[i]   = rrv_B_NR->getVal();
    BrhoP[i] = rrv_B_rhoP->getVal();
    Brho0[i] = rrv_B_rho0->getVal();
    BrhoM[i] = rrv_B_rhoM->getVal();
    BF[i]    = rrv_B_F->getVal();

    C32[i] = BNR[i];
    C21[i] = (1.0/sqrt(6.0))*(BrhoP[i] - 2.0*Brho0[i] + BrhoM[i]);
    C11[i] = (1.0/sqrt(2.0))*(BrhoP[i] - BrhoM[i]);
    C01[i] = (1.0/sqrt(3.0))*(BrhoP[i] + Brho0[i] + BrhoM[i]);
    C10[i] = BF[i];

    cout << "***************************************************" << endl;
    cout << "Isospin Coefficients : " << endl;
    cout << "B_NR["   << i << "]   =  " << BNR[i]    <<endl;
    cout << "B_rhoP[" << i << "]   =  " << BrhoP[i]  <<endl;
    cout << "B_rho0[" << i << "]   =  " << Brho0[i]  <<endl;
    cout << "B_rhoM[" << i << "]   =  " << BrhoM[i]  <<endl;
    cout << "B_F["    << i << "]   =  " << BF[i]     <<endl;
    cout << "C_32["   << i << "]   =  " <<  C32[i]   <<endl;
    cout << "C_21["   << i << "]   =  " <<  C21[i]   <<endl;
    cout << "C_11["   << i << "]   =  " <<  C11[i]   <<endl;
    cout << "C_01["   << i << "]   =  " <<  C01[i]   <<endl;
    cout << "C_10["   << i << "]   =  " <<  C10[i]   <<endl;
    cout << "***************************************************" << endl; 
 
    z[0][i] = C32[i] - C_32;
    z[1][i] = C21[i] - C_21;
    z[2][i] = C11[i] - C_11;
    z[3][i] = C01[i] - C_01;
    z[4][i] = C10[i] - C_10;
  }

  // Now calculate correlation coefficients V_mn between the isospin vectors
  // C_m and C_n
  // identification: 0 == C_32, 1 == C_21, 2 == C_11, 3 == C_01, 4 == C_10

  TMatrixDSym mtx(5);
  TMatrixDSym mtxScaled(5);
  for(int i=0; i<5; i++) {
    for(int j=0; j<=i; j++) {
	mtx(i,j) = getErrCoeff(z[i], z[j]);
	mtxScaled(i,j) = getErrCoeff(z[i], z[j])/(C_01*C_01);
	mtx(j,i) = mtx(i,j);
	mtxScaled(j,i) = mtxScaled(i,j);
    }
  }

  mtxScaled(3,3) = 0.0; // since this element is fixed
  
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "Isospin Error matrix : " << endl;
  cout << "index: 0==C_32, 1==C_21, 2==C_11, 3==C_01, 4==C_10" << endl;
  mtx.Print();
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

  cout << "Normalized Isospin Error matrix : " << endl;
  mtxScaled.Print();
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
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
