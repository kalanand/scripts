// $Id: mcStatSyst.cc,v 1.9 2007/04/02 14:15:34 fwinkl Exp $
// Systematic error due to MC statistics
// run with syst/mcStatSyst.csh for all event types

void mcStatSyst(const BdkEvtTypes::Type t,
                const char* resultDir = "./",
                Bool_t fit = kTRUE,
		  TString parFile = "all.par",
                Bool_t recalcDDbarNorm = kFALSE)
{
  cout << "-----------------------------------------------------"<<endl;
  cout << "     MC statistics error for event type "<<t<<endl;
  cout << "-----------------------------------------------------"<<endl;

  // Read parameters with original setting of fixed/floating
  useBothFitVars();
  readOnResDKPar(false);
  pdfOnResDK.parameters().readFromFile(parFile);
  
  fixAll(pdfOnResDK.cpParams());
  fixVar(pdfOnResDK.parameters(),"blindMode");
    
  if (t<0 || t>=BdkEvtTypes::NTYPES) {
    cout << "Invalid event type "<<t<<endl;
    return;
  }

  TList list;
  list.AddAll(doPdf(*pdfOnResDK.prodN(t)->getVarPdf(BdkPdfProdAll::DELTAE),parFile));
  list.AddAll(doPdf(*pdfOnResDK.prodN(t)->getVarPdf(BdkPdfProdAll::NNCONT),parFile));
  list.AddAll(doPdf(*pdfOnResDK.prodN(t)->getVarPdf(BdkPdfProdAll::NNCOMB),parFile));
  list.AddAll(doPdf(*pdfOnResDK.prodN(t)->getVarPdf(BdkPdfProdAll::DALITZ),parFile));

  cout << list.GetEntries() << " different parameter configurations."<<endl;

  readCut = cutSigReg;
  if (fit) data = read(dataTree);
  
  // Read parameters with only numEvts.par floating
  readOnResDKPar();
  RooArgSet* initParams = (RooArgSet*)pdfOnResDK.parameters().snapshot(false);
  initParams->setName("initParams");
  initParams->Print("v");

  TList fitResults;

  // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with original parameters"<<endl;
  cout << "-----------------------------------------------------"<<endl;

  if (fit) {
    fit(pdfOnResDK, *data, true);
    RooArgSet* result = pdfOnResDK.fitResult();
    result->setName("ref");  
    pdfOnResDK.yieldFitResult()->Print();
    pdfOnResDK.xyFitResult()->Print();
    fitResults.Add(result);
  }

  for (int i=0; i<list.GetEntries(); i++) {
    // reset initial parameters
    useBothFitVars();
    pdfOnResDK.parameters() = *initParams;

    RooArgList* newpars = (RooArgList*)list.At(i);
    cout << "----------------------------------------------------------------"<<endl;
    cout << "    Fit "<< i << ": " << newpars->GetName() << endl;
    cout << "----------------------------------------------------------------"<<endl;
    pdfOnResDK.parameters().selectCommon(*newpars)->Print("v");
    newpars->Print("v");
  
    if (fit) {
      // set new parameters
      pdfOnResDK.parameters() = *newpars;
      setDalitzP(*newpars);

      // This needs to be done if amplitudes and phases are being changed
      if (recalcDDbarNorm) {
	dalitzHolderN.sigGoodD0Type().dalitzAmp()->calDDbarNorm(1e7);
	// Calculate and show what the new x0 is (but don't use it)
	Double_t newx0 = -dalitzHolderN.sigGoodD0Type().dalitzAmp()->normReDDbar() / 
	  dalitzHolderN.sigGoodD0Type().dalitzAmp()->normDSqr();
	cout << "Recalculated x0 = " << newx0 << endl;
      }
      
      // fit
      fit(pdfOnResDK, *data, true);
      pdfOnResDK.yieldFitResult()->Print();
      pdfOnResDK.xyFitResult()->Print();
      
      RooArgSet* result = pdfOnResDK.fitResult();
      result->setName(newpars->GetName());
      fitResults.Add(result);    
    }
  }
    
  // Save all fit results
  TString filename = TString(resultDir)+"/mcStatSyst-";
  filename += (int)t;
  filename += ".root";

  TFile f(filename,"recreate");
  fitResults.Write();
  f.Close();
}


// Do the systematics for one pdf
TList* doPdf(const BdkPdfAbsBase& pdf, const char* parFile)
{
  RooArgList pars(pdf.parametersFree());
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

  cout << endl << "Error matrix for "<<pdf.GetName()<<endl;
  m.Print();
  cout << endl << "Correlation matrix for "<<pdf.GetName()<<endl;
  mcorr.Print();

  return list;
}


// Set dalitzP parameters from dalitzN parameters in pars
void setDalitzP(const RooArgList& pars)
{
  TIterator* iter = pars.createIterator();
  RooRealVar* r;
  while (r = (RooRealVar*)iter->Next()) {
    TString name(r->GetName());
    if (name.Contains("dalitzHolderN")) {
      name.ReplaceAll("dalitzHolderN","dalitzHolderP");
      RooRealVar* rP = (RooRealVar*)pdfOnResDK.parameters().find(name);
      if (rP) rP->setVal(r->getVal());
    }
  }
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

  /*
  // For D Dalitz systematics START
  // Set some eror matrix elements to zero  
  int setZero[] = {0,2,4,8,16,18,20,22,24,26,28};
  for (int k=0; k<sizeof(setZero)/sizeof(int); k++) {
    cout << "Setting error matrix element to zero for "<<pars.at(setZero[k])->GetName()<<endl;
    for (int i=0; i<m.GetNcols(); i++) {
      m(setZero[k],i) = 0;
      m(i,setZero[k]) = 0;
    }
  }
  // For D Dalitz systematics END
  */
  
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
  sm = SM_READ;
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
