// $Id: asymEffSyst.cc,v 1.1 2006/07/11 21:48:51 fwinkl Exp $
// Script to evaluate the systematic error due to an asymmetric efficiency


void asymEffSyst()
{

  const Double_t scale = 100;    // times data luminosity

  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with asymmetric efficiency"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  useBothFitVars();
  readOnResDKPar();
  setNoCP();

  // Read efficiencies with asymmetric coefficients
  eff.parameters().readFromFile("../BToDKTo3piK/params/eff-asym.par");
  effOther.parameters().readFromFile("../BToDKTo3piK/params/eff-asym.par");

  // Generate B-
  BdkDDalitzAmp::normalizeAll(1e7);

  RooDataSet* dataN = genData(-1,scale);

  // Generate B+ with flipped asymmetric coefficients
  flipEffAsym();
  BdkDDalitzAmp::normalizeAll(1e7);
  data = genData(+1,scale);
  
  // merge
  data->append(*dataN);

  pdfOnResDK.nBB()->setVal(scale*pdfOnResDK.nBB()->getVal());
  fit(pdfOnResDK, *data, true);
  RooArgSet* asym = pdfOnResDK.fitResult();
  asym->setName("asym");


  // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with symmetric efficiency"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  useBothFitVars();
  readOnResDKPar();
  setNoCP();

  const char* effPar = "../BToDKTo3piK/params/eff.par";
  eff.parameters().readFromFile(effPar);
  effOther.parameters().readFromFile(effPar);

  pdfOnResDK.nBB()->setVal(scale*pdfOnResDK.nBB()->getVal());
  fit(pdfOnResDK, *data, true);
  RooArgSet* sym = pdfOnResDK.fitResult();
  sym->setName("sym");

  asym->Print("v");
  sym->Print("v");

  removeAllRanges(*sym);
  sub(*sym, *asym);
  sym->Print("v");
  
}



// Generate data for B+ or B- only
RooDataSet* genData(int charge, Double_t scale)
{
  useBothFitVars();
  RooDataSet* d = pdfOnResDK.generate(scale*pdfOnResDK.totalNumEvts());

  // Throw away other charge 
  if (charge<0) d = (RooDataSet*)d->reduce("Hdtrkchge>0");
  else d = (RooDataSet*)d->reduce("Hdtrkchge<0");
}



// set CP parameters to no CP violation
void setNoCP()
{
  pdfOnResDK.rhoMinus()->setVal(dalitzHolderN.sigGoodD0Type().x0()->getVal());
  pdfOnResDK.rhoPlus()->setVal(dalitzHolderP.sigGoodD0Type().x0()->getVal());
  pdfOnResDK.thetaMinus()->setVal(180);
  pdfOnResDK.thetaPlus()->setVal(180);
}

// flip sign of asymmetric efficiency coefficients
void flipEffAsym() {
  for (int i=1; i<5; i++) {
    eff.a(i)->setVal(-1*eff.a(i)->getVal());
    effOther.a(i)->setVal(-1*effOther.a(i)->getVal());
  }
}
