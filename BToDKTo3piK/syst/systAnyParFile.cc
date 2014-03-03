// $Id: systAnyParFile.cc,v 1.3 2007/03/02 09:49:21 fwinkl Exp $
//
// Obtain the systematic from an arbitray set of .par files
// Run this in batch with:  
//   bbrroot -q -b setup.cc '.x systAnyParFile.cc("../BToDKTo3piK/params/syst/syst1Dshapes.parlist")'


void systAnyParFile(const char* parListFile, TString postfix = "")
{
  RooArgList* parList = readParList(parListFile);
  parList->Print("v");

  // Read data
  readCut = cutSigReg;
  data = read(dataTree);

   // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with original parameters"<<endl;
  cout << "-----------------------------------------------------"<<endl;
  useBothFitVars();
  readOnResDKPar();
  fit(pdfOnResDK, *data, true);
  
  TList fitResults;
  RooArgSet* result = pdfOnResDK.fitResult();
  result->setName("ref");
  fitResults.Add(result);

  TIterator* iter = parList->createIterator();
  RooRealVar* r;
  while (r = (RooRealVar*)iter->Next()) {

    cout << "-----------------------------------------------------"<<endl;
    cout << "   Fitting with "<<r->GetTitle()<<" parameters"<<endl;
    cout << "-----------------------------------------------------"<<endl;

    // reset PDF and read par file  
    useBothFitVars();
    readOnResDKPar();
    cout << "Reading "<<r->GetName()<<endl;
    pdfOnResDK.parameters().readFromFile(r->GetName());

    // re-normalize if requested
    //    BdkDDalitzAmp::setUseFixedWidth(true);  // for fixed BW width systematic
    if (r->getVal()>0) BdkDDalitzAmp::normalizeAll();

    // fit
    fit(pdfOnResDK, *data, true);
    RooArgSet* result = pdfOnResDK.fitResult();
    result->setName(r->GetTitle());
    fitResults.Add(result);
    result->Print("v");
  }

  TFile f(TString("systAnyParFile")+postfix+".root","recreate");
  fitResults.Write();
  f.Close();

  delete parList;
}



// return a list of RooRealVars from a file with this format:
//
// ../BToDKTo3piK/params/syst/dstar_PPP_Group0.Bdk.par    Dalitz0     1
// ../BToDKTo3piK/params/syst/dstar_PPP_Group1.Bdk.par    Dalitz1     0
//
// The first colums becomes the name, the second the title, the third the value.
// value indicates if the PDF needs to be renormalized before fitting.
RooArgList* readParList(const char* parListFile)
{
  RooArgList* list = new RooArgList(parListFile);
  ifstream ifs;
  ifs.open(parListFile);
  while (!ifs.eof()) {
    string path, title;
    int norm;
    ifs >> path >> title >> norm;
    if (path=="") continue;
    RooRealVar* r = new RooRealVar(path.c_str(),title.c_str(),norm);
    list->addOwned(*r);
  }

  return list;
}



void printSystAnyParFile(const char* file, TString fmt = "%.4f")
{
  TFile f(file);
  TList* keys = f.GetListOfKeys();
  TIterator* iter = keys->MakeIterator();
  TKey* k;
  TList setList;
  while (k = (TKey*)iter->Next()) {
    RooArgSet* set = (RooArgSet*)k->ReadObj();
    setList.Add(set);
  }
  delete iter;
  f.Close();

  printSystAnyParFile(setList, fmt);
}



void printSystAnyParFile(TList& list, TString fmt = "%.4f")
{

  // Variables to evaluate systematics for
  RooArgSet syst;
  syst.add(*pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));
  syst.add(*pdfOnResDK.sigGoodD0NumEvts());
  syst.add(pdfOnResDK.cpParams());
  
  RooArgSet* ref = (RooArgSet*)list.FindObject("ref");

  cout << "Fit results:"<<endl;
  for (int i=0; i<list.GetEntries(); i++) {
    RooArgSet* set = (RooArgSet*)list.At(i);  
    cout << setw(20) << set->GetName() << " "
         << printOneLine(*set->selectCommon(syst),fmt) << endl;
  }

  TMatrixDSym V(syst.getSize());
  RooArgList systList(syst);

  cout << "Difference to reference:"<<endl;
  for (int i=0; i<list.GetEntries(); i++) {
    RooArgSet* set = (RooArgSet*)((RooArgSet*)list.At(i))->selectCommon(syst);
    RooAbsCollection* diff = ref->selectCommon(syst)->snapshot(false);
    removeAllRanges(*diff);
    sub(*diff,*set);

    cout << setw(20) << set->GetName() << " "
         << printOneLine(*diff->selectCommon(syst),fmt) << endl;

    for (int n=0; n<V.GetNrows(); n++) {    
      for (int m=0; m<=n; m++) {
        V(n,m) += ((RooArgSet*)diff)->getRealValue(systList[n].GetName()) *
                  ((RooArgSet*)diff)->getRealValue(systList[m].GetName());
        V(m,n) = V(n,m);
      }
    }
    delete diff;
  }
  
  cout << "Covariance matrix:"<<endl;
  V.Print();
}
