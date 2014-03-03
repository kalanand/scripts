// $Id: yieldsNllSyst.cc,v 1.4 2006/07/11 21:44:04 fwinkl Exp $
// Determine systematic errors from the BdkOnResNLLYields term
// Run this in batch with:  bbrroot -q -b setup.cc ".x yieldsNllSyst.cc"

void yieldsNllSyst(TTree* tree = dataTree,
                   TCut addCut = "")
{

  TList fitResults;

  // Read data
  readCut = cutSigReg+addCut;
  data = read(tree);

  // reference fit
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting with original parameters"<<endl;
  cout << "-----------------------------------------------------"<<endl;

  readOnResDKPar();
  pdfOnResDK.setNLLYieldsSystBit(-1);  // use all bits
  fit(pdfOnResDK, *data, true);
  RooFitResult* yieldFit = pdfOnResDK.yieldFitResult();

  RooArgSet* result = pdfOnResDK.fitResult();
  result->setName("ref");
  fitResults.Add(result);

  fitResults.Add(doSystBitFit(BdkOnResNLLYields::EFF,"EFF",*yieldFit));
  fitResults.Add(doSystBitFit(BdkOnResNLLYields::BRDK,"BRDK",*yieldFit));
  fitResults.Add(doSystBitFit(BdkOnResNLLYields::BR3PIoverK2PI,"BR3PIoverK2PI",*yieldFit));
  fitResults.Add(doSystBitFit(BdkOnResNLLYields::BRK2PI,"BRK2PI",*yieldFit));
  fitResults.Add(doSystBitFit(BdkOnResNLLYields::NBB,"NBB",*yieldFit));
  fitResults.Add(doSystBitFit(BdkOnResNLLYields::ALL,"ALL",*yieldFit));

  TFile f("yieldsNllSyst.root","recreate");
  fitResults.Write();
  f.Close();
}


RooArgSet* doSystBitFit(Int_t removeBit, TString title, RooFitResult& yieldFit)
{
  cout << "-----------------------------------------------------"<<endl;
  cout << "   Fitting without "<<title<<" systematics"<<endl;
  cout << "-----------------------------------------------------"<<endl;

  readOnResDKPar();
  pdfOnResDK.setNLLYieldsSystBit(BdkOnResNLLYields::ALL - removeBit);
  fit(pdfOnResDK, *data, yieldFit);
  
  RooArgSet* result = pdfOnResDK.fitResult();
  result->setName(title);
  return result;
}


void printYieldsNllSyst(const char* file, const char* fmt = "%.4f")
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

  printYieldsNllSyst(setList, fmt);
}

void printYieldsNllSyst(TList& list, const char* fmt = "%.4f")
{

  // Variables to evaluate systematics for
  RooArgSet syst;
  //  syst.add(*pdfOnResDK.sigGoodD0NumEvts());
  //  syst.add(*pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));
  syst.add(pdfOnResDK.cpParams());

  cout << "Fit results:"<<endl;

  for (int i=0; i<list.GetEntries(); i++) {
    RooArgSet* set = ((RooArgSet*)list.At(i))->selectCommon(syst);
    cout << setw(30) << set->GetName() << " "
         << printOneLine(*set, fmt, true) << endl;
  }

  cout << "Change in errors:"<<endl;

  RooArgSet* ref = (RooArgSet*)list.FindObject("ref");
  removeAllRanges(*ref);
  copyErrorToVal(*ref);
  sqr(*ref);
  for (int i=0; i<list.GetEntries(); i++) {
    RooArgSet* set = (RooArgSet*)list.At(i);
    if (TString(set->GetName())=="ref") continue;
    removeAllRanges(*set);
    copyErrorToVal(*set);
    sqr(*set);
    RooAbsCollection* diff = ref->snapshot(false);
    removeAllRanges(*diff);
    diff->setName(set->GetName());
    sub(*diff, *set);
    sqrt(*diff);
    
    RooArgSet* printSet = diff->selectCommon(syst);
    cout << setw(20) << printSet->GetName() << " "
         << printOneLine(*printSet, fmt, false) << endl;
  }
}
