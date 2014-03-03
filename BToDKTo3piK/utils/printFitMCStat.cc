// $Id: printFitMCStat.cc,v 1.3 2006/04/28 01:56:18 fwinkl Exp $
// Script to print stats created by analysis/fitMC.cc

// Print the average initial and final fit values
void printFitMCStat(const char* fitFile, Bool_t printLatex = kFALSE)
{

  // Print some information from the log file
  TString logFile(fitFile);
  logFile.ReplaceAll(".root",".log");
  TString cmd = TString("grep 'todo = ' ")+logFile;
  gSystem->Exec(cmd);
  cmd = TString("grep 'toreplace =' ")+logFile;
  gSystem->Exec(cmd);

  TFile f(fitFile);

  int i = 0;
  while (1) {
    TString name = "initValues";
    name += i;
    RooDataSet* initValues = (RooDataSet*)f.Get(name);
    name = "finalValues";
    name += i;
    RooDataSet* finalValues = (RooDataSet*)f.Get(name);
    
    if (finalValues==0 || initValues==0) return;
  
    if ((initValues->numEntries() != finalValues->numEntries()) ||
        (initValues->get()->getSize() != finalValues->get()->getSize())) {
      cout << "Initial and final fit result datasets are not compatible."<<endl;
      return;
    }
    
    // Calculate mean values
    TIterator* iter = initValues->get()->createIterator();
    RooRealVar* r;
    cout << "--------------------------------------------------------------------"<<endl;
    cout << "Average over "<<finalValues->numEntries()<<" fits:"<<endl;

 
    RooArgSet initSet;
    RooArgSet finalSet;
    while ((r = (RooRealVar*)iter->Next())) {
      Double_t initMean = getMean(*initValues,*r);
      
      Double_t finalMean = getMean(*finalValues,*r);
      Double_t finalMeanErr = getMean((TTree&)finalValues->tree(),
                                      TString(r->GetName())+"_err");
     
      Double_t change = (finalMean-initMean)/finalMeanErr;

      RooRealVar* rtmp = new RooRealVar(r->GetName(),"",initMean);
      initSet.addOwned(*rtmp);
      rtmp = new RooRealVar(r->GetName(),"",finalMean);
      rtmp->setError(finalMeanErr);
      finalSet.addOwned(*rtmp);
      
      if (!printLatex) printf("%-30s%15.4f%15.4f%10.4f%7.1f\n",
                              r->GetName(),initMean,finalMean,finalMeanErr,change);
    }
    if (printLatex)
      initSet.printLatex(Sibling(finalSet),Format("E",AutoPrecision(3)));
    delete initValues;
    delete finalValues;
    i++;
  }
  f.Close();
}





// Print all final fit values for var
void printFitMCStat(const char* fitFile, const char* var)
{
 
  TFile f(fitFile);

  int i = 0;
  while (1) {
    TString name = "finalValues";
    name += i;
    RooDataSet* finalValues = (RooDataSet*)f.Get(name);
    
    if (finalValues==0) return;

    for (int j=0; j<finalValues->numEntries(); j++) {
      RooRealVar* r = (RooRealVar*)finalValues->get(j)->find(var);
      if (r) r->Print();
    }

    delete finalValues;
    i++;
  }
  f.Close();
}
