#include <fstream>

void fitToy(const char* resultFile, 
	    UInt_t Seed, 
	    const char * parInputFile = "toy/defaultParInputFile.cc",
	    const char * floatFile = "toy/defaultFloatFile.cc") {
  
  // set seed:
  setRandomGenSeed(Seed);

  // Execute commands before generating events:
  cout << "reading parInputFile " << parInputFile << endl;
  gROOT->Macro(parInputFile);

  // figure out how many events to generate:
  int localNumEvtsToGenerate = numEvtsToGenerate;
  if (0 == localNumEvtsToGenerate) {
    localNumEvtsToGenerate = pdfOnResDK.totalNumEvts();
  }

  //generate MC toy study data
  cout << "Generating " << numEvtsToGenerate << " events" << endl;
  data = pdfOnResDK.generate(localNumEvtsToGenerate);

  // Fix everything, then execute the file that determines what's floating:
  pdfOnResDK.fixAll();
  cout << "executing floatFile " << floatFile << endl;
  gROOT->Macro(floatFile);

  // Store initial (true) values:
  RooArgSet paramsFree = pdfOnResDK.parametersFree();
  double initValues[300];

  TIterator *iter= paramsFree.createIterator();
  RooRealVar *arg ;
  int i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    initValues[i] = arg->getVal();
  }
  
  // Fit:
  cout<<"Fitting the variables"<<endl;
  paramsFree.Print("V");
  RooFitResult * result = fit(pdfOnResDK, *data);

  // Report the results:
  ofstream ofs(resultFile);
  Int_t status = -100;
  Double_t minNll = 100;
  if(0 != result ) {
    status = result->status(); 
    minNll = result->minNll();
  }
  ofs << "fitStatus " << status << " " << minNll 
      << " " << paramsFree.getSize();

  iter= paramsFree.createIterator();
  i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    double finalValue = arg->getVal();
    double error = arg->getError();
    ofs << " " << initValues[i] << " " << finalValue << " " << error << " ";
  }
  ofs << endl;
  cout << "exiting fitToy" << endl;
}



