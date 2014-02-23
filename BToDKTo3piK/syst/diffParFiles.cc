// $Id: diffParFiles.cc,v 1.1 2006/07/11 21:48:52 fwinkl Exp $
// Script that simply calculates the differenced between two .par files
// Only floating parameters are considered

void diffParFiles(const char* parFile,
                  const char* refFile)
{
  pdfOnResDK.parameters().readFromFile(refFile);  
  RooAbsCollection* ref = pdfOnResDK.parametersFree().snapshot(false);
  removeAllRanges(*ref);

  pdfOnResDK.parameters().readFromFile(parFile);  
  RooAbsCollection* par = pdfOnResDK.parametersFree().snapshot(false);
  removeAllRanges(*par);

  cout << "Reference fit:"<<endl;
  ref->Print("v");

  cout << "Alternative fit:"<<endl;
  par->Print("v");

  if (ref->getSize() != par->getSize()) {
    cout << "WARNING: The number of floating parameters is different between "
         << "the reference and alternative fit. Only showing differences for floating "
         << "parameters in reference fit."<<endl;    
    par = par->selectCommon(*ref);
  }

  cout << "Difference:"<<endl;
  sub(*ref, *par);
  ref->Print("v");
}
