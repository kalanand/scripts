// calculates the true charge asymmetry in a data set

#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooAbsData.hh"

#include "../BToDKTo3piK/globals/globals.hh"

double trueAsym(RooAbsData * dat) {
  if (0 == dat) {
    return 0;
  }

  int nPos = 0;
  int nNeg = 0;

  for(int i=0; i< dat->numEntries(); ++i) {
    const RooArgSet * event = dat->get(i);
    RooCategory * kCat = (RooCategory *) event->find(Hdtrkchge->GetName());

    if(1 == kCat->getIndex()) {
      ++nPos;
    }
    if(-1 == kCat->getIndex()) {
      ++nNeg;
    }
  }
 
  if (0 == nPos + nNeg) {
    return 0;
  }
 
  return ((double)(nNeg - nPos)) / (nPos + nNeg);
}
  
