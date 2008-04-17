// Replaces var in each event in the data with the pdf-generated values.
// Caller owns the return value.

#ifndef RANDREP_CC
#define RANDREP_CC

#include "TChain.h"
#include "TCut.h"

#include "RooFitCore/RooDataSet.hh"

#include "../BToDKTo3piK/utils/read.cc"


RooDataSet * randRep(RooDataSet * origData, TCut cut) { 
  RooDataSet* OtData = (RooDataSet*)origData->reduce(cut);

  cout<<" Randomly Selected "<< OtData->numEntries() 
      << " out of initial "<< origData->numEntries()
      << " using cut " << cut.GetTitle() << endl; 

  return OtData;
}

//----------------------------
RooDataSet * randRep(TTree* tree, TCut cut, TCut cut1="") {
//cut is  the standard cut
//cut1 is the random selection cut.
  readCut = cut; 
  RooDataSet * Orig = read(tree);
  if("" == cut1) return Orig; 
  RooDataSet * OtData = (RooDataSet*)Orig->reduce(cut1);
  return OtData;
}       
 

//-----------------------------------------------------------
// give a numerical value and a bin (sample) number instead:
RooDataSet * randRep(RooDataSet * origData, double upperCut, int bin = 0) { 
  cout << "randRep: called with upperCut = " << upperCut << ", bin = " << bin << endl;
  if( upperCut>1) {
    cerr << "*** ERROR: upper cut should be less than 1 " << endl;
  }
  if (0 > bin) {
    cerr << "*** ERROR: randRep: bin has to be >= 0. not cutting" << endl;
    return 0;
  }
  else if (0 == bin) {
	TString cutString = "randAdd<=";
        cutString += upperCut;
        TCut cut(cutString.Data());
        return randRep(origData, cut);
  }
  else {
    if (upperCut * (bin + 1) > 1.0) {
      // can't provide bin bins without repeating the sequence:
      cerr << "*** WARNING randRep: cannot provide independent bin " 
	   << bin << " with upperCut " << upperCut
	   << ". repeating sequence with bin 0." << endl;
      return randRep(origData, upperCut, 0);
    }
    else {
      // calculate the cuts:
      double lowerCut = upperCut * bin;
      double upperCut = upperCut * (bin + 1);
      TString cutString = "randAdd<=";
      cutString += upperCut;
      cutString += "&&randAdd>";
      cutString += lowerCut;
      TCut cut(cutString.Data());
      return randRep(origData, cut);
    }
  }  
}    

//-----------------------------------------------------------
// give a numerical value and a bin (sample) number instead:
RooDataSet * randRep(RooDataSet * origData, double upperCut, double lowCut, int bin = 0 ) {
  cout << "randRep: called with upperCut = " << upperCut << ", lowCut = "<< lowCut <<" ,bin = " << bin << endl;

  if( lowCut <0 || upperCut <0 ) {
     cerr<< " *** ERROR, cut should be >0 "<<endl;
     return 0;
  }

  if( lowCut > 1 || upperCut > 1 ) {
     cerr<< " *** ERROR, cut should be <= 1 "<<endl;
     return 0;
  }

  if (0 > bin) {
    cerr << "*** ERROR: randRep: bin has to be >= 0. not cutting" << endl;
    return 0;
  }
  else if (0 == bin) {
       TString cutString = "randAdd<=";
       cutString += upperCut;
       cutString += "&&randAdd>";
       cutString += lowCut;
       TCut cut(cutString.Data());
       return randRep(origData, cut);
  }
  else {
    if (lowCut + (upperCut- lowCut) * (bin + 1) > 1.0) {
      // can't provide bin bins without repeating the sequence:
      cerr << "*** WARNING randRep: cannot provide independent bin "
           << bin << " with upperCut " << upperCut
           << ". repeating sequence with bin 0." << endl;
      return randRep(origData, upperCut, 0);
    }
    else {
      // calculate the cuts:
      double lowerCut = lowCut + (upperCut - lowCut) * bin;
      double upperCut = lowCut + (upperCut - lowCut) * (bin + 1);
      TString cutString = "randAdd<=";
      cutString += upperCut;
      cutString += "&&randAdd>";
      cutString += lowerCut;
      TCut cut(cutString.Data());
      return randRep(origData, cut);
    }
  }
}

#endif
